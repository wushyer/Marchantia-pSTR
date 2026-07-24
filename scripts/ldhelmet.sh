#!/usr/bin/env bash
set -euo pipefail

VCF="test.recode.vcf.gz"
REF="ref.fa"
CHR_LIST="chromosomes.txt"
OUTDIR="ldhelmet_test"
THREADS=8

WINDOW=50
THETA=0.01

BURNIN_ITER=1000000
MCMC_ITER=10000000
BLOCK_PENALTY=50

R_MIN=0
R_MAX=100
R_N=101

mkdir -p "${OUTDIR}"/{logs,fa,conf,post,txt,tmp}

echo "[`date`] Preparing sample list"
SAMPLES=$(bcftools query -l "${VCF}")
NSAMPLES=$(echo ${SAMPLES} | wc -w)
echo "[`date`] Samples: ${NSAMPLES}"

FIRST_CHR=$(head -n 1 "${CHR_LIST}")
echo "[`date`] Using ${FIRST_CHR} to generate lookup configuration"

FCHR_VCF_RAW="${OUTDIR}/tmp/${FIRST_CHR}.raw.vcf.gz"
FCHR_VCF="${OUTDIR}/tmp/${FIRST_CHR}.filtered.vcf.gz"
FCHR_FA_DIR="${OUTDIR}/fa/${FIRST_CHR}"
FCHR_MSA="${OUTDIR}/fa/${FIRST_CHR}.fa"
LOOKUP_CONF="${OUTDIR}/conf/lookup.conf"

mkdir -p "${FCHR_FA_DIR}"

if [[ ! -f "${FCHR_VCF}" ]]; then
  bcftools view -r "${FIRST_CHR}" "${VCF}" -Oz -o "${FCHR_VCF_RAW}"
  bcftools index "${FCHR_VCF_RAW}"

  bcftools view \
    -i 'QUAL>30' \
    -m2 -M2 -v snps \
    "${FCHR_VCF_RAW}" \
    -Oz -o "${FCHR_VCF}"

  bcftools index "${FCHR_VCF}"
fi

if [[ ! -f "${FCHR_MSA}" ]]; then
  for S in ${SAMPLES}; do
    OUT_FA="${FCHR_FA_DIR}/${FIRST_CHR}_${S}.fa"
    if [[ ! -f "${OUT_FA}" ]]; then
      samtools faidx "${REF}" "${FIRST_CHR}" \
      | bcftools consensus \
          -s "${S}" \
          -f - \
          "${FCHR_VCF}" \
      | sed "s/^>.*/>${S}/" > "${OUT_FA}"
    fi
  done
  cat "${FCHR_FA_DIR}"/*.fa > "${FCHR_MSA}"
fi

if [[ ! -f "${LOOKUP_CONF}" ]]; then
  ldhelmet find_confs \
    --num_threads ${THREADS} \
    -w ${WINDOW} \
    -o "${LOOKUP_CONF}" \
    "${FCHR_MSA}"
fi

LK_TABLE="${OUTDIR}/lookup.table"
PADE_TABLE="${OUTDIR}/pade.table"

if [[ ! -f "${LK_TABLE}" ]]; then
  ldhelmet table_gen \
    --num_threads ${THREADS} \
    -c "${LOOKUP_CONF}" \
    -t ${THETA} \
    -r ${R_MIN} ${R_MAX} ${R_N} \
    -o "${LK_TABLE}"
fi

if [[ ! -f "${PADE_TABLE}" ]]; then
  ldhelmet pade \
    --num_threads ${THREADS} \
    -c "${LOOKUP_CONF}" \
    -t ${THETA} \
    -x 11 \
    -o "${PADE_TABLE}"
fi

echo "[`date`] Running rjmcmc per chromosome"

while read -r CHR; do
  [[ -z "${CHR}" ]] && continue
  echo "[`date`] Processing ${CHR}"

  CHR_VCF_RAW="${OUTDIR}/tmp/${CHR}.raw.vcf.gz"
  CHR_VCF="${OUTDIR}/tmp/${CHR}.filtered.vcf.gz"
  CHR_FA_DIR="${OUTDIR}/fa/${CHR}"
  CHR_MSA="${OUTDIR}/fa/${CHR}.fa"
  CHR_POST="${OUTDIR}/post/${CHR}.post"
  CHR_RHO="${OUTDIR}/txt/${CHR}.rho.txt"

  mkdir -p "${CHR_FA_DIR}"

  if [[ ! -f "${CHR_VCF}" ]]; then
    bcftools view -r "${CHR}" "${VCF}" -Oz -o "${CHR_VCF_RAW}"
    bcftools index "${CHR_VCF_RAW}"

    bcftools view \
      -i 'QUAL>30' \
      -m2 -M2 -v snps \
      "${CHR_VCF_RAW}" \
      -Oz -o "${CHR_VCF}"

    bcftools index "${CHR_VCF}"
  fi

  if [[ ! -f "${CHR_MSA}" ]]; then
    for S in ${SAMPLES}; do
      OUT_FA="${CHR_FA_DIR}/${CHR}_${S}.fa"
      if [[ ! -f "${OUT_FA}" ]]; then
        samtools faidx "${REF}" "${CHR}" \
        | bcftools consensus \
            -s "${S}" \
            -f - \
            "${CHR_VCF}" \
        | sed "s/^>.*/>${S}/" > "${OUT_FA}"
      fi
    done
    cat "${CHR_FA_DIR}"/*.fa > "${CHR_MSA}"
  fi

  if [[ ! -f "${CHR_POST}" ]]; then
    ldhelmet rjmcmc \
      --num_threads ${THREADS} \
      -w ${WINDOW} \
      -l "${LK_TABLE}" \
      -p "${PADE_TABLE}" \
      --burn_in ${BURNIN_ITER} \
      -n ${MCMC_ITER} \
      -b ${BLOCK_PENALTY} \
      -s "${CHR_MSA}" \
      -o "${CHR_POST}"
  fi

  if [[ ! -f "${CHR_RHO}" ]]; then
    ldhelmet post_to_text \
      -o "${CHR_RHO}" \
      "${CHR_POST}"
  fi

  echo "[`date`] Finished ${CHR}"

done < "${CHR_LIST}"

echo "[`date`] LDhelmet analysis completed"
