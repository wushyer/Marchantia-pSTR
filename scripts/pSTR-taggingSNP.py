#!/usr/bin/env python3
"""
pc_adjusted_pstr_tagging.py

Compute PC-adjusted STR–SNP tagging:
For each STR, scan SNPs within +/- max_dist and take max(r^2) using residual
genotypes after regressing out the top PCs (derived from SNP PCA).

Supports STR input as:
  1) HipSTR VCF (use FORMAT/GB), or
  2) a matrix TSV like: #CHROM  POS  sample1 sample2 ...

Outputs a thresholds TSV:
  group  threshold  linked_pstrs  total_pstrs_tested  fraction_linked  n_samples_used  n_pc  max_dist ...

Example:
  python pc_adjusted_pstr_tagging.py \
    --snp mp.snp.pass.biallelic.vcf.gz \
    --str mp.str.matrix.all.tsv --str_format matrix_tsv \
    --pc at.pca.tsv --n_pc 5 \
    --out mp.pstr_tagging.pc_adjusted.all.tsv \
    --group_name all \
    --thresholds 0.1,0.2,0.3,0.4,0.5 \
    --max_dist 50000 --snp_thin_bp 200 --max_snp_per_str 200 --min_n 10 --seed 1
"""

import argparse
import sys
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict


# -----------------------------
# Utilities
# -----------------------------
def read_sample_list(path: str):
    out = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                out.append(s)
    return out


def gt_to_dosage(gt):
    """Count ALT alleles in GT tuple. Works for diploid/haploid."""
    if gt is None or any(a is None for a in gt):
        return np.nan
    return float(sum(1 for a in gt if a > 0))


def parse_hipstr_gb(gb):
    """
    HipSTR FORMAT/GB: base-pair diff(s) from reference.
    In your Marchantia HipSTR data it's typically a scalar string like "0" or "-2".
    Be tolerant to "0|2" or tuples.
    """
    if gb is None:
        return np.nan
    if isinstance(gb, (tuple, list)):
        vals = []
        for x in gb:
            if x is None:
                continue
            try:
                vals.append(float(str(x).split("|")[0].split("/")[0]))
            except Exception:
                pass
        return float(np.mean(vals)) if vals else np.nan

    s = str(gb).strip()
    if s == "." or s.lower() == "none" or s == "":
        return np.nan
    try:
        s = s.split("|")[0].split("/")[0]
        return float(s)
    except Exception:
        return np.nan


def parse_matrix_cell(x):
    """Parse STR matrix cell (numeric bp diff), treating empty/None as missing."""
    if x is None:
        return np.nan
    s = str(x).strip()
    if s == "" or s.lower() == "none" or s == "." or s == "NA" or s == "nan":
        return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def load_pc_table(pc_path: str, n_pc: int):
    """
    Expect a TSV with columns: sample, PC1..PCk (header required).
    Returns:
      samples (list)
      PC matrix shape (n_samples, n_pc)
    """
    pc = pd.read_csv(pc_path, sep="\t", dtype=str)
    if pc.shape[1] < 2:
        raise RuntimeError("PC file must have at least 2 columns: sample and PC1..")

    sample_col = pc.columns[0]
    want_cols = [f"PC{i}" for i in range(1, n_pc + 1)]
    missing = [c for c in want_cols if c not in pc.columns]
    if missing:
        raise RuntimeError(f"PC file missing columns: {missing}")

    samples = pc[sample_col].astype(str).tolist()
    X = pc[want_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return samples, X


def subset_order(samples_all, keep_samples):
    keep_set = set(keep_samples)
    idx = [i for i, s in enumerate(samples_all) if s in keep_set]
    kept = [samples_all[i] for i in idx]
    return idx, kept


# -----------------------------
# PC residualization + r2
# -----------------------------
def residualize_with_pcs(vec, pc_mat, min_n=10):
    """
    vec: shape (n,), can have NaN
    pc_mat: shape (n, k), can have NaN (rare)
    Returns residual vec_resid (n,) with NaN where vec was NaN.
    Regression: vec ~ intercept + PCs, using only non-missing rows.
    """
    m = np.isfinite(vec)
    if not np.any(m):
        return None
    Xm = pc_mat[m, :]
    ym = vec[m]

    # also require PCs finite
    m2 = np.all(np.isfinite(Xm), axis=1)
    if m2.sum() < min_n:
        return None
    Xm = Xm[m2, :]
    ym = ym[m2]

    # add intercept
    A = np.column_stack([np.ones(Xm.shape[0]), Xm])

    # solve least squares
    try:
        beta, *_ = np.linalg.lstsq(A, ym, rcond=None)
    except Exception:
        return None

    yhat = A @ beta
    resid_m = ym - yhat

    resid = np.full_like(vec, np.nan, dtype=float)
    # put residuals back into the corresponding indices
    idx_m = np.where(m)[0][m2]
    resid[idx_m] = resid_m
    return resid


def pairwise_r2(x, y, min_n=10):
    """
    Pairwise-complete r^2 (squared Pearson correlation) on finite pairs.
    x,y: arrays length n with NaN allowed.
    """
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < min_n:
        return None
    a = x[m]
    b = y[m]
    sa = a.std(ddof=1)
    sb = b.std(ddof=1)
    if not np.isfinite(sa) or not np.isfinite(sb) or sa == 0 or sb == 0:
        return None
    a = (a - a.mean()) / sa
    b = (b - b.mean()) / sb
    r = float(np.dot(a, b) / n)
    return r * r


# -----------------------------
# Readers
# -----------------------------
def read_snp_by_chr(vcf_path, sample_idx, thin_bp=200, max_snps=None):
    """
    Load SNP genotypes per chromosome.
    Keeps PASS + biallelic.
    Returns:
      pos_by_chr: dict chrom -> np.array(pos)
      G_by_chr:   dict chrom -> np.array(shape=(m_sites, n_samples_kept))
    """
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)

    pos_by_chr = defaultdict(list)
    G_by_chr = defaultdict(list)

    seen = set()
    total = 0

    for rec in vcf.fetch():
        # filters
        if rec.filter.keys() and ("PASS" not in rec.filter.keys()):
            continue
        if len(rec.alts or []) != 1:
            continue

        # thinning
        if thin_bp is not None:
            key = (rec.contig, rec.pos // thin_bp)
            if key in seen:
                continue
            seen.add(key)

        g = np.empty(len(sample_idx), dtype=float)
        for k, i in enumerate(sample_idx):
            s = samples[i]
            g[k] = gt_to_dosage(rec.samples[s].get("GT"))

        pos_by_chr[rec.contig].append(rec.pos)
        G_by_chr[rec.contig].append(g)

        total += 1
        if max_snps and total >= max_snps:
            break

    vcf.close()

    for c in list(pos_by_chr.keys()):
        order = np.argsort(pos_by_chr[c])
        pos_by_chr[c] = np.array(pos_by_chr[c], dtype=int)[order]
        G_by_chr[c] = np.stack(G_by_chr[c], axis=0)[order]

    return pos_by_chr, G_by_chr


def read_str_hipstr_vcf_by_chr(vcf_path, sample_idx, max_strs=None):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)

    pos_by_chr = defaultdict(list)
    G_by_chr = defaultdict(list)

    total = 0
    for rec in vcf.fetch():
        if rec.filter.keys() and ("PASS" not in rec.filter.keys()):
            continue

        g = np.empty(len(sample_idx), dtype=float)
        for k, i in enumerate(sample_idx):
            s = samples[i]
            gb = rec.samples[s].get("GB")
            g[k] = parse_hipstr_gb(gb)

        pos_by_chr[rec.contig].append(rec.pos)
        G_by_chr[rec.contig].append(g)

        total += 1
        if max_strs and total >= max_strs:
            break

    vcf.close()

    for c in list(pos_by_chr.keys()):
        order = np.argsort(pos_by_chr[c])
        pos_by_chr[c] = np.array(pos_by_chr[c], dtype=int)[order]
        G_by_chr[c] = np.stack(G_by_chr[c], axis=0)[order]

    return pos_by_chr, G_by_chr


def read_str_matrix_tsv_by_chr(matrix_tsv_path, keep_samples=None):
    """
    Matrix TSV:
      #CHROM  POS  sample1 sample2 ...
      chr1    1064  ... values ...

    Returns:
      samples(list), pos_by_chr, G_by_chr
    """
    with open(matrix_tsv_path) as f:
        header = f.readline().rstrip("\n").split("\t")

    if len(header) < 3:
        raise RuntimeError("STR matrix header must have at least 3 columns: #CHROM, POS, samples...")

    sample_names = header[2:]
    # subset columns if requested
    if keep_samples is not None:
        keep_set = set(keep_samples)
        col_idx = [i for i, s in enumerate(sample_names) if s in keep_set]
        kept_samples = [sample_names[i] for i in col_idx]
    else:
        col_idx = list(range(len(sample_names)))
        kept_samples = sample_names[:]

    # load file (strings) then parse numeric
    df = pd.read_csv(matrix_tsv_path, sep="\t", dtype=str)
    # rename potential first column "#CHROM" to "CHROM" for convenience
    if df.columns[0].startswith("#"):
        df = df.rename(columns={df.columns[0]: "CHROM"})
    else:
        df = df.rename(columns={df.columns[0]: "CHROM"})
    df = df.rename(columns={df.columns[1]: "POS"})

    pos_by_chr = defaultdict(list)
    G_by_chr = defaultdict(list)

    # iterate rows (ok for ~30k)
    for _, row in df.iterrows():
        chrom = str(row["CHROM"])
        pos = int(float(row["POS"]))  # safe for string POS
        # extract sample columns
        vals = []
        for i in col_idx:
            sname = sample_names[i]
            vals.append(parse_matrix_cell(row.get(sname, np.nan)))
        g = np.array(vals, dtype=float)

        pos_by_chr[chrom].append(pos)
        G_by_chr[chrom].append(g)

    # sort per chr
    for c in list(pos_by_chr.keys()):
        order = np.argsort(pos_by_chr[c])
        pos_by_chr[c] = np.array(pos_by_chr[c], dtype=int)[order]
        G_by_chr[c] = np.stack(G_by_chr[c], axis=0)[order]

    return kept_samples, pos_by_chr, G_by_chr


# -----------------------------
# Tagging (PC-adjusted)
# -----------------------------
def compute_max_r2_per_str_pc_adjusted(
    str_pos, str_G,
    snp_pos, snp_G,
    pc_mat,
    max_dist=50000,
    max_snp_per_str=200,
    min_n=10,
    seed=1
):
    """
    For each STR, compute max r^2 with SNPs within +/- max_dist, using PC-adjusted residuals.
    Returns max_r2 array length n_str (NaN if no valid comparisons).
    """
    rng = np.random.default_rng(seed)
    n_str = str_G.shape[0]
    max_r2 = np.full(n_str, np.nan, dtype=float)

    # cache SNP residuals (optional but saves a lot when many STRs hit same SNPs)
    snp_resid_cache = {}

    for i, p in enumerate(str_pos):
        left = np.searchsorted(snp_pos, p - max_dist, side="left")
        right = np.searchsorted(snp_pos, p + max_dist, side="right")
        if right <= left:
            continue

        idx = np.arange(left, right)
        if len(idx) > max_snp_per_str:
            idx = rng.choice(idx, size=max_snp_per_str, replace=False)

        x = str_G[i]
        x_resid = residualize_with_pcs(x, pc_mat, min_n=min_n)
        if x_resid is None:
            continue

        best = None
        for j in idx:
            # SNP residualize with cache
            if j in snp_resid_cache:
                y_resid = snp_resid_cache[j]
            else:
                y = snp_G[j]
                y_resid = residualize_with_pcs(y, pc_mat, min_n=min_n)
                snp_resid_cache[j] = y_resid
            if y_resid is None:
                continue

            r2 = pairwise_r2(x_resid, y_resid, min_n=min_n)
            if r2 is None:
                continue
            if best is None or r2 > best:
                best = r2

        if best is not None:
            max_r2[i] = best

    return max_r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--snp", required=True, help="SNP VCF.gz (PASS, biallelic recommended)")
    ap.add_argument("--str", required=True, help="STR input: HipSTR VCF.gz or matrix TSV")
    ap.add_argument("--str_format", required=True, choices=["hipstr_vcf", "matrix_tsv"])
    ap.add_argument("--pc", required=True, help="PC table TSV: sample + PC1..")
    ap.add_argument("--n_pc", type=int, default=5, help="Number of PCs to regress out")
    ap.add_argument("--out", required=True, help="Output TSV for thresholds summary")
    ap.add_argument("--out_maxr2", default=None, help="Optional: output per-STR max_r2 (TSV)")
    ap.add_argument("--group_name", default="all")
    ap.add_argument("--keep_samples", default=None, help="Optional: file of sample IDs to keep (one per line)")
    ap.add_argument("--thresholds", default="0.1,0.2,0.3,0.4")
    ap.add_argument("--max_dist", type=int, default=50000)
    ap.add_argument("--snp_thin_bp", type=int, default=200)
    ap.add_argument("--max_snp_per_str", type=int, default=200)
    ap.add_argument("--min_n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max_snps", type=int, default=None, help="Optional cap for debugging")
    ap.add_argument("--max_strs", type=int, default=None, help="Optional cap for debugging (HipSTR VCF only)")
    args = ap.parse_args()

    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]

    # Load PCs
    pc_samples, pc_mat_all = load_pc_table(args.pc, args.n_pc)

    # Optional keep list
    if args.keep_samples:
        keep = read_sample_list(args.keep_samples)
        keep_set = set(keep)
    else:
        keep = None
        keep_set = None

    # Determine sample order + intersection across SNP + STR + PC
    # SNP samples from VCF header
    snp_vcf = pysam.VariantFile(args.snp)
    snp_samples_all = list(snp_vcf.header.samples)
    snp_vcf.close()

    # STR samples
    if args.str_format == "hipstr_vcf":
        str_vcf = pysam.VariantFile(args.str)
        str_samples_all = list(str_vcf.header.samples)
        str_vcf.close()
        if snp_samples_all != str_samples_all:
            raise RuntimeError("SNP and HipSTR VCF sample lists differ. Subset/reorder first.")
        str_samples_all = snp_samples_all
    else:
        # matrix TSV sample names will be read later; start with SNP list
        str_samples_all = None

    # Compute common samples in desired order (we use SNP order as canonical)
    pc_set = set(pc_samples)

    if args.str_format == "matrix_tsv":
        # Read matrix header to get samples, then intersect with SNP + PCs
        with open(args.str) as f:
            header = f.readline().rstrip("\n").split("\t")
        matrix_samples = header[2:]
        matrix_set = set(matrix_samples)

        common = [s for s in snp_samples_all if (s in matrix_set and s in pc_set)]
        if keep_set is not None:
            common = [s for s in common if s in keep_set]

        if len(common) == 0:
            raise RuntimeError("No overlapping samples among SNP VCF, STR matrix, and PC file (after keep_samples).")

        # SNP indices in VCF order
        snp_idx, common_order = subset_order(snp_samples_all, common)
        # STR matrix will be loaded already subset to these common_order later
        used_samples = common_order

    else:
        # hipstr_vcf: SNP==STR order; need also in PC and keep
        common = [s for s in snp_samples_all if s in pc_set]
        if keep_set is not None:
            common = [s for s in common if s in keep_set]
        if len(common) == 0:
            raise RuntimeError("No overlapping samples between VCF and PC file (after keep_samples).")

        snp_idx, used_samples = subset_order(snp_samples_all, common)

    # subset PC matrix to used_samples in same order
    pc_index = {s: i for i, s in enumerate(pc_samples)}
    pc_mat = np.vstack([pc_mat_all[pc_index[s], :] for s in used_samples]).astype(float)

    # Load SNPs (subset + thin)
    snp_pos_by_chr, snp_G_by_chr = read_snp_by_chr(
        args.snp, snp_idx, thin_bp=args.snp_thin_bp, max_snps=args.max_snps
    )

    # Load STRs
    if args.str_format == "hipstr_vcf":
        str_pos_by_chr, str_G_by_chr = read_str_hipstr_vcf_by_chr(
            args.str, snp_idx, max_strs=args.max_strs
        )
        str_samples_used = used_samples
    else:
        # matrix loader will subset to used_samples and preserve that order
        str_samples_used, str_pos_by_chr, str_G_by_chr = read_str_matrix_tsv_by_chr(
            args.str, keep_samples=used_samples
        )
        if str_samples_used != used_samples:
            raise RuntimeError("STR matrix sample order mismatch after subsetting. Please check headers/IDs.")

    # Compute max r2 per STR across overlapping chromosomes
    max_r2_all = []
    max_r2_records = []  # optional per-STR output

    overlap_chroms = sorted(set(str_pos_by_chr.keys()) & set(snp_pos_by_chr.keys()))
    if len(overlap_chroms) == 0:
        raise RuntimeError("No overlapping chromosomes between STR and SNP inputs.")

    for chrom in overlap_chroms:
        mr2 = compute_max_r2_per_str_pc_adjusted(
            str_pos_by_chr[chrom], str_G_by_chr[chrom],
            snp_pos_by_chr[chrom], snp_G_by_chr[chrom],
            pc_mat,
            max_dist=args.max_dist,
            max_snp_per_str=args.max_snp_per_str,
            min_n=args.min_n,
            seed=args.seed
        )
        max_r2_all.append(mr2)

        if args.out_maxr2 is not None:
            for pos, val in zip(str_pos_by_chr[chrom], mr2):
                max_r2_records.append((chrom, int(pos), val))

    max_r2_all = np.concatenate(max_r2_all)
    tested_mask = np.isfinite(max_r2_all)
    total_tested = int(tested_mask.sum())

    rows = []
    for t in thresholds:
        linked = int(np.sum(max_r2_all[tested_mask] >= t))
        frac = linked / total_tested if total_tested > 0 else np.nan
        rows.append({
            "group": args.group_name,
            "threshold": t,
            "linked_pstrs": linked,
            "total_pstrs_tested": total_tested,
            "fraction_linked": frac,
            "n_samples_used": len(used_samples),
            "n_pc": args.n_pc,
            "max_dist": args.max_dist,
            "max_snp_per_str": args.max_snp_per_str,
            "snp_thin_bp": args.snp_thin_bp,
            "min_n": args.min_n
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, sep="\t", index=False)

    if args.out_maxr2 is not None:
        mr2_df = pd.DataFrame(max_r2_records, columns=["CHROM", "POS", "max_r2_pc_adjusted"])
        mr2_df.to_csv(args.out_maxr2, sep="\t", index=False)

    print(f"[OK] Wrote thresholds: {args.out}", file=sys.stderr)
    if args.out_maxr2 is not None:
        print(f"[OK] Wrote per-STR max r2: {args.out_maxr2}", file=sys.stderr)


if __name__ == "__main__":
    main()

