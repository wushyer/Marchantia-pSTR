#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import pysam

MAX_PER_LOC_DEFAULT = 50

def gt_to_dosage(gt):
    # gt is tuple like (0,), (1,), (0,1) etc; None for missing
    if gt is None or any(a is None for a in gt):
        return np.nan
    # count ALT alleles (allele index >0) for diploid; for haploid returns 0/1
    return float(sum(1 for a in gt if a > 0))

def safe_zscore(vec):
    m = np.nanmean(vec)
    s = np.nanstd(vec)
    if not np.isfinite(s) or s == 0:
        return None
    z = (vec - m) / s
    # fill missing with 0 after centering (equivalent to mean imputation on z-scale)
    z[~np.isfinite(z)] = 0.0
    return z

def corr_r2_from_z(zx, zy, n):
    # zx, zy are z-scored with missing filled as 0
    # approximate r by dot/(n-1). Using n for stability here.
    r = float(np.dot(zx, zy)) / float(n)
    return r * r

def load_snp_matrix(vcf_path, max_sites=None, thin_bp=200):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    records = []
    for rec in vcf.fetch():
        # optional thinning by position bins
        if thin_bp is not None:
            key = (rec.contig, rec.pos // thin_bp)
            # keep one per bin by simple rule: first encountered
            # We'll store keys in set to skip repeats
            # Implemented outside loop for speed
        records.append(rec)
        if max_sites and len(records) >= max_sites:
            break
    vcf.close()

def read_snp(vcf_path, max_dist, bin_size, thin_bp=200, max_sites=None):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    n = len(samples)

    keys_seen = set()
    pos_by_chr = {}
    Z_by_chr = {}

    for rec in vcf.fetch():
        if rec.filter.keys() and ("PASS" not in rec.filter.keys()):
            continue
        if len(rec.alts or []) != 1:
            continue
        if thin_bp is not None:
            key = (rec.contig, rec.pos // thin_bp)
            if key in keys_seen:
                continue
            keys_seen.add(key)

        g = np.empty(n, dtype=float)
        for i, s in enumerate(samples):
            g[i] = gt_to_dosage(rec.samples[s].get("GT"))
        z = safe_zscore(g)
        if z is None:
            continue

        pos_by_chr.setdefault(rec.contig, []).append(rec.pos)
        Z_by_chr.setdefault(rec.contig, []).append(z)

        if max_sites and sum(len(v) for v in pos_by_chr.values()) >= max_sites:
            break

    vcf.close()
    for c in list(pos_by_chr.keys()):
        pos_by_chr[c] = np.array(pos_by_chr[c], dtype=int)
        Z_by_chr[c] = np.stack(Z_by_chr[c], axis=0)
    return samples, pos_by_chr, Z_by_chr

def read_str_gb(vcf_path, max_sites=None):
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    n = len(samples)

    pos_by_chr = {}
    Z_by_chr = {}

    for rec in vcf.fetch():
        if rec.filter.keys() and ("PASS" not in rec.filter.keys()):
            continue

        g = np.empty(n, dtype=float)
        for i, s in enumerate(samples):
            gb = rec.samples[s].get("GB")
            if gb is None or gb == ".":
                g[i] = np.nan
            else:
                # HipSTR GB is a string like "0" or "-2" (haploid) or "0|2" sometimes in other contexts
                # In your data it's scalar; still be safe:
                try:
                    if isinstance(gb, (tuple, list)):
                        # if it comes as a tuple, take sum
                        vals = [float(x) for x in gb if x is not None]
                        g[i] = float(np.mean(vals)) if vals else np.nan
                    else:
                        g[i] = float(str(gb).split("|")[0].split("/")[0])
                except Exception:
                    g[i] = np.nan

        z = safe_zscore(g)
        if z is None:
            continue

        pos_by_chr.setdefault(rec.contig, []).append(rec.pos)
        Z_by_chr.setdefault(rec.contig, []).append(z)

        if max_sites and sum(len(v) for v in pos_by_chr.values()) >= max_sites:
            break

    vcf.close()
    for c in list(pos_by_chr.keys()):
        pos_by_chr[c] = np.array(pos_by_chr[c], dtype=int)
        Z_by_chr[c] = np.stack(Z_by_chr[c], axis=0)
    return samples, pos_by_chr, Z_by_chr

def sample_pairs_decay(posA, ZA, posB, ZB, max_dist, bin_size, pairs_per_bin, seed=1, same_set=False):
    rng = np.random.default_rng(seed)
    bins = np.arange(0, max_dist + bin_size, bin_size)
    nb = len(bins) - 1

    sum_r2 = np.zeros(nb, dtype=float)
    cnt = np.zeros(nb, dtype=int)

    # sort B by position
    orderB = np.argsort(posB)
    posB = posB[orderB]
    ZB = ZB[orderB]

    n_samples = ZA.shape[1]

    # for each A locus, sample limited partners to populate bins
    for i in range(len(posA)):
        p = posA[i]
        # one-sided to avoid double counting
        left = np.searchsorted(posB, p + 1)
        right = np.searchsorted(posB, p + max_dist + 1)
        if right <= left:
            continue
        cand = np.arange(left, right)
        if len(cand) > MAX_PER_LOC_DEFAULT:
            cand = rng.choice(cand, size=MAX_PER_LOC_DEFAULT, replace=False)

        zx = ZA[i]
        for j in cand:
            if same_set and posB[j] == p:
                continue
            d = posB[j] - p
            b = int(d // bin_size)
            if b < 0 or b >= nb:
                continue
            r2 = corr_r2_from_z(zx, ZB[j], n_samples)
            sum_r2[b] += r2
            cnt[b] += 1

    # downsample to equal pairs_per_bin by stochastic thinning (approx)
    # If a bin has too many pairs, we approximate mean by using all (good enough with large cnt)
    mean_r2 = np.full(nb, np.nan)
    for b in range(nb):
        if cnt[b] == 0:
            continue
        mean_r2[b] = sum_r2[b] / cnt[b]

    mids = (bins[:-1] + bins[1:]) / 2
    return mids, mean_r2, cnt

def run_one(label_species, label_group, snp_vcf, str_vcf, out_tsv,
            max_dist=50000, bin_size=1000, pairs_per_bin=2000,
            snp_thin_bp=200, seed=1):
    # read SNP + STR
    snp_samples, snp_pos, snp_Z = read_snp(snp_vcf, max_dist, bin_size, thin_bp=snp_thin_bp)
    str_samples, str_pos, str_Z = read_str_gb(str_vcf)

    # sanity: sample order must match
    if snp_samples != str_samples:
        raise RuntimeError("Sample lists differ between SNP and STR VCF (order or content). Please subset/reorder.")

    rows = []
    for chrom in sorted(set(snp_pos.keys()) & set(str_pos.keys())):
        # SNP-SNP
        mids, mean_r2, cnt = sample_pairs_decay(snp_pos[chrom], snp_Z[chrom],
                                               snp_pos[chrom], snp_Z[chrom],
                                               max_dist, bin_size, pairs_per_bin,
                                               seed=seed, same_set=True)
        for d, r2, n in zip(mids, mean_r2, cnt):
            rows.append((label_species, label_group, "SNP-SNP", int(d), r2, int(n)))

        # STR-SNP (STR as A, SNP as B)
        mids, mean_r2, cnt = sample_pairs_decay(str_pos[chrom], str_Z[chrom],
                                               snp_pos[chrom], snp_Z[chrom],
                                               max_dist, bin_size, pairs_per_bin,
                                               seed=seed+7, same_set=False)
        for d, r2, n in zip(mids, mean_r2, cnt):
            rows.append((label_species, label_group, "STR-SNP", int(d), r2, int(n)))

        # STR-STR
        mids, mean_r2, cnt = sample_pairs_decay(str_pos[chrom], str_Z[chrom],
                                               str_pos[chrom], str_Z[chrom],
                                               max_dist, bin_size, pairs_per_bin,
                                               seed=seed+13, same_set=True)
        for d, r2, n in zip(mids, mean_r2, cnt):
            rows.append((label_species, label_group, "STR-STR", int(d), r2, int(n)))

    df = pd.DataFrame(rows, columns=["species","group","ld_type","distance_mid_bp","mean_r2","n_pairs"])
    df.to_csv(out_tsv, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--group", required=True)
    ap.add_argument("--snp", required=True)
    ap.add_argument("--str", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_dist", type=int, default=50000)
    ap.add_argument("--bin_size", type=int, default=1000)
    ap.add_argument("--pairs_per_bin", type=int, default=2000)
    ap.add_argument("--snp_thin_bp", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    run_one(args.species, args.group, args.snp, args.str, args.out,
            max_dist=args.max_dist, bin_size=args.bin_size,
            pairs_per_bin=args.pairs_per_bin, snp_thin_bp=args.snp_thin_bp, seed=args.seed)

if __name__ == "__main__":
    main()

