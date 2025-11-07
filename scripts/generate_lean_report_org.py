#!/usr/bin/env python3
"""
LEAN Excel Report Generator
Builds a comprehensive variant analysis report with QC metrics and coverage gaps.
"""
import sys
import os
import re
import argparse
from pathlib import Path
from typing import Dict, Optional, Set, Tuple, List
from urllib.parse import quote_plus

import pandas as pd
from cyvcf2 import VCF


# =============================================================================
# CONFIGURATION & CONSTANTS
# =============================================================================

ACMG_SHEET_COLS = [
    "Gene", "Variant", "HGVSc", "HGVSp", "MANE_ID", "Zygosity", "GT",
    "AD_Ref", "AD_Alt", "DP", "GQ", "QUAL", "Consequence", "Exon", "Intron",
    "ClinVar", "ClinVar_ReviewStatus", "ClinVar_Stars", "ClinVar_StarsGlyph",
    "ClinVar_Link", "gnomAD_AF", "REVEL", "SpliceAI_DS_max", "SpliceAI_Event",
    "BayesDel_score", "AM_Pathogenicity", "AM_Class", "HGVS_full"
]

PASS_VARIANT_COLS = [
    "Variant", "Gene", "HGVSc", "HGVSp", "MANE_ID", "Transcript", "Consequence",
    "GT", "Zygosity", "AD_Ref", "AD_Alt", "DP", "GQ", "QUAL", "FILTER", "VAF",
    "gnomAD_AF", "ClinVar", "ClinVar_ReviewStatus", "ClinVar_Stars",
    "ClinVar_StarsGlyph", "REVEL", "SpliceAI_DS_max", "SpliceAI_Event",
    "BayesDel_score", "AM_Pathogenicity", "AM_Class", "HGVS_full"
]

COVERAGE_GAP_COLS = [
    "Gene", "MANE_ID", "Exon", "Chrom", "ExonStart", "ExonEnd", "ExonLen",
    "Pct>=20x", "Pct>=30x", "Gaps<20x_n", "Gaps<20x_bp", "Gaps<20x_intervals",
    "Gaps<30x_n", "Gaps<30x_bp", "Gaps<30x_intervals"
]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def safe_float(x) -> Optional[float]:
    """Safely convert to float, return None on failure."""
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def safe_int(x) -> Optional[int]:
    """Safely convert to int, return None on failure."""
    try:
        return int(x)
    except (ValueError, TypeError):
        return None


def sample_from_vcf(vcf_path: str) -> str:
    """Extract sample ID from VCF, fallback to filename."""
    try:
        v = VCF(vcf_path)
        if v.samples:
            return v.samples[0]
    except Exception:
        pass
    return Path(vcf_path).stem


def load_gene_list(path: Optional[str]) -> Set[str]:
    """Load gene symbols from a file (one per line)."""
    genes = set()
    if path and Path(path).exists():
        with open(path) as f:
            for line in f:
                gene = line.strip().split()[0]
                if gene:
                    genes.add(gene)
    return genes


# =============================================================================
# QC FILE PARSERS
# =============================================================================

def parse_flagstat(path: Optional[str]) -> Dict[str, any]:
    """Parse samtools flagstat output."""
    out = {}
    if not path or not Path(path).exists():
        return out
    
    with open(path) as f:
        for line in f:
            s = line.strip()
            if " in total " in s and "QC-passed" in s:
                n = s.split()[0]
                if n.isdigit():
                    out["Total_Reads"] = int(n)
            
            elif re.search(r"\bmapped\s*\(", s):
                n = s.split()[0]
                out["Mapped_Reads"] = safe_int(n)
                m = re.search(r"\(([\d\.]+)%", s)
                if m:
                    out["Mapped_Percent"] = safe_float(m.group(1))
            
            elif re.search(r"\bduplicates\b", s) and "(" in s:
                n = s.split()[0]
                out["Duplicate_Reads"] = safe_int(n)
                m = re.search(r"\(([\d\.]+)%", s)
                if m:
                    out["Duplicate_Percent"] = safe_float(m.group(1))
    
    return out


def parse_samtools_stats(path: Optional[str]) -> Dict[str, any]:
    """Parse samtools stats for Ti/Tv, Het/Hom, insert size."""
    out = {}
    if not path or not Path(path).exists():
        return out
    
    het = homalt = mean_insert_size = 0
    
    with open(path) as f:
        for line in f:
            if not line.startswith("SN"):
                continue
            
            if "number of heterozygous sites" in line:
                het = int(line.strip().split()[-1])
            elif "number of homozygous-ALT sites" in line:
                homalt = int(line.strip().split()[-1])
            elif "TSTV" in line or "ts/tv" in line:
                try:
                    out["TiTv"] = float(line.strip().split()[-1])
                except (ValueError, IndexError):
                    pass
            elif "insert size average" in line:
                mean_insert_size = float(line.strip().split()[-1])
    
    if het and homalt:
        out["Het_Hom_Ratio"] = round(het / max(homalt, 1), 3)
    
    if mean_insert_size:
        out["Mean_Insert_Size"] = mean_insert_size
    
    return out


def parse_mosdepth_summary_regions(path: Optional[str]) -> Dict[str, float]:
    """Parse mosdepth summary for mean coverage."""
    out = {}
    if not path or not Path(path).exists():
        return out
    
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return out
        
        # Focus on region rows
        region_df = df[df["chrom"].astype(str).str.endswith("_region")].copy()
        if region_df.empty:
            region_df = df.copy()
        
        # Weighted mean coverage
        cov_mean = (region_df["mean"] * region_df["length"]).sum() / region_df["length"].sum()
        out["Mean_Coverage"] = round(float(cov_mean), 2)
        
    except Exception as e:
        print(f"[WARN] Failed to parse mosdepth summary {path}: {e}", file=sys.stderr)
    
    return out


def parse_sexcheck(path: Optional[str]) -> Dict[str, str]:
    """Parse sex check output."""
    if not path or not Path(path).exists():
        return {}
    
    txt = open(path).read().strip()
    m = re.search(r"sex\s*[:=]\s*(\w+)", txt, re.IGNORECASE)
    return {"Sex_Check": m.group(1)} if m else {"Sex_Check": txt}


def acmg_pct_regions_covered(th_path: Optional[str], min_depth: int = 20) -> Dict[str, any]:
    """Compute % of ACMG regions fully covered at >= min_depth."""
    out = {}
    if not th_path or not Path(th_path).exists():
        return out
    
    try:
        df = pd.read_csv(th_path, sep="\t", compression="infer", header=0)
        df.columns = [str(c).lstrip("#") for c in df.columns]
        cols = list(df.columns)
        
        def pick_col(names, default=None):
            for n in cols:
                if str(n).lower() in names:
                    return n
            return default
        
        chrom = pick_col({"chrom", "chr"}, cols[0] if cols else None)
        start = pick_col({"start"}, cols[1] if len(cols) > 1 else None)
        end = pick_col({"end"}, cols[2] if len(cols) > 2 else None)
        region_col = pick_col({"region", "target", "name"}, None)
        
        if not all([chrom, start, end]):
            return out
        
        # Find threshold column
        preferred = f"{int(min_depth)}X"
        thr_col = preferred if preferred in df.columns else None
        if not thr_col:
            for c in cols:
                digits = "".join(ch for ch in str(c) if ch.isdigit())
                if digits and int(digits) == int(min_depth):
                    thr_col = c
                    break
        
        if not thr_col:
            return out
        
        # Calculate coverage
        s = pd.to_numeric(df[start], errors="coerce")
        e = pd.to_numeric(df[end], errors="coerce")
        thr = pd.to_numeric(df[thr_col], errors="coerce")
        df = df.assign(_len=(e - s), _thr=thr).dropna(subset=["_len"]).copy()
        
        if region_col and region_col in df.columns:
            grp = df.groupby(region_col, as_index=False).agg(
                len_sum=("_len", "sum"), thr_sum=("_thr", "sum")
            )
        else:
            df["_id"] = (df[chrom].astype(str) + ":" + 
                        s.astype("Int64").astype(str) + "-" + 
                        e.astype("Int64").astype(str))
            grp = df.groupby("_id", as_index=False).agg(
                len_sum=("_len", "sum"), thr_sum=("_thr", "sum")
            )
        
        total_regions = len(grp)
        if total_regions == 0:
            return out
        
        covered = int((grp["thr_sum"] >= grp["len_sum"]).sum())
        pct = round(covered / total_regions * 100.0, 2)
        
        out["ACMG_PctRegionsCovered"] = pct
        out["ACMG_RegionsCovered"] = covered
        out["ACMG_RegionsTotal"] = total_regions
        
    except Exception as e:
        print(f"[WARN] Failed to parse ACMG thresholds {th_path}: {e}", file=sys.stderr)
    
    return out


# =============================================================================
# COVERAGE DATA LOADERS
# =============================================================================

class CoverageData:
    """Handles loading and querying coverage data."""
    
    def __init__(self, coverage_path: str, r1r2_path: str, fr_path: str):
        self.coverage = self._load_coverage(coverage_path)
        self.r1r2 = self._load_r1r2(r1r2_path)
        self.fr = self._load_fr(fr_path)
    
    def _load_coverage(self, path: str) -> Dict[str, Dict[str, str]]:
        """Load exon coverage summary."""
        data = {}
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                region = parts[0]
                cov = {}
                for p in parts[1:]:
                    key, val = p.split(":")
                    cov[key.strip()] = val.strip().replace("%", "")
                data[region] = cov
        return data
    
    def _load_r1r2(self, path: str) -> Dict[Tuple, Tuple]:
        """Load R1/R2 ratio data."""
        data = {}
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                chrom, start, end, r1, r2, gene, ratio = parts
                data[(chrom, int(start), int(end))] = (r1, r2, gene, ratio)
        return data
    
    def _load_fr(self, path: str) -> Dict[Tuple, Tuple]:
        """Load forward/reverse ratio data."""
        data = {}
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                chrom, start, end, gene, fwd, rev, ffrac, rfrac = parts
                data[(chrom, int(start), int(end))] = (gene, fwd, rev, ffrac, rfrac)
        return data
    
    def get_coverage(self, chrom: str, pos: int) -> Dict[str, str]:
        """Get coverage metrics for a position."""
        for region, cov in self.coverage.items():
            rchrom, coords = region.split(":")
            start, end = map(int, coords.split("-"))
            if str(chrom) == rchrom and start <= pos <= end:
                return cov
        return {">=20x": "NA", ">=30x": "NA", ">=50x": "NA", ">=100x": "NA"}
    
    def get_r1r2(self, chrom: str, pos: int) -> Tuple:
        """Get R1/R2 metrics for a position."""
        for (c, s, e), vals in self.r1r2.items():
            if str(chrom) == c and s <= pos <= e:
                return vals
        return ("NA", "NA", "NA", "NA")
    
    def get_fr(self, chrom: str, pos: int) -> Tuple:
        """Get forward/reverse metrics for a position."""
        for (c, s, e), vals in self.fr.items():
            if str(chrom) == c and s <= pos <= e:
                return vals
        return ("NA", "NA", "NA", "NA", "NA")


# =============================================================================
# CLINVAR UTILITIES
# =============================================================================

def clinvar_stars_from_revstat(revstat: str) -> Optional[int]:
    """Map ClinVar review status to 0-4 stars."""
    if not revstat:
        return None
    
    status_map = {
        "practice guideline": 4,
        "reviewed by expert panel": 3,
        "criteria provided, multiple submitters, no conflicts": 2,
        "criteria provided, conflicting": 1,
        "criteria provided, single submitter": 1,
        "no assertion criteria provided": 0,
        "no classification provided": 0,
        "no classification for the individual variant": 0,
    }
    
    def norm(s: str) -> str:
        return s.lower().replace("_", " ").strip()
    
    toks = re.split(r"[|,;]+", revstat)
    scores = []
    for t in toks:
        t = norm(t)
        for key, score in status_map.items():
            if t.startswith(key):
                scores.append(score)
                break
        else:
            scores.append(0)
    
    return max(scores) if scores else 0


def clinvar_star_glyph(n: Optional[int]) -> str:
    """Convert star count to glyph."""
    if isinstance(n, (int, float)) and n >= 0:
        return "★" * int(n)
    return ""


def clinvar_url_from_row(row: pd.Series) -> str:
    """Build ClinVar URL from variant row."""
    def first_token(x):
        if x is None:
            return ""
        return re.split(r"[|,;]", str(x))[0].strip()
    
    alleleid = first_token(row.get("ClinVar_ALLELEID") or row.get("ALLELEID"))
    if alleleid:
        m = re.search(r"\d+", alleleid)
        if m:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{int(m.group(0))}/"
    
    return ""


# =============================================================================
# VCF PARSING
# =============================================================================

def select_csq_entry(var, csq_format: List[str], prefer_clinvar: bool = True) -> Optional[Dict]:
    """Select best CSQ annotation for variant."""
    csq = var.INFO.get("CSQ")
    if not csq:
        return None
    
    alts = [str(a) for a in (var.ALT or [])]
    alt0 = alts[0] if alts else None
    
    rows = []
    for entry in csq.split(","):
        cols = entry.split("|")
        d = dict(zip(csq_format, cols + [""] * (len(csq_format) - len(cols))))
        rows.append(d)
    
    # Filter by allele match
    if alt0:
        rows_alt = [d for d in rows if d.get("Allele") == alt0]
    else:
        rows_alt = rows[:]
    
    # Prefer ClinVar entries
    if prefer_clinvar:
        rows_cv = [d for d in rows_alt if d.get("ClinVar_ALLELEID") or d.get("ALLELEID")]
        if rows_cv:
            return rows_cv[0]
    
    return rows_alt[0] if rows_alt else (rows[0] if rows else None)


def parse_vcf_to_dataframe(vcf_path: str, coverage_data: CoverageData, 
                          sample_id: str) -> pd.DataFrame:
    """Parse VCF file into DataFrame with all annotations."""
    vcf = VCF(vcf_path)
    
    # Extract CSQ format
    csq_format = []
    for h in vcf.raw_header.split("\n"):
        if h.startswith("##INFO=<ID=CSQ"):
            try:
                csq_format = h.split("Format: ")[1].rstrip('">').split("|")
            except IndexError:
                pass
            break
    
    records = []
    for var in vcf:
        record = parse_variant(var, csq_format, coverage_data, sample_id)
        records.append(record)
    
    return pd.DataFrame(records)


def parse_variant(var, csq_format: List[str], coverage_data: CoverageData, 
                 sample_id: str) -> Dict:
    """Parse a single variant into a record dict."""
    chrom, pos, ref = var.CHROM, var.POS, var.REF
    alt = ",".join(var.ALT) if var.ALT else ""
    filt = var.FILTER or "PASS"
    qual = safe_float(getattr(var, "QUAL", None))
    
    # Genotype info
    dp = var.format("DP")[0][0] if "DP" in var.FORMAT else var.INFO.get("DP")
    vaf = var.format("VAF")[0][0] if "VAF" in var.FORMAT else None
    gq = var.format("GQ")[0][0] if "GQ" in var.FORMAT else None
    
    gt = var.genotypes[0][:2]
    zyg = "het" if gt[0] != gt[1] else ("hom" if gt[0] == 1 else "ref")
    
    ad_ref = ad_alt = None
    if "AD" in var.FORMAT:
        try:
            ad_ref, ad_alt = var.format("AD")[0]
        except (IndexError, ValueError):
            pass
    
    # Annotations
    ann = select_csq_entry(var, csq_format) if csq_format else {}
    
    # Coverage metrics
    cov = coverage_data.get_coverage(chrom, pos)
    r1r2 = coverage_data.get_r1r2(chrom, pos)
    fr = coverage_data.get_fr(chrom, pos)
    
    # Calculate ratios
    try:
        r1, r2 = int(r1r2[0]), int(r1r2[1])
        r1r2_abs = round(r1 / r2, 2) if r2 > 0 else None
    except (ValueError, IndexError):
        r1r2_abs = None
    
    try:
        fwd, rev = int(fr[1]), int(fr[2])
        fr_abs = round(fwd / rev, 2) if rev > 0 else None
    except (ValueError, IndexError):
        fr_abs = None
    
    # ClinVar
    revstat = ann.get("ClinVar_CLNREVSTAT") or ann.get("CLNREVSTAT") or ""
    stars = clinvar_stars_from_revstat(revstat)
    
    # SpliceAI
    spliceai_fields = {
        "DS_AG": ann.get("SpliceAI_pred_DS_AG"),
        "DS_AL": ann.get("SpliceAI_pred_DS_AL"),
        "DS_DG": ann.get("SpliceAI_pred_DS_DG"),
        "DS_DL": ann.get("SpliceAI_pred_DS_DL"),
    }
    spliceai_scores = {k: float(v) for k, v in spliceai_fields.items() 
                      if v not in [None, ""]}
    spliceai_event, spliceai_ds = (
        max(spliceai_scores.items(), key=lambda kv: kv[1]) 
        if spliceai_scores else (None, None)
    )
    
    # HGVS
    hgvsc = ann.get("HGVSc", "")
    hgvsp = ann.get("HGVSp", "")
    if hgvsc and ":" in hgvsc:
        hgvsc = hgvsc.split(":")[1]
    if hgvsp and ":" in hgvsp:
        hgvsp = hgvsp.split(":")[1]
    
    gene = ann.get("SYMBOL")
    hgvs_full = f"{gene or ''} {hgvsc or ''} {hgvsp or ''}".strip() if any([gene, hgvsc, hgvsp]) else None
    
    return {
        "Sample": sample_id,
        "Chrom": chrom, "Pos": pos, "Ref": ref, "Alt": alt,
        "Variant": f"{chrom}:{pos}:{ref}:{alt}",
        "FILTER": filt, "QUAL": qual,
        "Gene": gene,
        "Transcript": ann.get("Feature"),
        "Consequence": ann.get("Consequence"),
        "Impact": ann.get("IMPACT"),
        "Exon": ann.get("EXON"),
        "Intron": ann.get("INTRON"),
        "HGVSc": hgvsc,
        "HGVSp": hgvsp,
        "MANE_ID": ann.get("MANE_SELECT"),
        "GT": f"{gt[0]}/{gt[1]}",
        "AD_Ref": ad_ref, "AD_Alt": ad_alt,
        "DP": dp, "GQ": gq, "VAF": vaf,
        "Zygosity": zyg,
        "ExonCov20": cov.get(">=20x", "NA"),
        "ExonCov30": cov.get(">=30x", "NA"),
        "ExonCov50": cov.get(">=50x", "NA"),
        "ExonCov100": cov.get(">=100x", "NA"),
        "R1": r1r2[0], "R2": r1r2[1], "R1R2_frac": r1r2[2], "R1R2_abs": r1r2_abs,
        "FWD": fr[1], "REV": fr[2], "FWD_frac": fr[3], "REV_frac": fr[4], "FR_abs": fr_abs,
        "gnomAD_AF": ann.get("gnomADg_AF") or ann.get("gnomADe_AF") or ann.get("AF"),
        "REVEL": ann.get("REVEL"),
        "SpliceAI_DS_max": spliceai_ds,
        "SpliceAI_Event": spliceai_event,
        "BayesDel_score": ann.get("BayesDel"),
        "CADD": ann.get("CADD_PHRED") or ann.get("CADDraw"),
        "AM_Pathogenicity": ann.get("am_pathogenicity"),
        "AM_Class": ann.get("am_class"),
        "ClinVar": ann.get("ClinVar_CLNSIG") or ann.get("CLIN_SIG"),
        "ClinVar_ReviewStatus": revstat,
        "ClinVar_Stars": stars,
        "ALLELEID": ann.get("ClinVar_ALLELEID") or ann.get("ALLELEID"),
        "ClinVar_ALLELEID": ann.get("ClinVar_ALLELEID") or ann.get("ALLELEID"),
        "ClinVar_StarsGlyph": clinvar_star_glyph(stars),
        "HGVS_full": hgvs_full
    }


# =============================================================================
# COVERAGE ANALYSIS
# =============================================================================

def read_acmg_thresholds(th_path: Optional[str], sf_genes: Optional[Set[str]] = None) -> pd.DataFrame:
    """Read mosdepth thresholds for ACMG regions."""
    empty_df = pd.DataFrame(columns=[
        "RegionLabel", "Chrom", "ExonStart", "ExonEnd", "ExonLen",
        "Pct>=20x", "Pct>=30x", "Gene", "MANE_ID", "Exon"
    ])
    
    if not th_path or not Path(th_path).exists():
        return empty_df
    
    try:
        df = pd.read_csv(th_path, sep="\t", comment="#", low_memory=False)
        cols = {c.lower(): c for c in df.columns}
        
        def col(name, default=None):
            return cols.get(name.lower(), default)
        
        chrom = col("chrom") or "chrom"
        start = col("start") or "start"
        end = col("end") or "end"
        region = col("region")
        
        # Synthesize RegionLabel if missing
        if not region or region not in df.columns:
            df["RegionLabel"] = (df[chrom].astype(str) + ":" + 
                                df[start].astype(int).astype(str) + "-" + 
                                df[end].astype(int).astype(str))
        else:
            df["RegionLabel"] = df[region].astype(str)
            mask_empty = df["RegionLabel"].isin(["", ".", "unknown", "UNKNOWN", "None"])
            df.loc[mask_empty, "RegionLabel"] = (
                df.loc[mask_empty, chrom].astype(str) + ":" +
                df.loc[mask_empty, start].astype(int).astype(str) + "-" +
                df.loc[mask_empty, end].astype(int).astype(str)
            )
        
        # Parse gene info
        parts = df["RegionLabel"].astype(str).str.split("|", n=2, expand=True)
        df["Gene"] = parts[0] if parts.shape[1] > 0 else "."
        df["MANE_ID"] = parts[1] if parts.shape[1] > 1 else "."
        df["Exon"] = parts[2] if parts.shape[1] > 2 else "."
        
        if sf_genes:
            df = df[df["Gene"].isin(sf_genes)]
        
        df["ExonLen"] = df[end] - df[start]
        
        # Find threshold columns
        if "20X" not in df.columns or "30X" not in df.columns:
            for t in ("20", "30"):
                cand = [c for c in df.columns if "".join(ch for ch in str(c) if ch.isdigit()) == t]
                if cand:
                    df[f"{t}X"] = df[cand[0]]
        
        df["Pct>=20x"] = (df["20X"] / df["ExonLen"] * 100).round(2)
        df["Pct>=30x"] = (df["30X"] / df["ExonLen"] * 100).round(2)
        
        out = df.rename(columns={chrom: "Chrom", start: "ExonStart", end: "ExonEnd"})[
            ["RegionLabel", "Gene", "MANE_ID", "Exon", "Chrom", "ExonStart", "ExonEnd",
             "ExonLen", "Pct>=20x", "Pct>=30x"]
        ].drop_duplicates(subset=["RegionLabel", "Chrom", "ExonStart", "ExonEnd"])
        
        return out
    
    except Exception as e:
        print(f"[ERROR] Failed to read ACMG thresholds: {e}", file=sys.stderr)
        return empty_df


def load_gap_bed_for_agg(path: Optional[str], threshold_tag: str, 
                        sf_genes: Optional[Set[str]] = None) -> pd.DataFrame:
    """Load and aggregate gap BED file."""
    empty_df = pd.DataFrame(columns=[
        "RegionLabel", f"Gaps{threshold_tag}_n", 
        f"Gaps{threshold_tag}_bp", f"Gaps{threshold_tag}_intervals"
    ])
    
    if not path or not Path(path).exists():
        return empty_df
    
    try:
        df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
        
        if df.shape[1] >= 4:
            df.columns = ["Chrom", "Start", "End", "RegionLabel"] + [
                f"extra{i}" for i in range(df.shape[1] - 4)
            ]
        else:
            df.columns = ["Chrom", "Start", "End"]
            df["RegionLabel"] = (df["Chrom"].astype(str) + ":" + 
                                df["Start"].astype(int).astype(str) + "-" + 
                                df["End"].astype(int).astype(str))
        
        df["GapLen"] = pd.to_numeric(df["End"], errors="coerce") - pd.to_numeric(df["Start"], errors="coerce")
        df = df.dropna(subset=["Start", "End", "GapLen"])
        df = df.drop_duplicates(subset=["Chrom", "Start", "End", "RegionLabel"])
        
        def intervals_formatter(g):
            return ";".join(
                f"{r.Chrom}:{int(r.Start)}-{int(r.End)}" 
                for _, r in g.sort_values(["Chrom", "Start", "End"]).iterrows()
            )
        
        agg = df.groupby("RegionLabel", as_index=False).agg(**{
            f"Gaps{threshold_tag}_n": ("GapLen", "size"),
            f"Gaps{threshold_tag}_bp": ("GapLen", "sum"),
        })
        
        ints = (df.groupby("RegionLabel")
                .apply(intervals_formatter)
                .reset_index(name=f"Gaps{threshold_tag}_intervals"))
        
        out = agg.merge(ints, on="RegionLabel", how="left")
        out[f"Gaps{threshold_tag}_n"] = out[f"Gaps{threshold_tag}_n"].astype("Int64")
        out[f"Gaps{threshold_tag}_bp"] = out[f"Gaps{threshold_tag}_bp"].astype("Int64")
        
        return out
    
    except Exception as e:
        print(f"[ERROR] Failed to load gaps from {path}: {e}", file=sys.stderr)
        return empty_df


def build_genes_coverage(df_in: pd.DataFrame, sf_genes: Optional[Set[str]] = None) -> pd.DataFrame:
    """Build per-gene coverage summary."""
    required = ["Gene", "MANE_ID", "ExonLen", "Pct>=20x"]
    missing = [c for c in required if c not in df_in.columns]
    if missing:
        print(f"[WARN] Genes Coverage: missing columns: {', '.join(missing)}", file=sys.stderr)
        return pd.DataFrame(columns=["Gene", "MANE_ID", "Pct20"])
    
    df = df_in[required].copy()
    
    if sf_genes:
        df = df[df["Gene"].isin(sf_genes)]
    
    df["Gene"] = df["Gene"].astype(str).str.strip()
    df = df[df["Gene"].notna() & (~df["Gene"].isin(["", "NA", "."]))]
    
    df["MANE_ID"] = df["MANE_ID"].astype(str).str.strip()
    df["ExonLen"] = pd.to_numeric(df["ExonLen"], errors="coerce")
    
    p20 = (df["Pct>=20x"].astype(str)
           .str.replace("%", "", regex=False)
           .str.replace(",", ".", regex=False))
    df["Pct>=20x"] = pd.to_numeric(p20, errors="coerce")
    
    df = df.dropna(subset=["ExonLen", "Pct>=20x"])
    df = df[df["ExonLen"] > 0]
    
    if df.empty:
        return pd.DataFrame(columns=["Gene", "MANE_ID", "Pct20"])
    
    # Choose best transcript per gene
    tx_pref = (df[~df["MANE_ID"].isin(["", "NA", "."])]
               .groupby(["Gene", "MANE_ID"], as_index=False)["ExonLen"].sum()
               .sort_values(["Gene", "ExonLen"], ascending=[True, False]))
    
    best_tx = tx_pref.groupby("Gene", as_index=False).first()[["Gene", "MANE_ID"]]
    
    all_genes = pd.DataFrame({"Gene": sorted(df["Gene"].unique())})
    best_tx = all_genes.merge(best_tx, on="Gene", how="left").fillna({"MANE_ID": "NA"})
    
    # Length-weighted mean coverage
    df["_covered20_bp"] = df["ExonLen"] * (df["Pct>=20x"] / 100.0)
    agg = df.groupby("Gene", as_index=False).agg(
        total_len=("ExonLen", "sum"),
        covered20_bp=("_covered20_bp", "sum")
    )
    agg["Pct20"] = (agg["covered20_bp"] / agg["total_len"] * 100.0).round(2)
    
    out = (agg.merge(best_tx, on="Gene", how="left")
           [["Gene", "MANE_ID", "Pct20"]]
           .sort_values("Gene")
           .reset_index(drop=True))
    
    return out


# =============================================================================
# EXCEL REPORT BUILDER
# =============================================================================

def build_sample_summary(sample_id: str, args, df: pd.DataFrame) -> Dict:
    """Build sample summary metrics."""
    ss = {
        "Sample_ID": sample_id,
        "Assay": args.assay,
        "Build": args.build
    }
    
    # QC metrics
    ss.update(parse_flagstat(args.flagstat))
    ss.update(parse_samtools_stats(args.stats))
    ss.update(parse_mosdepth_summary_regions(args.mosdepth_summary))
    ss.update(acmg_pct_regions_covered(args.acmg_thresholds))
    ss.update(parse_sexcheck(args.sexcheck))
    
    # Derive Ti/Tv from variants if not provided
    if "TiTv" not in ss:
        try:
            snps = df[(df["Ref"].str.len() == 1) & (df["Alt"].str.len() == 1)]
            transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
            ti = sum((r.Ref, r.Alt) in transitions for _, r in snps.iterrows())
            tv = len(snps) - ti
            if tv > 0:
                ss["TiTv"] = round(ti / tv, 3)
        except Exception:
            pass
    
    # Het/Hom ratio from VCF
    if "Het_Hom_Ratio" not in ss:
        try:
            het = (df["Zygosity"] == "het").sum()
            hom = (df["Zygosity"] == "hom").sum()
            if hom > 0:
                ss["Het_Hom_Ratio"] = round(het / hom, 3)
        except Exception:
            pass
    
    return ss


def build_coverage_gaps_table(args, sf_genes: Optional[Set[str]]) -> pd.DataFrame:
    """Build combined coverage gaps table."""
    try:
        th_exon = read_acmg_thresholds(args.acmg_thresholds, sf_genes=sf_genes)
        g20_agg = load_gap_bed_for_agg(args.gaps20, "<20x", sf_genes=sf_genes)
        g30_agg = load_gap_bed_for_agg(args.gaps30, "<30x", sf_genes=sf_genes)
        
        combined = (th_exon
                   .merge(g20_agg, on="RegionLabel", how="left")
                   .merge(g30_agg, on="RegionLabel", how="left"))
        
        # Ensure all expected columns exist
        for c in ["Gaps<20x_n", "Gaps<20x_bp", "Gaps<20x_intervals",
                  "Gaps<30x_n", "Gaps<30x_bp", "Gaps<30x_intervals"]:
            if c not in combined.columns:
                combined[c] = pd.NA
        
        return combined[COVERAGE_GAP_COLS].sort_values(["Chrom", "ExonStart", "ExonEnd"])
    
    except Exception as e:
        print(f"[ERROR] Failed to build coverage gaps table: {e}", file=sys.stderr)
        return pd.DataFrame(columns=COVERAGE_GAP_COLS)


def write_excel_report(df: pd.DataFrame, args, sample_id: str, sf_genes: Set[str]):
    """Write all sheets to Excel file."""
    with pd.ExcelWriter(args.xlsx_out) as xw:
        # Sheet 1: Sample Summary
        ss = build_sample_summary(sample_id, args, df)
        pd.DataFrame([ss]).to_excel(xw, index=False, sheet_name="Sample Summary")
        
        # Sheet 2: ACMG SF (P/LP)
        acmg = df.copy()
        acmg["ClinVar_URL"] = acmg.apply(clinvar_url_from_row, axis=1)
        acmg["ClinVar_Link"] = acmg["ClinVar_URL"].apply(
            lambda u: f'=HYPERLINK("{u}","{u}")' if u else ""
        )
        
        if "ClinVar" in acmg.columns:
            acmg = acmg[acmg["ClinVar"].astype(str).str.contains(
                "pathogenic", case=False, na=False
            )]
        else:
            acmg = acmg.iloc[0:0]
        
        if sf_genes:
            acmg = acmg[acmg["Gene"].isin(sf_genes)]
        
        acmg[ACMG_SHEET_COLS].to_excel(xw, index=False, sheet_name="ACMG SF (P-LP)")
        
        # Sheet 3: Coverage gaps
        combined_df = build_coverage_gaps_table(args, sf_genes if sf_genes else None)
        combined_df.to_excel(xw, index=False, sheet_name="Coverage gaps")
        
        # Sheet 4: ACMG SF Genes Coverage
        try:
            genes_cov = build_genes_coverage(
                combined_df, 
                sf_genes=sf_genes if sf_genes else None
            )
        except Exception as e:
            print(f"[WARN] Failed to build genes coverage: {e}", file=sys.stderr)
            genes_cov = pd.DataFrame(columns=["Gene", "MANE_ID", "Pct20"])
        
        genes_cov.to_excel(xw, index=False, sheet_name="ACMG SF Genes Coverage")
        
        # Sheet 5: PASS variants
        pass_df = df[df["FILTER"] == "PASS"][PASS_VARIANT_COLS]
        pass_df.to_excel(xw, index=False, sheet_name="PASS variants")
    
    print(f"[INFO] Wrote Excel report → {args.xlsx_out}")


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    """Parse command line arguments."""
    p = argparse.ArgumentParser(
        description="Build LEAN Excel deliverable with variant QC and coverage analysis"
    )
    
    # Required inputs
    p.add_argument("vcf", help="VCF (bgzip+tabix) with VEP+ClinVar annotations")
    p.add_argument("coverage_summary", help="Exon coverage summary file")
    p.add_argument("r1r2_summary", help="R1/R2 ratio per exon TSV")
    p.add_argument("fr_summary", help="Forward/reverse strand ratio per exon TSV")
    p.add_argument("xlsx_out", help="Output Excel path")
    
    # Sample metadata
    p.add_argument("--sample-id", help="Sample ID (default: infer from VCF)")
    p.add_argument("--assay", default="", help="Assay type (WES/WGS)")
    p.add_argument("--build", default="GRCh38", help="Reference genome build")
    
    # QC inputs
    p.add_argument("--flagstat", help="samtools flagstat output")
    p.add_argument("--stats", help="samtools stats output")
    p.add_argument("--picard-align", help="Picard CollectAlignmentSummaryMetrics.txt")
    p.add_argument("--picard-insert", help="Picard CollectInsertSizeMetrics.txt")
    p.add_argument("--mosdepth-summary", help="mosdepth summary.txt")
    p.add_argument("--verifybamid", help="verifyBamID2 output (FREEMIX)")
    p.add_argument("--sexcheck", help="Sex check output file")
    
    # Gene filtering
    p.add_argument("--sf-genes", help="ACMG SF gene list (one per line)")
    
    # Coverage analysis
    p.add_argument("--acmg-thresholds", help="mosdepth thresholds.bed(.gz)")
    p.add_argument("--gaps20", help="Annotated gaps <20x BED")
    p.add_argument("--gaps30", help="Annotated gaps <30x BED")
    
    return p.parse_args()


def main():
    """Main execution."""
    args = parse_args()
    
    # Load gene list
    sf_genes = load_gene_list(args.sf_genes)
    
    # Determine sample ID
    sample_id = args.sample_id or sample_from_vcf(args.vcf)
    print(f"[INFO] Processing sample: {sample_id}")
    
    # Load coverage data
    print(f"[INFO] Loading coverage data...")
    coverage_data = CoverageData(
        args.coverage_summary,
        args.r1r2_summary,
        args.fr_summary
    )
    
    # Parse VCF
    print(f"[INFO] Parsing VCF: {args.vcf}")
    df = parse_vcf_to_dataframe(args.vcf, coverage_data, sample_id)
    print(f"[INFO] Found {len(df)} variants")
    
    # Write Excel report
    print(f"[INFO] Writing Excel report...")
    write_excel_report(df, args, sample_id, sf_genes)
    
    print(f"[INFO] Complete!")


if __name__ == "__main__":
    main()