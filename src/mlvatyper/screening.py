"""Summaries across many samples: call rate per sample and per locus, and loci that need a second look.

The tables are built from the status and allele tables of a run. They are written as
``<prefix>_sample_summary.tsv`` and ``<prefix>_locus_summary.tsv`` and shown in the report. A locus
is listed as a concern when more than half of its calls carry the same flag, or when an allele is
zero or below; either points at the scheme (reference values, units or patterns) more often than at
the strains.
"""

import re

import pandas as pd

ABSENT = ("NOT_DETECTED", "NEAR_CONTIG_END", "NOT_SEARCHED", "ONE_PRIMER_ONLY")
NUMBER = re.compile(r"^-?\d+(\.\d+)?$")
FLAG_SHARE = 0.5  # a flag on more than this share of the typed calls of a locus is reported as a concern
LOW_CALL_RATE = 0.75  # samples with a smaller share of OK loci are listed as low call rate
MIN_SAMPLES = 5  # the overview is meant for sets, not for single samples


def _is_number(text: str) -> bool:
    return bool(NUMBER.match(str(text)))


def sample_summary(alleles: pd.DataFrame, status: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n = alleles.shape[1]
    for sample in alleles.index:
        a, s = alleles.loc[sample], status.loc[sample]
        typed = sum(_is_number(v) for v in a)
        ok = sum(v == "OK" for v in s)
        missing = sum(v in ABSENT for v in s)
        rows.append({"sample": sample, "loci": n, "typed": typed, "ok": ok, "flagged": typed - ok,
                     "missing": missing, "ok_fraction": round(ok / n, 3) if n else 0.0})
    return pd.DataFrame(rows).set_index("sample")


def _simpson(values: list[str]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    counts = pd.Series(values).value_counts()
    return round(1 - float(((counts / n) ** 2).sum()), 3)


def locus_summary(alleles: pd.DataFrame, status: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for marker in alleles.columns:
        a, s = alleles[marker], status[marker]
        numeric = [str(v) for v in a if _is_number(v)]
        typed = len(numeric)
        ok = int((s == "OK").sum())
        missing = int(s.isin(ABSENT).sum())
        codes: dict[str, int] = {}
        for sample_status, value in zip(s, a, strict=True):
            if not _is_number(value):
                continue
            for code in str(sample_status).split(";"):
                if code != "OK":
                    codes[code] = codes.get(code, 0) + 1
        top, top_n = max(codes.items(), key=lambda kv: kv[1]) if codes else ("", 0)
        lowest = min((float(v) for v in numeric), default=float("nan"))
        concern = []
        if typed and top_n / typed > FLAG_SHARE:
            concern.append(f"{top} in {100 * top_n / typed:.0f}% of calls")
        if typed and lowest <= 0:
            concern.append(f"allele {lowest:g}")
        share = round(top_n / typed, 3) if typed else 0.0
        rows.append({"marker": marker, "samples": len(a), "typed": typed, "ok": ok, "missing": missing,
                     "flagged": typed - ok, "top_flag": top, "top_flag_share": share,
                     "distinct_alleles": len(set(numeric)), "simpson": _simpson(numeric), "lowest_allele": lowest,
                     "concern": "; ".join(concern)})
    return pd.DataFrame(rows).set_index("marker")


def missing_by_reason(missing: pd.DataFrame) -> pd.DataFrame:
    if missing is None or not len(missing):
        return pd.DataFrame(columns=["reason", "products", "markers"])
    grouped = missing.groupby("reason")
    return pd.DataFrame({"products": grouped.size(),
                         "markers": grouped["marker"].agg(lambda m: ", ".join(
                             f"{k} ({v})" for k, v in m.value_counts().head(5).items()))}).sort_values(
        "products", ascending=False).reset_index()


def overview(alleles: pd.DataFrame, status: pd.DataFrame, missing: pd.DataFrame | None) -> dict:
    """Everything the report shows for a set of samples."""
    samples = sample_summary(alleles, status)
    loci = locus_summary(alleles, status)
    low = samples[samples["ok_fraction"] < LOW_CALL_RATE]
    return {"samples": samples, "loci": loci, "low_call_rate": list(low.index),
            "concerns": loci[loci["concern"] != ""], "missing_by_reason": missing_by_reason(missing),
            "ok_fraction_mean": float(samples["ok_fraction"].mean()),
            "ok_fraction_min": float(samples["ok_fraction"].min())}
