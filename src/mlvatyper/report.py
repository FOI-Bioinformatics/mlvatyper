"""Self-contained HTML run report.

The report is built only from files in the output directory, so it can be
regenerated after editing definitions or rounding settings without a new
search. Styles and figures are embedded; the file needs no network access
and prints to PDF from a browser.
"""

import html
import io
import json
from dataclasses import dataclass, field
from datetime import datetime, UTC
from pathlib import Path

import pandas as pd

from . import __version__
from .log import get_logger

log = get_logger()

STATUS_GLOSSARY = {
    "OK": "amplicon found and allele within tolerance of an integer",
    "NOT_DETECTED": "no amplicon; neither primer matched exactly",
    "NEAR_CONTIG_END": "no amplicon, but a primer lies within the size window of a contig end",
    "NOT_SEARCHED": "no primer holds a hash word of the word size used, so merPCR did not search the marker",
    "ONE_PRIMER_ONLY": "one primer matches exactly, the other exceeds the mismatch limit or has an indel",
    "INTERMEDIATE": "repeat number not within tolerance of an integer",
    "NEGATIVE": "amplicon shorter than the flanking region of the definition",
    "BINNED": "allele taken from the bin table",
    "OFF_TABLE": "size outside every bin; the formula was used",
    "MULTIPLE": "several non-overlapping amplicons in the sample",
    "AMBIGUOUS": "overlapping hits for the locus were collapsed to the best one",
    "PRIMER_INDEL": "amplicon found by the primer alignment with one indel in a primer",
    "CENSORED": "array at least as long as --max-array-length; reported as a lower bound",
    "HAS_N": "amplicon contains N",
    "COUNT_DISAGREES": "tandem array count differs from the length based allele",
    "NO_DEFINITION": "no unit length or reference values for the marker",
}

# muted colours for flagged cells; OK stays white
STATUS_COLOURS = {
    "OK": "#ffffff", "MULTIPLE": "#fff3cd", "AMBIGUOUS": "#fff3cd", "INTERMEDIATE": "#ffe0b2",
    "COUNT_DISAGREES": "#ffe0b2", "PRIMER_INDEL": "#e3f2fd", "HAS_N": "#e3f2fd", "BINNED": "#e8f5e9",
    "NOT_DETECTED": "#eeeeee", "NEAR_CONTIG_END": "#eeeeee", "NOT_SEARCHED": "#eeeeee",
    "ONE_PRIMER_ONLY": "#eeeeee", "NO_DEFINITION": "#eeeeee", "OFF_TABLE": "#f8d7da", "NEGATIVE": "#f8d7da",
    "CENSORED": "#f8d7da",
}
ACCENT = "#1f5f8b"


@dataclass
class ReportData:
    outdir: Path
    prefix: str
    title: str
    run_info: dict
    samples: list[str]
    markers: list[str]
    alleles: pd.DataFrame
    status: pd.DataFrame
    profile: pd.DataFrame
    calls: pd.DataFrame
    missing: pd.DataFrame
    definitions: pd.DataFrame | None
    discovery: pd.DataFrame | None
    figures: dict[str, str] = field(default_factory=dict)

    @property
    def n_amplicons(self) -> int:
        return int((self.calls["copy"] != "").sum())


def _read(path: Path, **kwargs) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def load_report_data(outdir: Path, prefix: str, title: str | None = None) -> ReportData:
    """Load the tables of one run; raises FileNotFoundError when the calls are missing."""
    outdir = Path(outdir)
    calls_path = outdir / f"{prefix}_calls_long.tsv"
    if not calls_path.exists():
        raise FileNotFoundError(f"{calls_path} not found; run 'mlvatyper summarize' first")
    run_info_path = outdir / "run_info.json"
    run_info = json.loads(run_info_path.read_text(encoding="utf-8")) if run_info_path.exists() else {}
    calls = pd.read_csv(calls_path, sep="\t", dtype=str, keep_default_na=False)
    calls = calls.replace({"NA": ""})
    alleles = _read(outdir / f"{prefix}_alleles.tsv", index_col=0)
    status = _read(outdir / f"{prefix}_status.tsv", index_col=0)
    profile = _read(outdir / f"{prefix}_profile.tsv", index_col=0)
    if alleles is None or status is None or profile is None:
        raise FileNotFoundError(f"{prefix}_alleles.tsv, _status.tsv and _profile.tsv are needed in {outdir}")
    samples = [str(s) for s in alleles.index]
    markers = [str(m) for m in alleles.columns]
    missing = _read(outdir / "missing_markers.tsv")
    if missing is None:
        missing = pd.DataFrame(columns=["sample", "marker", "reason"])
    definitions = _read(outdir / "marker_defs_derived.tsv")
    if definitions is None:
        markers_file = (run_info.get("summarize", {}).get("parameters", {}) or {}).get("markers")
        if markers_file and Path(markers_file).exists():
            try:
                definitions = _read(Path(markers_file))
            except Exception:  # pragma: no cover - a legacy file without header
                definitions = None
    discovery = _read(outdir / "discover_report.tsv")
    scheme = scheme_name(run_info) or prefix
    return ReportData(outdir, prefix, title or f"MLVA typing report: {scheme}", run_info, samples, markers,
                      alleles, status, profile, calls, missing, definitions, discovery)


def scheme_name(run_info: dict) -> str:
    """The scheme used: its --scheme name, else the STS file stem."""
    search = run_info.get("search", {}).get("parameters", {}) or {}
    scheme = search.get("scheme")
    if scheme:
        return Path(str(scheme)).name
    sts = search.get("sts")
    return Path(str(sts)).stem if sts else ""


def status_counts(data: ReportData) -> pd.DataFrame:
    """Status code counts per sample, codes as columns in a fixed order."""
    counts: dict[str, dict[str, int]] = {}
    for sample in data.samples:
        counts[sample] = {}
        for marker in data.markers:
            for code in str(data.status.loc[sample, marker]).split(";"):
                counts[sample][code] = counts[sample].get(code, 0) + 1
    df = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    order = [c for c in STATUS_GLOSSARY if c in df.columns] + [c for c in df.columns if c not in STATUS_GLOSSARY]
    return df.reindex(columns=order)


def failures(data: ReportData) -> list[tuple[str, str, list[str]]]:
    """(code, explanation, 'sample: marker' items) for every flagged code present."""
    found: dict[str, list[str]] = {}
    for sample in data.samples:
        for marker in data.markers:
            for code in str(data.status.loc[sample, marker]).split(";"):
                if code != "OK":
                    found.setdefault(code, []).append(f"{sample}: {marker}")
    return [(code, STATUS_GLOSSARY.get(code, ""), items) for code, items in found.items()]


def cell_colour(status_text: str) -> str:
    codes = str(status_text).split(";")
    for code in codes:
        if code in STATUS_COLOURS and code != "OK":
            return STATUS_COLOURS[code]
    return STATUS_COLOURS["OK"]


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def figures(data: ReportData) -> dict[str, str]:
    """SVG figures: amplicon lengths, residuals and the status matrix."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8, "font.family": "sans-serif", "svg.fonttype": "none"})
    out: dict[str, str] = {}
    found = data.calls[data.calls["copy"] != ""].copy()
    found["length"] = _numeric(found["length"])
    found["expected_size"] = _numeric(found["expected_size"])
    found["residual_bp"] = _numeric(found["residual_bp"])
    markers = data.markers
    pos = {m: i for i, m in enumerate(markers)}
    width = max(6.0, 0.28 * len(markers))

    fig, ax = plt.subplots(figsize=(width, 3.2))
    exp = found.drop_duplicates("marker")
    ax.scatter([pos[m] for m in exp["marker"]], exp["expected_size"], marker="_", s=120, color="#999999",
               label="STS expected size")
    ax.scatter([pos[m] for m in found["marker"]], found["length"], s=12, color=ACCENT, alpha=0.7, label="amplicon")
    ax.set_xticks(range(len(markers)))
    ax.set_xticklabels(markers, rotation=90)
    ax.set_ylabel("length (bp)")
    ax.legend(frameon=False, loc="upper left", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    out["lengths"] = _svg(fig)

    fig, ax = plt.subplots(figsize=(width, 2.6))
    with_res = found[found["residual_bp"].notna()]
    ax.axhline(0, color="#999999", linewidth=0.8)
    ax.scatter([pos[m] for m in with_res["marker"]], with_res["residual_bp"], s=12, color=ACCENT, alpha=0.7)
    ax.set_xticks(range(len(markers)))
    ax.set_xticklabels(markers, rotation=90)
    ax.set_ylabel("residual (bp)")
    limit = max(2.0, float(with_res["residual_bp"].abs().max()) if len(with_res) else 2.0) * 1.2
    ax.set_ylim(-limit, limit)
    ax.spines[["top", "right"]].set_visible(False)
    out["residuals"] = _svg(fig)

    codes = list(STATUS_COLOURS)
    matrix = [[codes.index(_first_flag(data.status.loc[s, m])) for m in markers] for s in data.samples]
    fig, ax = plt.subplots(figsize=(width, max(1.5, 0.25 * len(data.samples) + 1.2)))
    cmap = matplotlib.colors.ListedColormap([STATUS_COLOURS[c] for c in codes])
    ax.imshow(matrix, cmap=cmap, vmin=0, vmax=len(codes) - 1, aspect="auto")
    ax.set_xticks(range(len(markers)))
    ax.set_xticklabels(markers, rotation=90)
    ax.set_yticks(range(len(data.samples)))
    ax.set_yticklabels(data.samples)
    for spine in ax.spines.values():
        spine.set_visible(False)
    out["status"] = _svg(fig)
    return out


def _first_flag(status_text: str) -> str:
    for code in str(status_text).split(";"):
        if code != "OK" and code in STATUS_COLOURS:
            return code
    return "OK"


def _svg(fig) -> str:
    import matplotlib.pyplot as plt
    fig.tight_layout()
    buffer = io.StringIO()
    fig.savefig(buffer, format="svg")
    plt.close(fig)
    text = buffer.getvalue()
    return text[text.index("<svg"):]


def methods_text(run_info: dict) -> str:
    """A methods paragraph generated from the recorded parameters."""
    search = run_info.get("search", {}).get("parameters", {}) or {}
    summ = run_info.get("summarize", {}).get("parameters", {}) or {}
    merpcr = run_info.get("merpcr_version", "unknown")
    version = run_info.get("mlvatyper_version", __version__)
    parts = [
        f"VNTR amplicons were located in the assemblies with mlvatyper {version} "
        f"(https://github.com/FOI-Bioinformatics/mlvatyper) using merPCR {merpcr}, a Python reimplementation "
        f"of me-PCR, allowing {search.get('mismatches', 2)} mismatches per primer outside the "
        f"{search.get('three_prime_match', 1)} protected 3' base(s), a hash word of {search.get('wordsize', 11)} "
        f"bases and a product size window of the expected size plus or minus {search.get('margin', 500)} bp.",
    ]
    if search.get("iupac", True):
        parts.append("IUPAC ambiguity codes in primers were matched to the bases they stand for.")
    if search.get("rescue_wordsize", 0):
        parts.append(
            f"Markers without a hit were searched again with a hash word of {search['rescue_wordsize']} bases, "
            f"one-sided search and a window of {search.get('rescue_margin') or 2 * search.get('margin', 500)} bp"
            + (", and then by anchored alignment allowing one indel in a primer." if search.get("rescue_indels")
               else "."))
    parts.append(
        "Overlapping hits for one locus were collapsed to the hit with the fewest primer mismatches. "
        "Repeat numbers were calculated as reference repeats plus the difference between the amplicon length "
        "and the reference size divided by the unit length"
        + (f", with reference values taken from sample {summ['reference_sample']}" if summ.get("reference_sample")
           else "")
        + f". Values within {summ.get('allele_tolerance', 0.25)} of an integer were assigned that integer; other "
        f"values were handled with the {summ.get('rounding', 'intermediate')} policy. Each call was cross-checked "
        f"against a tandem array count tolerating {summ.get('unit_mismatches', 1)} mismatch(es) per unit."
    )
    if summ.get("bins"):
        parts.append("Amplicon sizes were first mapped to alleles with the bin table given.")
    return " ".join(parts)


def _command_line(run_info: dict) -> str:
    search = run_info.get("search", {}).get("parameters", {}) or {}
    items = []
    for key, value in sorted(search.items()):
        if value in (None, False, "", []) or key in ("report", "title"):
            continue
        flag = "--" + key.replace("_", "-")
        if value is True:
            items.append(flag)
        elif isinstance(value, list):
            items.append(f"{flag} {' '.join(str(v) for v in value)}")
        else:
            items.append(f"{flag} {value}")
    return "mlvatyper run " + " ".join(items)


EVIDENCE_COLUMNS = {
    "marker": "Marker", "copy": "Copy", "contig": "Contig", "start": "Start", "end": "End", "strand": "Strand",
    "length": "Length (bp)", "expected_size": "Expected (bp)", "residual_bp": "Residual (bp)",
    "fwd_mismatches": "Fwd mism.", "rev_mismatches": "Rev mism.", "fwd_indels": "Fwd indels",
    "rev_indels": "Rev indels", "search_pass": "Search pass", "allele": "Allele", "tandem_copies": "Tandem copies",
    "tandem_exact": "Exact", "tandem_variant": "Variant", "status": "Status",
}


def _table_html(df: pd.DataFrame, index: bool = True, colour: pd.DataFrame | None = None,
                headers: dict[str, str] | None = None, row_status: str | None = None) -> str:
    """Plain HTML table.

    colour gives a background per cell from a status table of the same
    shape; row_status names a column whose status shades the whole row;
    headers renames columns for display.
    """
    cols = list(df.columns)
    names = headers or {}
    rows = ["<table>", "<thead><tr>" + ("<th></th>" if index else "")
            + "".join(f"<th>{html.escape(str(names.get(c, c)))}</th>" for c in cols) + "</tr></thead>", "<tbody>"]
    for idx, row in df.iterrows():
        cells = []
        row_style = ""
        if row_status is not None:
            shade = cell_colour(row[row_status])
            row_style = f' style="background:{shade}"' if shade != STATUS_COLOURS["OK"] else ""
        for c in cols:
            style = ""
            if colour is not None:
                style = f' style="background:{cell_colour(colour.loc[idx, c])}"'
            cells.append(f"<td{style}>{html.escape(str(row[c]))}</td>")
        rows.append(f"<tr{row_style}>" + (f"<th>{html.escape(str(idx))}</th>" if index else "")
                    + "".join(cells) + "</tr>")
    rows.extend(["</tbody>", "</table>"])
    return "\n".join(rows)


def _integer_text(value: str) -> str:
    """'0.0' as '0', other text unchanged."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return str(int(number)) if number == int(number) else str(value)


def legend_html() -> str:
    """Colour swatches for the flagged status codes."""
    groups = [
        ("several or ambiguous loci", "MULTIPLE"), ("intermediate or count disagreement", "INTERMEDIATE"),
        ("primer indel or N", "PRIMER_INDEL"), ("binned", "BINNED"), ("missing or no definition", "NOT_DETECTED"),
        ("outside bins, negative or censored", "OFF_TABLE"),
    ]
    items = "".join(
        f'<span class="swatch" style="background:{STATUS_COLOURS[code]}"></span>{html.escape(label)}'
        for label, code in groups
    )
    return f'<p class="legend"><span class="swatch" style="background:#ffffff"></span>OK {items}</p>'


SOFTWARE = [
    ("mlvatyper", "https://github.com/FOI-Bioinformatics/mlvatyper", "this report"),
    ("merPCR", "https://github.com/FOI-Bioinformatics/merpcr",
     "reimplementation of me-PCR: Schuler GD. Sequence mapping by electronic PCR. Genome Res 1997;7:541-550"),
    ("Tandem Repeats Finder", "https://tandem.bu.edu/trf/trf.html",
     "Benson G. Nucleic Acids Res 1999;27:573-580 (only when discover used the trf engine)"),
    ("MAFFT", "https://mafft.cbrc.jp/alignment/software/",
     "Katoh K, Standley DM. Mol Biol Evol 2013;30:772-780 (only when alignments were produced)"),
]


def _summary_table(data: ReportData) -> pd.DataFrame:
    """Per sample: markers typed (integer allele), flagged and missing, then the status counts."""
    counts = status_counts(data)
    typed = []
    missing = []
    for sample in data.samples:
        row = data.alleles.loc[sample]
        typed.append(int(sum(str(v).lstrip("-").isdigit() for v in row)))
        missing.append(int(sum(str(v) in STATUS_GLOSSARY for v in row)))
    counts.insert(0, "markers typed", typed)
    counts.insert(1, "missing", missing)
    return counts


def _discovery_summary(discovery: pd.DataFrame | None) -> str:
    if discovery is None or "category" not in discovery.columns:
        return ""
    found = discovery[discovery["category"] != ""]
    cats = found["category"].value_counts().to_dict()
    mixed = int(pd.to_numeric(discovery.get("mixed", pd.Series(dtype=str)), errors="coerce").fillna(0).sum())
    parts = [f"{len(found)} of {len(discovery)} markers have a candidate unit"]
    parts.append(", ".join(f"{cats.get(c, 0)} in category {c}" for c in ("A", "B", "C")))
    if mixed:
        parts.append(f"{mixed} mixed")
    return "; ".join(parts) + "."


def render(data: ReportData) -> str:
    from jinja2 import Environment, PackageLoader, select_autoescape

    env = Environment(loader=PackageLoader("mlvatyper", "templates"), autoescape=select_autoescape(["html"]))
    template = env.get_template("report.html.j2")
    if not data.figures:
        data.figures = figures(data)
    evidence_cols = list(EVIDENCE_COLUMNS)
    found = data.calls[data.calls["copy"] != ""].copy()
    found["residual_bp"] = found["residual_bp"].map(_integer_text)
    evidence = {s: _table_html(found[found["sample"] == s][evidence_cols], index=False, headers=EVIDENCE_COLUMNS,
                               row_status="status") for s in data.samples}
    missing = data.missing.copy()
    if len(missing):
        missing["explanation"] = missing["reason"].map(lambda r: STATUS_GLOSSARY.get(r, ""))
    transposed = len(data.markers) > len(data.samples)
    alleles, status = (data.alleles.T, data.status.T) if transposed else (data.alleles, data.status)
    search = data.run_info.get("search", {}).get("parameters", {}) or {}
    summ = data.run_info.get("summarize", {}).get("parameters", {}) or {}
    meta = [
        ("Scheme", scheme_name(data.run_info) or data.prefix),
        ("Samples", str(len(data.samples))), ("Markers", str(len(data.markers))),
        ("Amplicons", str(data.n_amplicons)),
        ("Reference sample", summ.get("reference_sample") or search.get("reference_sample") or "none"),
        ("Primers", Path(str(search.get("sts", ""))).name or "unknown"),
        ("Marker definitions", Path(str(summ.get("markers") or search.get("markers") or "")).name or "none"),
        ("Search date", (data.run_info.get("search", {}).get("timestamp") or "")[:16].replace("T", " ")),
        ("Output directory", str(data.outdir)),
    ]
    return template.render(
        data=data, generated=datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"), version=__version__,
        merpcr_version=data.run_info.get("merpcr_version", "unknown"), command=_command_line(data.run_info),
        meta=meta, summary_table=_table_html(_summary_table(data)),
        alleles_table=_table_html(alleles, colour=status), transposed=transposed, legend=legend_html(),
        profile_table=_table_html(data.profile[["profile"]]), evidence=evidence, failures=failures(data),
        missing_table=_table_html(missing, index=False) if len(missing) else "",
        discovery_table=_table_html(data.discovery, index=False) if data.discovery is not None else "",
        discovery_summary=_discovery_summary(data.discovery),
        definitions_table=_table_html(data.definitions, index=False) if data.definitions is not None else "",
        methods=methods_text(data.run_info), inputs=data.run_info.get("inputs", []), glossary=STATUS_GLOSSARY,
        software=SOFTWARE, prefix=data.prefix,
        parameters=json.dumps({k: v for k, v in data.run_info.items() if k not in ("inputs",)}, indent=2,
                              sort_keys=True),
    )


def write_report(outdir: Path, prefix: str, title: str | None = None) -> Path:
    data = load_report_data(outdir, prefix, title)
    path = Path(outdir) / f"{prefix}_report.html"
    path.write_text(render(data), encoding="utf-8")
    log.info("Wrote %s", path)
    return path
