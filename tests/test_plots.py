"""Clustering and the multi sample plots."""

import numpy as np
import pandas as pd

from mlvatyper.plots import (
    MIN_SHARED, allele_matrix, allele_spectra, call_rate_bars, clustered_heatmap, dendrogram_segments,
    distance_heatmap, distance_matrix, leaf_order, upgma,
)


def _alleles():
    loci = [f"L{i}" for i in range(MIN_SHARED + 2)]
    base = [str(5 + i) for i in range(len(loci))]
    rows = {
        "a1": base, "a2": base,                                    # identical pair
        "b1": [str(int(x) + 1) for x in base], "b2": [str(int(x) + 1) for x in base],
        "c": [str(int(x) + 3) for x in base],
    }
    rows["b2"][0] = "NOT_DETECTED"
    return pd.DataFrame(rows, index=loci).T


def test_distance_ignores_missing_and_counts_differing_loci():
    values = allele_matrix(_alleles())
    d = distance_matrix(values)
    assert d[0, 1] == 0.0 and d[0, 2] == 1.0  # a1 and a2 equal; every locus of b1 differs from a1
    cells = allele_matrix(pd.DataFrame({"x": ["3/4", "7", "ND", "-1"]}, index=list("abcd")))[:, 0]
    assert np.isnan(cells[[0, 2]]).all() and cells[[1, 3]].tolist() == [7.0, -1.0]
    assert np.isnan(allele_matrix(pd.DataFrame({"x": ["3/4"]}, index=["a"]))[0, 0])


def test_upgma_groups_the_close_pairs_and_keeps_every_leaf():
    alleles = _alleles()
    values = allele_matrix(alleles)
    merges = upgma(distance_matrix(values))
    assert len(merges) == len(alleles) - 1
    heights = [m[2] for m in merges]
    assert heights == sorted(heights) and heights[0] == 0.0
    order = leaf_order(merges, len(alleles))
    assert sorted(order) == list(range(len(alleles)))
    names = [alleles.index[i] for i in order]
    assert abs(names.index("a1") - names.index("a2")) == 1
    assert abs(names.index("b1") - names.index("b2")) == 1
    assert len(dendrogram_segments(merges, len(alleles), order)) == len(merges)


def test_figures_are_svg_and_small_sets_are_skipped():
    alleles = _alleles()
    svg, order = clustered_heatmap(alleles)
    assert svg.startswith("<svg") and sorted(order) == sorted(alleles.index)
    assert distance_heatmap(alleles, order).startswith("<svg")
    assert allele_spectra(alleles).startswith("<svg")
    status = pd.DataFrame("OK", index=alleles.index, columns=alleles.columns)
    status.iloc[0, 0] = "NOT_DETECTED"
    assert call_rate_bars(status).startswith("<svg")
    assert clustered_heatmap(alleles.iloc[:2]) is None and distance_heatmap(alleles, ["a1"]) is None


def test_tick_labels_stay_inside_the_canvas():
    """A short heatmap once lost the first characters of its rotated column labels."""
    import re

    alleles = _alleles()
    alleles.columns = [f"{2000 + 65 * i}kb-locus" for i in range(alleles.shape[1])]
    svg, _ = clustered_heatmap(alleles)
    height = float(re.search(r'viewBox="0 0 [\d.]+ ([\d.]+)"', svg).group(1))
    bottoms = [float(m) for m in re.findall(r'translate\([\d.]+ ([\d.]+)\) rotate\(-90\)', svg)]
    assert bottoms and max(bottoms) <= height + 0.5
