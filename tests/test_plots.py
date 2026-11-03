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


def test_minimum_spanning_tree_and_layout():
    from mlvatyper.plots import classical_mds, mds_plot, minimum_spanning_tree, mst_plot, profile_groups, \
        tree_distances, tree_layout
    alleles = _alleles()
    values = allele_matrix(alleles)
    dist = distance_matrix(values)
    edges = minimum_spanning_tree(dist)
    assert len(edges) == len(alleles) - 1
    # the pair of identical profiles is joined at zero distance
    assert any(w == 0.0 for *_, w in edges)
    # a chain a-b-c: path lengths add
    path = tree_distances(3, [(0, 1, 2.0), (1, 2, 3.0)])
    assert path[0, 2] == 5.0 and path[2, 0] == 5.0
    layout = tree_layout(3, [(0, 1, 2.0), (1, 2, 3.0)])
    assert abs(np.linalg.norm(layout[0] - layout[2]) - 5.0) < 0.2
    # identical profiles share a node
    groups = profile_groups(values, list(alleles.index))
    assert sorted(len(g) for g in groups) == [1, 1, 1, 2]
    coords, share = classical_mds(dist)
    assert coords.shape == (len(alleles), 2) and 0 < share.sum() <= 1.0 + 1e-9
    # distances along the first axes reproduce a Euclidean example exactly
    pts = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 4.0], [3.0, 4.0]])
    euclid = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    got, share = classical_mds(euclid)
    assert np.allclose(np.linalg.norm(got[:, None] - got[None], axis=2), euclid) and share.sum() > 0.999
    status = pd.DataFrame("OK", index=alleles.index, columns=alleles.columns)
    assert mst_plot(alleles, status).startswith("<svg") and mds_plot(alleles).startswith("<svg")
    assert mst_plot(alleles.iloc[:3], status.iloc[:3]) is None and mds_plot(alleles.iloc[:3]) is None


def test_the_tree_breaks_ties_by_single_locus_variants():
    """goeBURST rule: among edges of equal length the node with more SLVs wins, whatever the sample order."""
    from mlvatyper.plots import difference_counts, minimum_spanning_tree
    base = [0] * MIN_SHARED
    def profile(**changes):
        row = list(base)
        for locus, value in changes.items():
            row[int(locus[1:])] = value
        return row
    # hub H has SLVs X, Y, Z and T; T is also an SLV of X, so the edge T-H and the edge T-X are equal
    profiles = {"X": profile(l0=1), "T": profile(l0=2), "H": profile(), "Y": profile(l1=1), "Z": profile(l2=1)}
    for names in (list(profiles), list(reversed(list(profiles)))):
        values = np.array([profiles[n] for n in names], dtype=float)
        counts = difference_counts(values)
        dist = distance_matrix(values)
        edges = minimum_spanning_tree(dist, counts)
        linked = {frozenset((names[a], names[b])) for a, b, _ in edges}
        assert frozenset(("T", "H")) in linked and frozenset(("T", "X")) not in linked, names
        assert len(edges) == len(names) - 1
