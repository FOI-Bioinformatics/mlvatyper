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


def test_a_small_panel_still_gets_distances():
    """Panels of four to nine loci cannot share eight loci; the threshold is half of the panel."""
    from mlvatyper.plots import shared_threshold
    assert [shared_threshold(n) for n in (1, 4, 5, 9, 16, 42)] == [2, 2, 3, 5, 8, 8]
    rng = np.random.default_rng(1)
    alleles = pd.DataFrame(rng.integers(2, 6, size=(8, 5)).astype(str), index=[f"s{i}" for i in range(8)],
                           columns=list("abcde"))
    d = distance_matrix(allele_matrix(alleles))
    assert d.shape == (8, 8) and np.isfinite(d).all() and (d.diagonal() == 0).all()
    svg, order = clustered_heatmap(alleles)
    assert svg.startswith("<svg") and mds_plot_ok(alleles)


def mds_plot_ok(alleles):
    from mlvatyper.plots import mds_plot, mst_plot
    status = pd.DataFrame("OK", index=alleles.index, columns=alleles.columns)
    return mds_plot(alleles).startswith("<svg") and mst_plot(alleles, status) is not None


def test_no_shared_loci_gives_equal_distances_and_no_crash():
    values = np.full((6, 12), np.nan)
    for i in range(6):
        values[i, 2 * i:2 * i + 2] = 3  # every sample is typed at its own two loci only
    d = distance_matrix(values)
    assert np.isfinite(d).all() and np.allclose(d[np.triu_indices(6, 1)], d[0, 1])
    from mlvatyper.plots import classical_mds
    coords, share = classical_mds(d)
    assert np.isfinite(coords).all()


def test_loci_that_group_the_samples_alike_are_close():
    from mlvatyper.plots import locus_distance
    # A and B split the samples the same way with different repeat numbers; C splits them differently
    values = np.array([[1, 5, 1], [1, 5, 2], [2, 7, 1], [2, 7, 2], [1, 5, 1]], dtype=float)
    d = locus_distance(values)
    assert d[0, 1] == 0.0 and d[0, 2] > 0.3 and np.allclose(d, d.T) and (d.diagonal() == 0).all()
    # missing values reduce the pairs that are compared, and no pair leaves the mean distance
    values[:, 2] = np.nan
    d = locus_distance(values)
    assert np.isfinite(d).all()


def test_heatmap_orders_the_columns_by_their_clustering():
    import re
    loci = ["a", "b", "c", "d"]
    rows = {"s1": [1, 9, 3, 9], "s2": [1, 9, 4, 8], "s3": [2, 9, 3, 9], "s4": [2, 9, 4, 8], "s5": [1, 9, 3, 9]}
    alleles = pd.DataFrame({k: [str(x) for x in v] for k, v in rows.items()}, index=loci).T
    svg, _ = clustered_heatmap(alleles)
    # c and d split the samples identically, so they are next to each other among the column labels
    labels = re.findall(r'rotate\(-90\)">([a-d])</text>', svg)
    assert sorted(labels) == loci and abs(labels.index("c") - labels.index("d")) == 1


def test_half_and_censored_alleles_are_kept_in_the_matrix():
    """Bruce19 types as n.5 in nearly every strain; dropping half alleles emptied the locus."""
    alleles = pd.DataFrame({"half": ["20.5", "21.5", "18"], "censored": [">=50", "3", "4"],
                            "several": ["9/9", "2", "2"], "code": ["NOT_DETECTED", "1", "1"]},
                           index=list("abc"))
    values = allele_matrix(alleles)
    assert values[:, 0].tolist() == [20.5, 21.5, 18.0]
    assert values[:, 1].tolist() == [50.0, 3.0, 4.0]
    assert np.isnan(values[0, 2]) and np.isnan(values[0, 3])


def test_a_locus_of_half_alleles_still_clusters():
    rng = np.random.default_rng(3)
    base = rng.integers(2, 6, size=(6, 5))
    alleles = pd.DataFrame([[f"{v}.5" for v in row] for row in base], index=[f"s{i}" for i in range(6)],
                           columns=[f"m{j}" for j in range(5)])
    svg, order = clustered_heatmap(alleles)
    assert svg.startswith("<svg") and len(order) == 6
    assert "20.5" not in svg  # the cells show the alleles that are there
    assert ".5" in svg  # as half alleles, not truncated to integers
