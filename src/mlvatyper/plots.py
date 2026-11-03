"""Plots for a set of samples: clustered allele heatmap, distance matrix, allele spectra and call rates.

Clustering uses the categorical distance between profiles and average linkage (UPGMA), written
with NumPy so that the report needs no further dependency. The distance between two samples is
the share of loci that differ among the loci with a single integer allele in both; pairs that
share fewer than MIN_SHARED loci have no distance and receive the mean of the others, so that
a poorly typed sample does not pull the tree.
"""

import re

import numpy as np
import pandas as pd

NUMBER = re.compile(r"^-?\d+$")
MIN_SHARED = 8  # loci that two samples must share for a distance
ACCENT = "#1f5f8b"


def allele_matrix(alleles: pd.DataFrame) -> np.ndarray:
    """Integer alleles as floats; NaN where the cell has no single integer allele."""
    values = np.full(alleles.shape, np.nan)
    for i, (_, row) in enumerate(alleles.iterrows()):
        for j, cell in enumerate(row):
            if NUMBER.match(str(cell)):
                values[i, j] = int(str(cell))
    return values


def distance_matrix(values: np.ndarray) -> np.ndarray:
    """Share of differing loci between every pair of samples, NaN-aware."""
    n = values.shape[0]
    dist = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i, n):
            both = ~np.isnan(values[i]) & ~np.isnan(values[j])
            if both.sum() >= MIN_SHARED:
                dist[i, j] = dist[j, i] = float((values[i][both] != values[j][both]).mean())
    np.fill_diagonal(dist, 0.0)
    fill = np.nanmean(dist[~np.eye(n, dtype=bool)]) if n > 1 else 0.0
    return np.where(np.isnan(dist), 0.0 if n < 2 else fill, dist)


def upgma(dist: np.ndarray) -> list[tuple[int, int, float, int]]:
    """Average linkage merges as (left, right, height, size), clusters numbered from n upwards."""
    n = dist.shape[0]
    d = {(i, j): dist[i, j] for i in range(n) for j in range(i + 1, n)}
    size = {i: 1 for i in range(n)}
    active = set(range(n))
    merges = []
    next_id = n
    while len(active) > 1:
        (a, b), height = min(d.items(), key=lambda kv: (kv[1], kv[0]))
        new = next_id
        next_id += 1
        for c in active - {a, b}:
            da = d[(min(a, c), max(a, c))]
            db = d[(min(b, c), max(b, c))]
            d[(c, new)] = (da * size[a] + db * size[b]) / (size[a] + size[b])
        d = {k: v for k, v in d.items() if a not in k and b not in k}
        size[new] = size[a] + size[b]
        active -= {a, b}
        active.add(new)
        merges.append((a, b, float(height), size[new]))
    return merges


def leaf_order(merges: list[tuple[int, int, float, int]], n: int) -> list[int]:
    """Leaves left to right; the larger cluster of each merge goes first so that ties are stable."""
    if not merges:
        return list(range(n))
    children = {n + k: (a, b) for k, (a, b, _, _) in enumerate(merges)}
    order: list[int] = []
    stack = [n + len(merges) - 1]
    while stack:
        node = stack.pop()
        if node < n:
            order.append(node)
        else:
            a, b = children[node]
            stack.extend([b, a])
    return order


def dendrogram_segments(merges, n: int, order: list[int]):
    """Line segments (x, y) pairs for a dendrogram with leaves on the y axis, root at the right."""
    position = {leaf: i for i, leaf in enumerate(order)}
    coord = {i: (position[i], 0.0) for i in range(n)}
    segments = []
    for k, (a, b, height, _) in enumerate(merges):
        (ya, xa), (yb, xb) = coord[a], coord[b]
        segments.append(([xa, height, height, xb], [ya, ya, yb, yb]))
        coord[n + k] = ((ya + yb) / 2, height)
    return segments


def _short(svg_figure) -> str:
    import io
    import warnings

    import matplotlib.pyplot as plt
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # the dendrogram axes are switched off on purpose
        svg_figure.tight_layout()
    buffer = io.StringIO()
    # tight_layout skips the switched-off dendrogram axes and can leave tick labels outside the canvas, so the
    # canvas is also fitted to everything drawn
    svg_figure.savefig(buffer, format="svg", bbox_inches="tight", pad_inches=0.06)
    plt.close(svg_figure)
    text = buffer.getvalue()
    return text[text.index("<svg"):]


def clustered_heatmap(alleles: pd.DataFrame) -> tuple[str, list[str]] | None:
    """Samples ordered by clustering against loci; each locus is coloured on its own scale.

    Returns the SVG and the sample order, or None when there are too few samples or loci to cluster.
    """
    import matplotlib.pyplot as plt
    values = allele_matrix(alleles)
    n, m = values.shape
    if n < 3 or m < 2:
        return None
    dist = distance_matrix(values)
    merges = upgma(dist)
    order = leaf_order(merges, n)
    ordered = values[order]
    # per locus rank colours: alleles of one locus map to evenly spaced shades, so loci with a wide
    # range of repeat numbers do not hide those with a narrow one
    shade = np.full(ordered.shape, np.nan)
    for j in range(m):
        column = ordered[:, j]
        levels = np.unique(column[~np.isnan(column)])
        for k, level in enumerate(levels):
            shade[column == level, j] = k / max(len(levels) - 1, 1)
    height = max(3.0, 0.2 * n + 1.8)
    width = max(7.0, 0.3 * m + 3.5)
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, max(4, m * 0.5)], wspace=0.02)
    ax_tree = fig.add_subplot(grid[0])
    ax = fig.add_subplot(grid[1])
    for xs, ys in dendrogram_segments(merges, n, order):
        ax_tree.plot(xs, ys, color="#555555", linewidth=0.8)
    ax_tree.set_ylim(n - 0.5, -0.5)
    ax_tree.invert_xaxis()
    ax_tree.axis("off")
    cmap = plt.get_cmap("viridis").with_extremes(bad="#f0f0f0")
    ax.imshow(np.ma.masked_invalid(shade), cmap=cmap, aspect="auto", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(m))
    ax.set_xticklabels(list(alleles.columns), rotation=90)
    ax.set_yticks(range(n))
    ax.set_yticklabels([str(alleles.index[i]) for i in order])
    ax.yaxis.tick_right()
    ax.set_xticks([i - 0.5 for i in range(m + 1)], minor=True)
    ax.set_yticks([i - 0.5 for i in range(n + 1)], minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # number the alleles inside the cells when the matrix is small enough to read
    if n * m <= 900:
        for i in range(n):
            for j in range(m):
                if not np.isnan(ordered[i, j]):
                    ax.text(j, i, f"{int(ordered[i, j])}", ha="center", va="center", fontsize=5.5,
                            color="white" if shade[i, j] < 0.55 else "#222222")
    return _short(fig), [str(alleles.index[i]) for i in order]


def distance_heatmap(alleles: pd.DataFrame, order: list[str]) -> str | None:
    """Pairwise distances in the order of the clustered heatmap."""
    import matplotlib.pyplot as plt
    if len(order) < 3:
        return None
    values = allele_matrix(alleles.loc[order])
    dist = distance_matrix(values)
    n = len(order)
    fig, ax = plt.subplots(figsize=(max(4.5, 0.2 * n + 2.5), max(4.0, 0.2 * n + 2.0)))
    image = ax.imshow(dist, cmap="magma_r", vmin=0, vmax=max(float(dist.max()), 0.01), aspect="equal")
    ax.set_xticks(range(n))
    ax.set_xticklabels(order, rotation=90)
    ax.set_yticks(range(n))
    ax.set_yticklabels(order)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
    bar.set_label("share of loci that differ")
    return _short(fig)


def allele_spectra(alleles: pd.DataFrame) -> str | None:
    """For each locus, how many samples carry each allele, as small bar plots."""
    import matplotlib.pyplot as plt
    values = allele_matrix(alleles)
    m = values.shape[1]
    if m == 0 or np.isnan(values).all():
        return None
    columns = min(8, m)
    rows = -(-m // columns)
    fig, axes = plt.subplots(rows, columns, figsize=(columns * 1.35, rows * 1.15 + 0.4), squeeze=False)
    for j, marker in enumerate(alleles.columns):
        ax = axes[j // columns][j % columns]
        column = values[:, j][~np.isnan(values[:, j])]
        if len(column):
            levels, counts = np.unique(column, return_counts=True)
            ax.bar(levels, counts, color=ACCENT, width=0.8)
        ax.set_title(str(marker), fontsize=7, pad=2)
        ax.tick_params(labelsize=5.5, length=2)
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(x=0.08)
    for k in range(m, rows * columns):
        axes[k // columns][k % columns].axis("off")
    fig.supxlabel("repeat number", fontsize=7)
    fig.supylabel("samples", fontsize=7)
    return _short(fig)


def call_rate_bars(status: pd.DataFrame) -> str | None:
    """Share of OK calls per locus (top) and per sample (bottom), sorted."""
    import matplotlib.pyplot as plt
    if status.empty:
        return None
    per_locus = (status == "OK").mean(axis=0).sort_values()
    per_sample = (status == "OK").mean(axis=1).sort_values()
    n = len(per_sample)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(max(6.0, 0.28 * len(per_locus)), 2.4 + 0.1 * min(n, 60)),
                                 gridspec_kw={"height_ratios": [1, max(1.0, min(n, 60) / 20)]})
    a1.bar(range(len(per_locus)), 100 * per_locus.to_numpy(), color=ACCENT)
    a1.set_xticks(range(len(per_locus)))
    a1.set_xticklabels(per_locus.index, rotation=90)
    a1.set_ylabel("OK (% of samples)")
    a1.set_ylim(0, 100)
    a1.spines[["top", "right"]].set_visible(False)
    if n <= 60:
        a2.barh(range(n), 100 * per_sample.to_numpy(), color=ACCENT)
        a2.set_yticks(range(n))
        a2.set_yticklabels(per_sample.index)
        a2.set_xlabel("OK (% of loci)")
        a2.set_xlim(0, 100)
    else:
        a2.hist(100 * per_sample.to_numpy(), bins=20, range=(0, 100), color=ACCENT)
        a2.set_xlabel("OK (% of loci)")
        a2.set_ylabel("samples")
    a2.spines[["top", "right"]].set_visible(False)
    return _short(fig)
