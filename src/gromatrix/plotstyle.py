"""Shared matplotlib style for every gromatrix figure (lazy: matplotlib is an optional extra)."""

from __future__ import annotations

from pathlib import Path

# Reference categorical slots (light surface): protein = slot 1, target = slot 2.
INK, MUTED, GRID, C1, C2 = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#eb6834"
OK, BAD = "#2f9e63", "#d64545"  # pass / fail marks


def have_matplotlib() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def style(ax, title: str, xlabel: str, ylabel: str | None, grid: str = "both") -> None:
    ax.set_title(title, loc="left", fontsize=10, color=INK)
    ax.set_xlabel(xlabel, fontsize=9, color=MUTED)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9, color=MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.grid(True, color=GRID, linewidth=0.6, axis=grid)
    ax.set_axisbelow(True)


def figure(ncols: int = 1, size=(6.4, 3.6), **kw):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt.subplots(1, ncols, figsize=size, dpi=200, **kw)


def save(fig, path: Path) -> str:
    import matplotlib.pyplot as plt

    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path.name
