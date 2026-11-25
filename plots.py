"""Coverage figures for the Brain & Mind datasets slide, in Open Targets style."""

from __future__ import annotations

from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import FancyBboxPatch, Rectangle

DATA = Path(__file__).parent / "data"
RUN = DATA / "ppp2612"
FIGURES = Path(__file__).parent / "figures"

OT_BLUE = "#2C77B5"
CONTEXT_GRAY = "#C5CCD4"
TEXT = "#141414"
TEXT_MUTED = "#58595B"
GRID = "#E3E6EA"
SURFACE = "#FFFFFF"

BAR_PX = 18
END_RADIUS_PX = 4

BRAIN = ("UBERON_0000955", "CL_0002319")
PSYCHIATRIC_DISORDER = "MONDO_0002025"

PGC_DISEASE_LABELS = {
    "mental disorder": "Cross-disorder p-factor",
    "internalizing disorder": "Internalizing factor",
    "attention deficit-hyperactivity disorder": "ADHD",
    "post-traumatic stress disorder": "PTSD",
    "obsessive-compulsive disorder": "OCD",
}
BRAIN_PROJECT_LABELS = {
    "GTEx_v10": "GTEx v10 (12 regions)",
    "OTAR2065_neuroinflammation": "OTAR2065 microglia",
    "Walker_2019": "Walker 2019",
    "Schwartzentruber_2018": "Schwartzentruber 2018",
    "Jerber_2021": "Jerber 2021",
    "Aygun_2021": "Aygun 2021",
    "Young_2019": "Young 2019",
}

PSYCHAD_CLASSES = {
    "EN": "Excitatory neurons",
    "IN": "Inhibitory neurons",
    "Astro": "Astrocytes",
    "Oligo": "Oligodendrocytes",
    "OPC": "OPCs",
    "Endo": "Endothelial",
    "Immune": "Immune",
    "Mural": "Mural / vascular",
}
PSYCHAD_SUBCLASS_TO_CLASS = {
    "Micro": "Immune",
    "PVM": "Immune",
    "Adaptive": "Immune",
    "PC": "Mural",
    "SMC": "Mural",
    "VLMC": "Mural",
}


def main() -> None:
    """Query the PPP 26.12 run, print the coverage numbers and write every figure."""
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "ytick.labelsize": 10,
            "text.color": TEXT,
            "axes.labelcolor": TEXT_MUTED,
            "xtick.color": TEXT_MUTED,
            "ytick.color": TEXT,
            "svg.fonttype": "none",
        }
    )
    con = duckdb.connect()
    _register_run(con)
    _print_coverage(con)
    _plot_pgc_disease_pairs(con)
    _plot_brain_qtl_projects(con)
    _plot_psychad_cell_types(con)


def _register_run(con: duckdb.DuckDBPyConnection) -> None:
    con.sql(f"""
        create view study as from '{RUN}/study/*.parquet';
        create view credible_set as select id, studyId from '{RUN}/credible_set/*.parquet';
        create view study_disease as select studyId, unnest(diseaseIds) diseaseId from study;
        create view brain_biosample as
            select biosampleId from '{RUN}/biosample/*.parquet'
            where biosampleId in {BRAIN}
               or list_has_any(ancestors, {list(BRAIN)});
        create view psychiatric_disease as
            select id diseaseId, name from '{RUN}/disease/*.parquet'
            where list_contains(therapeuticAreas, '{PSYCHIATRIC_DISORDER}');
        create view psychiatric_pair as
            select targetId, diseaseId, name,
                   bool_or(projectId = 'pgc') from_pgc, bool_or(projectId <> 'pgc') from_other
            from '{RUN}/evidence_gwas_credible_sets/*.parquet' evidence
            join credible_set on credible_set.id = evidence.studyLocusId
            join study using (studyId)
            join psychiatric_disease using (diseaseId)
            group by all;
    """)


def _print_coverage(con: duckdb.DuckDBPyConnection) -> None:
    con.sql("""
        select projectId, studyType, count(distinct studyId) studies, count(distinct id) credible_sets,
               count(distinct geneId) genes, count(distinct biosampleId) biosamples
        from study join credible_set using (studyId)
        where projectId in ('BigBrain', 'pgc')
        group by all
    """).show()
    con.sql("""
        select count(distinct geneId) all_brain_qtl_genes,
               count(distinct geneId) filter (where projectId <> 'BigBrain') without_bigbrain
        from study join credible_set using (studyId) join brain_biosample using (biosampleId)
        where studyType <> 'gwas'
    """).show()
    con.sql("""
        select count(distinct diseaseId) pgc_diseases,
               count(distinct diseaseId) filter (where studyId in (select studyId from credible_set))
                   pgc_diseases_with_credible_sets
        from study_disease join study using (studyId)
        where projectId = 'pgc'
    """).show()
    con.sql("""
        select projectId, count(distinct id) psychiatric_credible_sets
        from study join credible_set using (studyId) join study_disease using (studyId)
        join psychiatric_disease using (diseaseId)
        where studyType = 'gwas'
        group by all
    """).show()
    con.sql("""
        select count(distinct targetId) filter (where from_pgc) pgc_genes,
               count(*) filter (where from_pgc) pgc_pairs,
               count(*) filter (where from_pgc and not from_other) pgc_only_pairs,
               count(*) all_psychiatric_pairs
        from psychiatric_pair
    """).show()


def _plot_pgc_disease_pairs(con: duckdb.DuckDBPyConnection, top: int = 12) -> None:
    rows = con.sql(f"""
        select name, count(*) filter (where from_pgc) pairs,
               count(*) filter (where from_pgc and not from_other) new
        from psychiatric_pair
        group by name
        order by pairs desc
        limit {top}
    """).fetchall()
    rows.reverse()
    names = [PGC_DISEASE_LABELS.get(name, _sentence_case(name)) for name, _, _ in rows]

    fig, ax = plt.subplots(figsize=(6.0, 3.9))
    _style_value_axis(ax, len(rows), max(pairs for _, pairs, _ in rows))
    ax.set_yticks(range(len(rows)), names)
    ax.set_xlabel("PGC gene–disease pairs (L2G ≥ 0.05)")
    _legend(ax, [(OT_BLUE, "Only from PGC"), (CONTEXT_GRAY, "Also from other GWAS")])
    fig.tight_layout()
    for y, (_, pairs, new) in enumerate(rows):
        _bar(ax, y, 0, new, OT_BLUE, rounded=pairs == new)
        if pairs > new:
            _bar(ax, y, new, pairs - new, CONTEXT_GRAY, rounded=True, gap=new > 0)
        ax.text(pairs + 12, y, f"{pairs:,}", va="center", fontsize=10, color=TEXT)
    _save(fig, "pgc_disease_pairs")


def _plot_brain_qtl_projects(con: duckdb.DuckDBPyConnection) -> None:
    rows = con.sql("""
        select projectId, count(distinct id) credible_sets
        from study join credible_set using (studyId) join brain_biosample using (biosampleId)
        where studyType <> 'gwas'
        group by all
        order by credible_sets
    """).fetchall()

    fig, ax = plt.subplots(figsize=(5.0, 2.9))
    _style_value_axis(ax, len(rows), max(n for _, n in rows))
    ax.set_yticks(range(len(rows)), [BRAIN_PROJECT_LABELS.get(p, p) for p, _ in rows])
    ax.xaxis.set_major_formatter(lambda value, _: f"{value / 1000:.0f}k" if value else "0")
    ax.set_xlabel("Brain QTL credible sets")
    fig.tight_layout()
    for y, (project, n) in enumerate(rows):
        color = OT_BLUE if project == "BigBrain" else CONTEXT_GRAY
        _bar(ax, y, 0, n, color, rounded=True, bar_px=11)
        ax.text(n + 6000, y, _thousands(n), va="center", fontsize=9, color=TEXT_MUTED)
    _save(fig, "brain_qtl_projects")


def _plot_psychad_cell_types(con: duckdb.DuckDBPyConnection) -> None:
    subclasses = con.sql(f"""
        select split_part(sample_group, '_', 1) prefix, sample_size::int
        from read_csv('{DATA}/eqtl_catalogue_upcoming.tsv', delim='\t')
        where study_label = 'PsychAD' and sample_group like '%_subclass'
    """).fetchall()
    counts: dict[str, list[int]] = {cls: [] for cls in PSYCHAD_CLASSES}
    for prefix, sample_size in subclasses:
        counts[PSYCHAD_SUBCLASS_TO_CLASS.get(prefix, prefix)].append(sample_size)
    rows = sorted(counts.items(), key=lambda item: len(item[1]))

    fig, ax = plt.subplots(figsize=(5.0, 2.5))
    _style_value_axis(ax, len(rows), max(len(sizes) for _, sizes in rows))
    ax.set_yticks(range(len(rows)), [PSYCHAD_CLASSES[cls] for cls, _ in rows])
    ax.set_xlim(0, 13.5)
    ax.set_xticks([0, 3, 6, 9])
    ax.set_xlabel("Cell subtypes with an eQTL dataset")
    fig.tight_layout()
    for y, (_, sizes) in enumerate(rows):
        _bar(ax, y, 0, len(sizes), OT_BLUE, rounded=True, bar_px=12)
        donors = (
            f"{_thousands(min(sizes))}–{_thousands(max(sizes))}"
            if min(sizes) != max(sizes)
            else _thousands(sizes[0])
        )
        ax.text(len(sizes) + 0.2, y, f"{donors} donors", va="center", fontsize=9, color=TEXT_MUTED)
    _save(fig, "psychad_cell_types")


def _style_value_axis(ax: Axes, n_rows: int, max_value: float) -> None:
    ax.set_ylim(-0.6, n_rows - 0.4)
    ax.set_xlim(0, max_value * 1.12)
    ax.xaxis.grid(True, color=GRID, linewidth=1)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(GRID)


def _bar(
    ax: Axes,
    y: float,
    left: float,
    width: float,
    color: str,
    *,
    rounded: bool,
    gap: bool = False,
    bar_px: int = BAR_PX,
) -> None:
    """Draw one horizontal bar segment, square at its start and with a 4px rounded data-end."""
    if width <= 0:
        return
    fig = ax.figure
    fig.canvas.draw()
    to_data = ax.transData.inverted()
    origin = to_data.transform((0, 0))
    px = to_data.transform((1, 1)) - origin
    height = bar_px * px[1]
    gap_x = 2 * px[0] if gap else 0
    left, width = left + gap_x, width - gap_x
    radius = min(END_RADIUS_PX * px[0], width / 2) if rounded else 0
    bottom = y - height / 2
    if radius:
        ax.add_patch(
            FancyBboxPatch(
                (left, bottom),
                width,
                height,
                boxstyle=f"round,pad=0,rounding_size={radius}",
                mutation_aspect=px[1] / px[0],
                facecolor=color,
                edgecolor="none",
            )
        )
    ax.add_patch(
        Rectangle((left, bottom), width - radius, height, facecolor=color, edgecolor="none")
    )


def _legend(ax: Axes, entries: list[tuple[str, str]]) -> None:
    handles = [Rectangle((0, 0), 1, 1, facecolor=color) for color, _ in entries]
    ax.legend(
        handles,
        [label for _, label in entries],
        loc="lower right",
        frameon=False,
        fontsize=10,
        handlelength=1,
        handleheight=1,
    )


def _thousands(value: int) -> str:
    if value >= 10_000:
        return f"{value / 1000:.0f}k"
    return f"{value / 1000:.1f}k" if value >= 1000 else str(value)


def _sentence_case(name: str) -> str:
    return name[0].upper() + name[1:]


def _save(fig: plt.Figure, name: str) -> None:
    for suffix in ("svg", "png"):
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=200, transparent=suffix == "svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
