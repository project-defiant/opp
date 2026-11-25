"""Coverage figures for the Brain & Mind datasets slide, in Open Targets style."""

from __future__ import annotations

import json
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import FancyBboxPatch, Rectangle

DATA = Path(__file__).parent / "data"
RUN = DATA / "ppp2612"
FIGURES = Path(__file__).parent / "figures"
TABLES = Path(__file__).parent / "tables"

OT_BLUE = "#2C77B5"
CONTEXT_GRAY = "#C5CCD4"
TEXT = "#141414"
TEXT_MUTED = "#58595B"
GRID = "#E3E6EA"
SURFACE = "#FFFFFF"

QTL_COLORS = {"eQTL": OT_BLUE, "sQTL": "#EB6834", "tuQTL": "#1BAF7A"}
EQTL_CATALOGUE_EXCLUDED = ("BigBrain", "OTAR2065_neuroinflammation")
EQTL_CATALOGUE_LABELS = {
    "GTEx_v10": "GTEx v10 · 11 brain regions",
    "BrainSeq": "BrainSeq · DLPFC",
    "ROSMAP": "ROSMAP · DLPFC",
    "CommonMind": "CommonMind · DLPFC",
    "Walker_2019": "Walker 2019 · neocortex",
    "Schwartzentruber_2018": "Schwartzentruber 2018 · sensory neurons",
    "Jerber_2021": "Jerber 2021 · iPSC neural cells (single-cell)",
    "Braineac2": "Braineac2 · putamen, substantia nigra",
    "Aygun_2021": "Aygun 2021 · neurons",
    "Young_2019": "Young 2019 · microglia",
}
QTL_OFFSETS = {"eQTL": -0.2, "sQTL": 0.0, "tuQTL": 0.2}

BAR_PX = 18
END_RADIUS_PX = 4

BRAIN = ("UBERON_0000955", "CL_0002319")
PITUITARY = "UBERON_0000007"
PSYCHIATRIC_DISORDER = "MONDO_0002025"

OT_NAVY = "#163A5F"
FINE_MAPPING_COLORS = {"SuSiE-inf": OT_BLUE, "PICS": "#1BAF7A"}

PGC_DISEASE_LABELS = {
    "cross-disorder (multiple)": "Other cross-disorder",
    "anxiety measurement": "Anxiety (quantitative)",
    "mental disorder": "Cross-disorder p-factor",
    "internalizing disorder": "Internalizing factor",
    "attention deficit-hyperactivity disorder": "ADHD",
    "post-traumatic stress disorder": "PTSD",
    "obsessive-compulsive disorder": "OCD",
}
BRAIN_PROJECT_LABELS = {
    "GTEx_v10": "GTEx v10 (11 regions)",
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
    _plot_bigbrain_summary(con)
    _plot_eqtl_catalogue_brain(con)
    _plot_pgc_summary(con)
    _write_pgc_publications(con)
    _plot_pgc_disease_pairs(con)
    _plot_brain_qtl_projects(con)
    _plot_psychad_cell_types(con)


def _register_run(con: duckdb.DuckDBPyConnection) -> None:
    con.sql(f"""
        create view study as from '{RUN}/study/*.parquet';
        create view credible_set_method as
            select id, studyId, fineMappingMethod from '{RUN}/credible_set/*.parquet';
        create view credible_set as select id, studyId from credible_set_method;
        create view study_disease as select studyId, unnest(diseaseIds) diseaseId from study;
        create view brain_biosample as
            select biosampleId from '{RUN}/biosample/*.parquet'
            where (biosampleId in {BRAIN} or list_has_any(ancestors, {list(BRAIN)}))
              and biosampleId <> '{PITUITARY}';
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
        create view pgc_study_group as
            select studyId,
                   case when len(diseaseIds) = 1
                        then (select name from '{RUN}/disease/*.parquet' where id = diseaseIds[1])
                        else 'cross-disorder (multiple)' end disease_group
            from study where projectId = 'pgc';
        create view gwas_pair as
            select distinct targetId, diseaseId, projectId, studyId
            from '{RUN}/evidence_gwas_credible_sets/*.parquet' evidence
            join credible_set on credible_set.id = evidence.studyLocusId
            join study using (studyId);
        create view pgc_l2g as
            select disease_group, targetId, diseaseId,
                   exists (select 1 from gwas_pair other
                           where other.projectId <> 'pgc' and other.targetId = pgc.targetId
                             and other.diseaseId = pgc.diseaseId) from_other
            from gwas_pair pgc join pgc_study_group using (studyId)
            where pgc.projectId = 'pgc';
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


def _plot_bigbrain_summary(con: duckdb.DuckDBPyConnection) -> None:
    rows = con.sql("""
        select studyType, count(distinct id) credible_sets, count(distinct geneId) genes
        from study join credible_set using (studyId)
        where projectId = 'BigBrain'
        group by all
        order by studyType desc
    """).fetchall()
    labels = [f"{study_type[0]}QTL" for study_type, _, _ in rows]
    panels = [("Credible sets", [cs for _, cs, _ in rows]), ("Genes", [g for _, _, g in rows])]

    fig, axes = plt.subplots(1, 2, figsize=(6.4, 1.9), sharey=True)
    for ax, (title, values) in zip(axes, panels):
        _style_value_axis(ax, len(rows), max(values))
        ax.set_xlim(0, max(values) * 1.28)
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=TEXT)
        ax.xaxis.set_major_formatter(lambda value, _: f"{value / 1000:.0f}k" if value else "0")
    axes[0].set_yticks(range(len(rows)), labels)
    fig.tight_layout(w_pad=2.5)
    for ax, (_, values) in zip(axes, panels):
        for y, value in enumerate(values):
            _bar(ax, y, 0, value, QTL_COLORS[labels[y]], rounded=True, bar_px=20)
            ax.text(value * 1.02, y, f"{value:,}", va="center", fontsize=10, color=TEXT)
    _save(fig, "bigbrain_summary")


def _plot_eqtl_catalogue_brain(con: duckdb.DuckDBPyConnection) -> None:
    """Grouped bars of credible sets and genes per eQTL Catalogue brain project and QTL type."""
    rows = con.sql(f"""
        select projectId,
               case studyType when 'sceqtl' then 'eqtl' else studyType end study_type,
               count(distinct id) credible_sets, count(distinct geneId) genes,
               sum(count(distinct id)) over (partition by projectId) project_total
        from study join credible_set using (studyId) join brain_biosample using (biosampleId)
        where studyType <> 'gwas' and projectId not in {EQTL_CATALOGUE_EXCLUDED}
        group by projectId, study_type
        order by project_total
    """).fetchall()
    projects = list(dict.fromkeys(project for project, *_ in rows))
    values = {
        (project, "tuQTL" if study_type == "tuqtl" else f"{study_type[0]}QTL"): (credible_sets, genes)
        for project, study_type, credible_sets, genes, _ in rows
    }

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 6.4), sharey=True)
    for ax, (title, column) in zip(axes, [("Credible sets", 0), ("Genes", 1)]):
        maximum = max(v[column] for v in values.values())
        _style_value_axis(ax, len(projects), maximum)
        ax.set_xlim(0, maximum * 1.18)
        ax.xaxis.set_major_formatter(lambda value, _: f"{value / 1000:.0f}k" if value else "0")
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=TEXT)
    axes[0].set_yticks(range(len(projects)), [EQTL_CATALOGUE_LABELS[p] for p in projects])
    _legend(axes[1], [(QTL_COLORS[q], q) for q in QTL_COLORS])
    fig.tight_layout(w_pad=2.5)
    for ax, column in zip(axes, (0, 1)):
        for y, project in enumerate(projects):
            for qtl in QTL_COLORS:
                if (project, qtl) not in values:
                    continue
                value = values[project, qtl][column]
                bar_y = y + QTL_OFFSETS[qtl] * -1.3
                _bar(ax, bar_y, 0, value, QTL_COLORS[qtl], rounded=True, bar_px=10)
                ax.text(
                    value + ax.get_xlim()[1] * 0.01, bar_y, _thousands(value),
                    va="center", fontsize=7.5, color=TEXT_MUTED,
                )
    _save(fig, "eqtl_catalogue_brain")


def _plot_pgc_summary(con: duckdb.DuckDBPyConnection, top: int = 12) -> None:
    """Credible sets by fine-mapping method and L2G genes by novelty, per PGC disease group."""
    groups = con.sql("""
        with credible_sets as (
            select disease_group,
                   count(distinct id) filter (where fineMappingMethod = 'SuSiE-inf') susie,
                   count(distinct id) filter (where fineMappingMethod = 'PICS') pics
            from pgc_study_group
            left join (select id, studyId, fineMappingMethod from credible_set_method) using (studyId)
            group by all
        ),
        genes as (
            select disease_group, count(distinct targetId) genes,
                   count(distinct targetId) filter (where not from_other) new_genes
            from pgc_l2g group by all
        )
        select disease_group, susie, pics, coalesce(genes, 0), coalesce(new_genes, 0)
        from credible_sets left join genes using (disease_group)
        order by susie + pics desc, disease_group
    """).fetchall()
    shown, rest = groups[:top], groups[top:]
    other = (f"Other ({len(rest)} traits)", *[sum(row[i] for row in rest) for i in range(1, 5)])
    rows = [other, *reversed(shown)]
    labels = [PGC_DISEASE_LABELS.get(row[0], _sentence_case(row[0])) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.6), sharey=True)
    panels = [
        ("Credible sets", [(r[1], r[2]) for r in rows], [OT_BLUE, FINE_MAPPING_COLORS["PICS"]],
         ["SuSiE-inf", "PICS"]),
        ("L2G genes", [(r[4], r[3] - r[4]) for r in rows], [OT_NAVY, CONTEXT_GRAY],
         ["New gene–disease pair", "Also in other GWAS"]),
    ]
    for ax, (title, values, _, _) in zip(axes, panels):
        _style_value_axis(ax, len(rows), max(sum(v) for v in values))
        ax.set_xlim(0, max(sum(v) for v in values) * 1.2)
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", color=TEXT)
    axes[0].set_yticks(range(len(rows)), labels)
    fig.tight_layout(rect=(0, 0.1, 1, 1), w_pad=2.5)
    for ax, (_, values, colors, legend) in zip(axes, panels):
        for y, (first, second) in enumerate(values):
            _bar(ax, y, 0, first, colors[0], rounded=second == 0, bar_px=14)
            _bar(ax, y, first, second, colors[1], rounded=True, gap=first > 0, bar_px=14)
            total = first + second
            ax.text(total + ax.get_xlim()[1] * 0.015, y, f"{total:,}", va="center", fontsize=9, color=TEXT)
        handles = [Rectangle((0, 0), 1, 1, facecolor=c) for c in colors]
        ax.legend(handles, legend, loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=2,
                  frameon=False, fontsize=9, handlelength=1, handleheight=1)
    _save(fig, "pgc_summary")


def _write_pgc_publications(con: duckdb.DuckDBPyConnection) -> None:
    """Join the cached PubMed summaries to PGC study and credible set counts."""
    result = json.loads((DATA / "pgc_pubmed.json").read_text())["result"]
    rows = []
    for pmid in result["uids"]:
        record = result[pmid]
        doi = next((a["value"] for a in record["articleids"] if a["idtype"] == "doi"), None)
        author = record["sortfirstauthor"].split(". Electronic address")[0]
        rows.append((pmid, author, record["pubdate"][:4], record["source"],
                     record["title"], doi))
    con.execute(
        "create or replace table pubmed "
        "(pmid varchar, first_author varchar, year varchar, journal varchar, title varchar, doi varchar)"
    )
    con.executemany("insert into pubmed values (?, ?, ?, ?, ?, ?)", rows)
    TABLES.mkdir(exist_ok=True)
    con.sql("""
        select coalesce(pmid, '(no PMID)') pmid, first_author, year, journal, title, doi,
               count(distinct studyId) studies, count(distinct id) credible_sets,
               string_agg(distinct disease_group, '; ') diseases
        from study
        join pgc_study_group using (studyId)
        left join credible_set using (studyId)
        left join pubmed on pmid = pubmedId
        where projectId = 'pgc'
        group by all
        order by credible_sets desc, studies desc
    """).write_csv(str(TABLES / "pgc_publications.tsv"), sep="\t")


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
