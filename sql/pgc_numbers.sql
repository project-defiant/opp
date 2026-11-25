-- PGC numbers from the PPP 26.12 run. Run after sql/load_pgc.sql:
-- duckdb opp.duckdb < sql/pgc_numbers.sql

-- Studies, and how many have credible sets: 192 / 107
select count(*) studies,
       count(*) filter (where studyId in (select studyId from credible_set)) studies_with_credible_sets
from study;

-- Credible sets by fine-mapping method: 1,605 SuSiE-inf + 2,249 PICS = 3,854
select fineMappingMethod, count(*) credible_sets
from credible_set
group by all
union all
select 'total', count(*) from credible_set;

-- Diseases: 47 mapped, 27 with credible sets
select count(distinct diseaseId) diseases,
       count(distinct diseaseId) filter (where studyId in (select studyId from credible_set))
           diseases_with_credible_sets
from (select studyId, unnest(diseaseIds) diseaseId from study);

-- Publications: 51 PMIDs, 11 studies without one
select count(distinct pubmedId) publications,
       count(*) filter (where pubmedId is null) studies_without_pmid
from study;

-- Ancestry, counted by study (EUR 155, AFR 30, EAS 29, AMR 12)
select ld.ldPopulation, count(*) studies
from study, unnest(ldPopulationStructure) t(ld)
group by all
order by studies desc;

-- Sample size: min 527, median 70,876, max 2,598,930
select min(nSamples), median(nSamples), max(nSamples) from study;

-- Psychiatric GWAS credible sets by source: GCST 8,752, pgc 3,768, FINNGEN_R12 723
select projectId, count(distinct id) credible_sets
from psychiatric_gwas_credible_set
group by all
order by credible_sets desc;

-- PGC L2G on psychiatric diseases: 1,998 genes, 5,034 pairs, 1,933 found by no other GWAS
with psychiatric_pair as (
    select targetId, diseaseId,
           bool_or(projectId = 'pgc') from_pgc,
           bool_or(projectId <> 'pgc') from_other
    from gwas_pair
    join disease on disease.id = gwas_pair.diseaseId
    where list_contains(disease.therapeuticAreas, 'MONDO_0002025')
    group by all
)
select count(distinct targetId) filter (where from_pgc) pgc_genes,
       count(*) filter (where from_pgc) pgc_pairs,
       count(*) filter (where from_pgc and not from_other) pgc_only_pairs
from psychiatric_pair;

-- Credible sets and L2G genes per disease group (the pgc_summary plot).
-- Studies mapped to several diseases form the 'cross-disorder (multiple)' group.
with study_group as (
    select studyId,
           case when len(diseaseIds) = 1
                then (select name from disease where id = diseaseIds[1])
                else 'cross-disorder (multiple)' end disease_group
    from study
),
credible_sets as (
    select disease_group,
           count(distinct id) filter (where fineMappingMethod = 'SuSiE-inf') susie,
           count(distinct id) filter (where fineMappingMethod = 'PICS') pics
    from study_group left join credible_set using (studyId)
    group by all
),
genes as (
    select disease_group,
           count(distinct pgc.targetId) l2g_genes,
           count(distinct pgc.targetId) filter (where not exists (
               select 1 from gwas_pair other
               where other.projectId <> 'pgc'
                 and other.targetId = pgc.targetId
                 and other.diseaseId = pgc.diseaseId
           )) genes_with_new_pair
    from gwas_pair pgc join study_group using (studyId)
    where pgc.projectId = 'pgc'
    group by all
)
select disease_group, susie, pics, susie + pics credible_sets,
       coalesce(l2g_genes, 0) l2g_genes, coalesce(genes_with_new_pair, 0) genes_with_new_pair
from credible_sets left join genes using (disease_group)
order by credible_sets desc, disease_group;
