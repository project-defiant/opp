-- Build opp.duckdb with the PGC slice of the PPP 26.12 run (ppp-2612-1).
-- Run from the repository root: duckdb opp.duckdb < sql/load_pgc.sql

create or replace table study as
select *
from 'data/ppp2612/study/*.parquet'
where projectId = 'pgc';

create or replace table credible_set as
select id, studyId, chromosome, leadVariantId, fineMappingMethod, qualityControls
from 'data/ppp2612/credible_set/*.parquet'
where studyId in (select studyId from study);

create or replace table disease as
select id, name, therapeuticAreas
from 'data/ppp2612/disease/*.parquet';

-- One row per GWAS gene–disease pair (L2G ≥ 0.05) and source project, across all
-- GWAS in the run, so PGC pairs can be compared with the other sources.
create or replace table gwas_pair as
select distinct evidence.targetId, evidence.diseaseId, all_study.projectId, all_study.studyId
from 'data/ppp2612/evidence_gwas_credible_sets/*.parquet' evidence
join 'data/ppp2612/credible_set/*.parquet' all_credible_set
    on all_credible_set.id = evidence.studyLocusId
join 'data/ppp2612/study/*.parquet' all_study
    on all_study.studyId = all_credible_set.studyId;

-- GWAS credible sets mapped to a psychiatric disease, for every GWAS source.
create or replace table psychiatric_gwas_credible_set as
select distinct all_study.projectId, all_credible_set.id, study_disease.diseaseId
from 'data/ppp2612/study/*.parquet' all_study
join 'data/ppp2612/credible_set/*.parquet' all_credible_set using (studyId)
join (
    select studyId, unnest(diseaseIds) diseaseId from 'data/ppp2612/study/*.parquet'
) study_disease using (studyId)
join disease on disease.id = study_disease.diseaseId
where all_study.studyType = 'gwas'
  and list_contains(disease.therapeuticAreas, 'MONDO_0002025');  -- psychiatric disorder
