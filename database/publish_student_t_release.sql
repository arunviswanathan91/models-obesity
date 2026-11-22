-- Publish the Student-t intercept collections, the model comparison and simulation part 2.
-- Run only after the Atlas branch has been reviewed. Run align_manuscript_figure_numbers_r10.sql in the same release.
-- The orthogonal patient-effect collections stay published and unchanged.
begin;
update public.atlas_web_datasets set published = true
 where cellular_model in ('student_t_intercept', 'comparison');
update public.atlas_runs set publication_status = 'published'
 where run_key like 'st\_%' or run_key like 'compare\_%' or run_key = 'simulation_part2';
-- Cell-fraction multiplicity table withdrawn from the public catalogue.
update public.atlas_web_datasets set published = false
 where id = 'reference/theta_bmi_multiplicity_recalculated';
commit;

-- Check what the public role sees.
begin;
set local role anon;
select coalesce(cellular_model, 'other') as cellular_model, count(*) as collections
from public.atlas_web_datasets group by 1 order by 1;
rollback;
