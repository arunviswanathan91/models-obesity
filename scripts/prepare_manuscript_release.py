"""Prepare a verified, aggregate-only manuscript release. No model fitting.
Usage: python prepare_manuscript_release.py FIGURE_ZIP ASAH1_CSV PAIRED_ASAH1_CSV OUTPUT_DIR
SQL batches stage unpublished datasets; publish only after count/value checks.
"""
import csv,io,json,hashlib,re,sys,zipfile
from pathlib import Path
archive,asah,paired,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
web=Path(__file__).resolve().parents[1]/'web'
z=zipfile.ZipFile(archive); run='manuscript_ptm_20261004';collections=[];figures={};manifest=[]
def typed(v):
 if v=='':return None
 if v in ('True','False'):return v=='True'
 if re.fullmatch(r'-?\d+',v):return int(v)
 try:return float(v)
 except ValueError:return v

def add(path,key,title,description,kind='table',section='manuscript',raw=None,transform=None):
 raw=raw if raw is not None else z.read(path); original=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))));rows=[{k:typed(v) for k,v in r.items()} for r in original]
 if transform:rows=transform(rows)
 item=dict(id=section+'/'+key,title=title,section=section,kind=kind,description=description,run_key=run,row_count=len(rows),provisional=False,published=False)
 collections.append((item,original,rows));manifest.append(dict(id=item['id'],source_file=path,sha256=hashlib.sha256(raw).hexdigest(),rows=len(rows),source_rows=len(original)))
 return item['id']
def effects(rows):
 for r in rows:
  r.update(cell_type=None,feature=r['feature_id'],label=r.get('display_label') or r['gene']+' | '+r['feature_id'],component=r['quantity'],contrast=r['term'],estimate=r['posterior_mean'],hdi_lower=r['hdi_95_lower'],hdi_upper=r['hdi_95_upper'],units=r['effect_units'],entity_level='feature',interval_excludes_zero=r['hdi_excludes_zero'],interval_probability=.95)
  r.pop('source',None)
 return rows
ids={}
for key,path,title,desc in [
 ('ptm_heatmap_lookup','selection/PTM_full_label_lookup.csv','Figure 6 · PTM heatmap identities and estimates','All 14 phospho and top 20 of 67 glyco conditional candidates. Gene, site/peptide/glycan, both effects and pointwise HDIs; magnitude-ranked without multiplicity control.'),
 ('phospho_labels','selection/ptm_phospho_label_map.csv','Figure 6c · phosphosite labels','The 14 displayed phosphosites with complete feature identities and both effect columns.'),
 ('glyco_labels','selection/ptm_glyco_label_map.csv','Figure 6d · glycofeature labels','The 20 displayed glycofeatures with complete peptide and glycan identities and both effect columns.'),
 ('rna_protein_labels','selection/rna_protein_label_map.csv','Figure 5e · RNA–protein programme labels','RP01–RP16 link displayed labels to unique measured programmes and source signatures.'),
 ('rna_heatmap_values','heatmaps/paper_rna_top_heatmap_displayed_values.csv','Figure 3 · RNA heatmap values','Displayed continuous RNA coefficients for the 24 selected signature rows. Unmodeled combinations remain missing.'),
 ('rna_protein_heatmap_values','heatmaps/paper_rp_top_heatmap_displayed_values.csv','Figure 5e · RNA–protein heatmap values','Displayed paired-model RNA and protein slopes and direct protein-minus-RNA contrasts.'),
 ('candidate_sensitivity_coverage','updated_ptm_tables/candidate_sensitivity_availability.csv','Figure S18 · candidate sensitivity coverage','Covariate-sensitivity availability for 27 of 81 primary conditional candidates: 24 selected features plus three ASAH1 glycoforms.'),
 ('sensitivity_feature_labels','PTM_extended/tables/feature_identity_map.csv','Figure S25 · sensitivity feature identities','E01–E24 map to the 24 post hoc selected features from CTSD, HLA-DRA, LAMC1, POSTN and GSDMD.'),
 ('sensitivity_comparisons','PTM_extended/tables/matched_conditional_comparisons.csv','Figure S25 · matched sensitivity comparisons','192 comparisons. Stage and purity adjustment use the matched unadjusted subset; restriction and other variants use the primary fit. Classification changes are descriptive, not posterior contrast tests.'),
 ('sensitivity_classifications','PTM_extended/tables/interval_classification_summary.csv','Figure S25b · interval classifications','Counts retained or changed in each comparison among 24 post hoc features. Both compared fits must be numerically adequate.'),
 ('ptm_calibration_summary','PTM_extended/tables/calibration_summary_recomputed.csv','Figure S24 · PTM calibration summary','Six quantities, five scenarios and two scopes: all completed versus numerically adequate fits. Coverage bars are Wilson confidence intervals, not HDIs. No screen-wide false-discovery control.'),
 ('ptm_calibration_templates','PTM_extended/tables/calibration_by_template.csv','Figure S24 · calibration by template','Separate phospho and glyco template results; aggregate summaries across simulated replicates.'),
 ('ptm_calibration_replicates','PTM_extended/tables/calibration_replicate_metrics.csv','Figure S24 · per-replicate calibration metrics','250 simulated datasets × six quantities. Literal null denotes the null scenario; these are simulated quantities, not patient records.'),
 ('ptm_calibration_completion','PTM_extended/tables/calibration_completion.csv','Figure S24e · calibration sample counts','250 completed simulated datasets; 249 numerically adequate.'),
 ('ptm_sbc_bins','PTM_extended/tables/sbc_rank_bins.csv','Figure S24d · SBC rank histogram counts','Pointwise 95% binomial reference limits are not simultaneous bands or a formal calibration test.'),
 ('ptm_sbc_values','PTM_extended/tables/sbc_rank_values.csv','Figure S24d · SBC ranks','Prior-predictive SBC ranks for 50 simulated datasets × six quantities. Randomized ties use seed 20261004.')]:ids[key]=add(path,key,title,desc)
ids['asah1_labels']=add(asah.name,'asah1_labels','Figure S19b · ASAH1 glycofeature identities','A01–A20 identify distinct peptide-plus-glycan features. All 20 conditional estimates are retained.',raw=asah.read_bytes())
ids['asah1_paired']=add(paired.name,'asah1_paired','ASAH1 · paired model quantities','Primary protein, marginal PTM and conditional PTM effects, standardized slope differences, coupling and residual correlation. All 20 parent-protein BMI intervals include zero.',kind='effects',section='ptm',raw=paired.read_bytes(),transform=effects)
ids['sensitivity_effects']=add('PTM_extended/tables/all_effects_95HDI.csv','extended_sensitivity_effects','Figure S25 · expanded PTM sensitivity estimates','24 post hoc selected features, primary references and eight sensitivity variants. 192 refits; all numerically adequate. Use matched subsets for stage and purity adjustment.',kind='effects',section='ptm',transform=effects)
ids['sensitivity_predictive']=add('PTM_extended/tables/sensitivity_predictive_checks.csv','extended_sensitivity_predictive','Expanded PTM sensitivity · predictive checks','Saved predictive checks for the 192 refits. Sampling adequacy and predictive fit assess different aspects of a model.',kind='diagnostics',section='ptm')
# Literal null must not be converted to a missing value.
assert {r['scenario'] for m,_,rs in collections if m['id']==ids['ptm_calibration_summary'] for r in rs}=={'null','protein_only','conditional_signal','heavy_tails','sbc'}
(web/'manuscript-manifest.json').write_text(json.dumps({'release':run,'source_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'datasets':manifest},indent=2))
(out/'collections.json').write_text(json.dumps([dict(**m,rows=rows) for m,_,rows in collections]))
q=lambda x:"'"+str(x).replace("'","''")+"'"
j=lambda x:q(json.dumps(x,allow_nan=False,separators=(',',':')))+'::jsonb'
queries=[]
queries.append('insert into atlas_runs(run_key,modality,exposure_model,source_manifest,publication_status) values('+q(run)+",'ptm','continuous',"+j({'source_archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'datasets':manifest})+",'staged') on conflict(run_key) do nothing;")
for (m,rawrows,rows),source in zip(collections,manifest):
 sk=m['id']; fields=list(m)
 queries.append('begin; insert into atlas_sources(source_key,run_key,source_path,sha256,expected_rows,columns_json,import_status) values('+','.join([q(sk),q(run),q(source['source_file']),q(source['sha256']),str(len(rawrows)),j(list(rawrows[0])),"'verified'"])+') on conflict(source_key) do nothing; insert into atlas_web_datasets('+','.join(fields)+') values('+','.join(j(v) if isinstance(v,(dict,list)) else str(v).lower() if isinstance(v,(bool,int)) else q(v) for v in m.values())+') on conflict(id) do nothing; commit;')
 # Chunk below tool transport limits. Each transaction retains raw and public projections together.
 for start in range(0,len(rows),100):
  a=rawrows[start:start+100];b=rows[start:start+100]
  queries.append('begin; insert into atlas_records(source_key,row_number,payload) values '+','.join('('+q(sk)+','+str(start+i+1)+','+j(r)+')' for i,r in enumerate(a))+' on conflict(source_key,row_number) do nothing; insert into atlas_web_rows(dataset_id,row_number,data) values '+','.join('('+q(sk)+','+str(start+i+1)+','+j(r)+')' for i,r in enumerate(b))+' on conflict(dataset_id,row_number) do nothing; commit;')
(out/'queries.json').write_text(json.dumps(queries))
print(json.dumps({'collections':len(collections),'rows':sum(m['row_count'] for m,_,_ in collections),'queries':len(queries),'ids':ids},indent=2))
