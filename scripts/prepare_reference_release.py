"""Build aggregate-only release SQL from verified Drive exports. No model fitting.
Usage: python prepare_reference_release.py INVENTORY_JSON ANNOTATION_JSON OUTPUT_DIR
INVENTORY_JSON entries contain title and local.path; annotation entries contain title/text.
The input manifests are private working files and are never copied into the website.
"""
import csv,hashlib,io,json,re,sys
from pathlib import Path
import pandas as pd
inventory,annotations,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
meta={
'pca_pc_bmi_associations':('PC–BMI associations','summary','Descriptive saved PC–BMI associations with HC3 standard errors and within-compartment/set q-values. This exploratory PCA QC is not a primary BMI model.'),
'pca_component_contributions':('PCA cell-type contributions','contributions','Saved squared-axis contributions from deconvolved expression, grouped by cell type. Select one compartment, QC set and PC.'),
'pca_explained_variance':('PCA explained variance','scree','Saved PCA of within-component logCPM z-scores capped at 5. Equal total squared weight per retained correlation cluster and cell type. Tiered and strict QC sets are separate.'),
'pca_tiered_vs_strict_stability':('PCA QC sensitivity and review flags','summary','Tiered versus strict QC geometry. Non-immune results retain the recommendation to review strict QC as a possible primary choice; this is not a resolved validation result.'),
'pca_input_counts':('PCA input counts','summary','Saved cohort, feature and cell-type counts. Broad primary family: coarse immune plus non-immune. Fine immune is a nested secondary family.'),
 'theta_bmi_multiplicity_recalculated':('Deconvolution fractions and BMI','summary','Aggregate BayesPrism theta summaries and saved BMI associations. Original within-compartment q-values and recalculated global q-values are all retained; these adjustment families are not interchangeable.'),
'component_inclusion_tiers':('Deconvolution component inclusion and QC','summary','All component QC flags and inclusion tiers. Tiered QC excludes zero-total or median zero-gene-fraction >0.20 components; strict QC additionally excludes low-theta/high-theta-CV components.'),
'normalized_dependence_summary':('Normalized score dependence','summary','Saved dependence summaries across normalization variants and component scopes; all review outcomes retained.'),
'normalized_dependence_by_celltype':('Score dependence by cell type','summary','Saved high-correlation counts before and after residual adjustment. Select a normalization variant and compartment.'),
'raw_vs_normalized_summary':('Raw versus normalized scores','summary','Saved correlations and descriptive slope changes between raw and normalized endpoints. These are endpoint checks, not new fits.'),
'cap3_vs_cap5_summary':('Score clipping sensitivity: cap 3 versus cap 5','summary','Saved score correlations and descriptive slope changes under clipping at 3 or 5.'),
'component_reliability_all':('Deconvolution component reliability','summary','Full aggregate reliability table, including low-theta, high-CV and sparse-expression flags.')}
collections=[];manifest=[]
for f in json.loads(inventory.read_text()):
 p=Path(f['local']['path']);raw=p.read_bytes();sha=hashlib.sha256(raw).hexdigest()
 if p.suffix=='.json':
  assert sha=='3a6cc8e21208e1f6fe0af605da077229e7647d935ddb38a102351e62925ed456', 'Signature catalogue differs from audited version'
  d=json.loads(raw);rows=[dict(cell_type=c,signature=s,gene_count=len(g),genes=g) for c,v in d.items() for s,g in v.items()]
  assert all(isinstance(g,str) for r in rows for g in r['genes'])
  section='signatures';key='audited_gene_sets';kind='signatures';title='Audited gene signature catalogue'
  desc='2,242 gene sets under 65 original cell labels. Exact catalogue hash matches the September 2026 analysis audit. Gene sets score deconvolved expression; this is not the BayesPrism single-cell reference matrix. Original labels and aliases are retained; this catalogue is not a list of statistically significant findings.'
 else:
  key=p.stem;title,kind,desc=meta[key];section='reference'
  df=pd.read_csv(p);rows=json.loads(df.to_json(orient='records',double_precision=15))
  assert not any(re.search(r'(^|_)(sample_id|patient_id|barcode|path|filename)($|_)',c,re.I) for c in df.columns)
 payload=json.dumps(rows);assert not re.search(r'C3[NL]-\d|/content/|MyDrive|@',payload)
 collections.append(dict(id=section+'/'+key,section=section,title=title,kind=kind,description=desc,run_key='reference_audit_20260906',row_count=len(rows),provisional=False,rows=rows))
 manifest.append(dict(collection_id=section+'/'+key,source_file=f['title'],sha256=sha,rows=len(rows)))
ann=[]
for f in json.loads(annotations.read_text()):
 m=re.fullmatch(r'annotation_summary_(\w+)_(fine|main)\.csv',f['title']);assert m
 rows=list(csv.DictReader(io.StringIO(f['text'])));assert rows and set(rows[0])=={'Cell_Type','Count'}
 ann.extend(dict(reference=m[1],resolution=m[2],annotation=r['Cell_Type'],cell_count=int(r['Count'])) for r in rows)
 manifest.append(dict(collection_id='single_cell/annotation_counts',source_file=f['title'],sha256=hashlib.sha256(f['text'].encode()).hexdigest(),rows=len(rows)))
collections.append(dict(id='single_cell/annotation_counts',section='single_cell',title='Single-cell reference annotation counts',kind='annotation',description='Archived SingleR reference-based annotation summaries from the BayesPrism annotation outputs: Blueprint, Novershtern, Monaco and DICE, each with main and fine labels. These are alternative reference annotations, not independent cohorts or validated final labels. Do not sum across panels. No barcodes or patient metadata are published.',run_key='single_cell_reference_annotation_archive',row_count=len(ann),provisional=False,rows=ann))
def q(x):return "'"+str(x).replace("'","''")+"'"
for i,c in enumerate(collections):
 fields=['id','section','title','kind','description','run_key','row_count','provisional']
 sql='begin;\ninsert into public.atlas_web_datasets('+','.join(fields)+',published) values('+','.join(q(c[k]) for k in fields)+',false) on conflict(id) do nothing;\n'
 sql+='insert into public.atlas_web_rows(dataset_id,row_number,data) values\n'+',\n'.join('('+q(c['id'])+','+str(n+1)+','+q(json.dumps(r,allow_nan=False))+'::jsonb)' for n,r in enumerate(c['rows']))+' on conflict(dataset_id,row_number) do nothing;\ncommit;'
 (out/f'{i:02}.sql').write_text(sql)
(out/'collections.json').write_text(json.dumps(collections))
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'collections':len(collections),'rows':sum(c['row_count'] for c in collections),'counts':{c['id']:c['row_count'] for c in collections}},indent=2))
