"""Convert saved aggregate deconvolution QC tables to public rows; verify purity against crosstabs.
Usage: python prepare_purity_release.py SOURCES_JSON OUTPUT_DIR
Source entries: compartment, title, text. No fitting or patient-level export.
"""
import csv,hashlib,io,json,math,sys
from pathlib import Path
sources=json.load(open(sys.argv[1]));out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
def convert(v):
 try:return float(v) if any(c in v for c in '.eE') else int(v)
 except:return v
rows={}
for f in sources:
 rs=[{k:convert(v) for k,v in r.items()} for r in csv.DictReader(io.StringIO(f['text']))]
 rows.setdefault(f['title'],[]).extend({'compartment':f['compartment'],**r} for r in rs)
for r in rows['celltype_purity.csv']:
 x=next(x for x in rows['confusion_crosstab.csv'] if x['compartment']==r['compartment'] and x['CellType']==r['CellType']);counts=[v for k,v in x.items() if k.isdigit()];assert math.isclose(max(counts)/sum(counts),r['max_cluster_purity'],abs_tol=1e-12);r['n_profiles']=sum(counts)
meta={
'celltype_purity.csv':('Deconvolution cell-type clustering purity','purity','Saved within-cell-type concentration in the dominant K-means cluster: max cluster count / total profiles for that type. Verified against the saved crosstab. This is not tumour purity, expression contamination, or independent accuracy validation. The purity CSV does not record raw versus row-normalized mode; no mode is inferred.'),
'confusion_crosstab.csv':('Deconvolution cluster assignment counts','summary','Saved cell-type by K-means-cluster counts. Cluster numbers are arbitrary within each compartment and cannot be matched between compartments. Types can share a cluster even when their individual purity equals 1.'),
'separation_metrics.csv':('Deconvolution PCA separation metrics','summary','Saved silhouette and adjusted Rand index for raw and row-normalized expression geometry. Separate from the later signature-score PCA. Raw-PC1 correlation with total signal is retained as a magnitude diagnostic.'),
'per_celltype_patient_spread.csv':('Deconvolution within-type PC spread','summary','Saved standard deviations of PC1, PC2 and PC3 across patient profiles within each cell type. PC1–PC3 columns are spreads, not PCA coordinates. This archive table does not specify raw versus row-normalized mode.')}
q=lambda x:"'"+str(x).replace("'","''")+"'"
cs=[];manifest=[]
for filename,(title,kind,desc) in meta.items():
 rs=rows[filename];c=dict(id='reference/'+filename[:-4],section='reference',title=title,kind=kind,description=desc,run_key='deconvolution_pca_qc_archive',row_count=len(rs),provisional=False,rows=rs);cs.append(c)
 fields=['id','section','title','kind','description','run_key','row_count','provisional'];sql='begin;insert into atlas_web_datasets('+','.join(fields)+',published) values('+','.join(q(c[k]) for k in fields)+',false) on conflict(id) do nothing;'
 sql+='insert into atlas_web_rows(dataset_id,row_number,data) values'+','.join('('+q(c['id'])+','+str(i+1)+','+q(json.dumps(r))+'::jsonb)' for i,r in enumerate(rs))+' on conflict(dataset_id,row_number) do nothing;commit;'
 (out/(filename[:-4]+'.sql')).write_text(sql)
for f in sources:manifest.append(dict(collection_id='reference/'+f['title'][:-4],source_file=f['compartment']+'/'+f['title'],sha256=hashlib.sha256(f['text'].encode()).hexdigest()))
(out/'collections.json').write_text(json.dumps(cs));(out/'manifest.json').write_text(json.dumps(manifest,indent=2));print([(c['id'],c['row_count']) for c in cs])
