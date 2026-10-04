import {CONFIG} from './config.js';
import {PROGRAMME_LABELS} from './programme-labels.js';
const $=id=>document.getElementById(id);
const sections=[['rna','RNA programmes'],['rna_protein','RNA–protein'],['ptm','Protein & PTM'],['manuscript','Manuscript data'],['simulation','Simulations'],['reference','Deconvolution & PCA'],['signatures','Gene signatures'],['single_cell','Single-cell reference']];
const state={catalog:[],section:'rna',dataset:null,rows:[],filtered:[],facets:{},page:0,view:'plot',mode:'biology',version:0,figure:null,worker:null,plotId:0};
const cache=new Map();
const number=new Intl.NumberFormat('en',{maximumSignificantDigits:5});
const human=s=>String(s).replaceAll('_',' ');
function programmeLabel(r){const key=r.feature||r.programme;const names=PROGRAMME_LABELS[key];return names?names.map(n=>human(n.replace(/_Signature$/,''))).join(' / ')+' ['+key.slice(-6)+']':r.label||r.feature||r.programme;}
function searchable(r){return Object.values(r).concat(programmeLabel(r)||'');}
function el(tag,text,attrs={}){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);return e;}
function status(message,error=false){$('status').textContent=message;$('status').classList.toggle('error',error);}
async function api(table,query){const url=new URL(CONFIG.url+'/rest/v1/'+table);for(const[k,v]of Object.entries(query))url.searchParams.set(k,v);const res=await fetch(url,{headers:{apikey:CONFIG.key},signal:AbortSignal.timeout(30000)});if(!res.ok)throw Error('Results could not be loaded ('+res.status+'). Please try again.');return res.json();}
function options(select,values,selected){select.replaceChildren();for(const [value,label]of values){const o=el('option',label,{value});select.append(o);}if(selected!==undefined)select.value=selected;}
function safeValue(v){if(v===null||v===undefined||v==='')return '—';if(typeof v==='number')return number.format(v);if(typeof v==='object')return JSON.stringify(v);return String(v);}
function requestDownload(){ $('download-contact').showModal(); }
function save(){ requestDownload(); }
function csvCell(v){let s=v==null?'':typeof v==='object'?JSON.stringify(v):String(v);if(/^[=+@\t\r]/.test(s)||(/^-.+/.test(s)&&!Number.isFinite(Number(s))))s="'"+s;return '"'+s.replaceAll('"','""')+'"';}
function filteredCSV(rows){const cols=[...new Set(rows.flatMap(Object.keys))];return [cols.map(csvCell).join(','),...rows.map(r=>cols.map(c=>csvCell(r[c])).join(','))].join('\r\n');}
function clearFigure(){state.plotId++;pendingPlot=null;$('cancel-plot').hidden=true;state.figure=null;$('svg').disabled=$('png').disabled=true;$('plot-caption').textContent='';$('figure-container').replaceChildren(el('div','Choose a figure type and draw the current selection.',{class:'plot-empty'}));}
const runLabel={immune_coarse:'Coarse immune',immune_fine:'Fine immune',nonimmune:'Non-immune'};

function brief(d){
 if(d.run_key==='manuscript_ptm_20261004')return d.description;
 if(d.kind==='signatures')return 'Gene sets used to score biological programmes within cell types.';
 if(d.kind==='pca_original')return 'Original PCA figures for three deconvolution compartments.';
 if(d.kind==='diagnostics')return 'Sampling precision, chain agreement and predictive checks.';
 if(d.kind==='effects')return 'BMI associations with pointwise 95% posterior intervals.';
 return explain(d);
}
function explain(d){
 if(d.run_key==='manuscript_ptm_20261004')return d.description;
 if(d.kind==='signatures')return 'A signature is a named group of genes representing a biological process. Browse its member genes and compare the number of genes across signatures within a cell type.';
 if(d.kind==='pca_original')return 'Principal component analysis (PCA) places similar expression profiles near one another. Choose a compartment to inspect separation among its cell types.';
 if(d.kind==='effects')return 'The estimate describes the direction and size of the BMI association. The 95% highest-density interval (HDI) describes posterior uncertainty. Read the measurement scale and BMI contrast before comparing effects; exclusion of zero alone does not control multiple testing.';
 if(d.kind==='diagnostics')return 'R-hat assesses agreement between chains; bulk and tail ESS describe effective information for central estimates and posterior tails. BFMI assesses energy exploration for the fitted model. These checks concern sampling reliability, not biological significance.';
 if(d.section==='simulation')return 'Simulated datasets test estimation error, interval coverage and decision rules under known conditions. Compare scenarios within the same model, geometry and contrast. Coverage is the fraction of intervals containing the true effect; false discovery proportion and power answer different questions.';
 if(d.id.includes('celltype_purity'))return 'Clustering purity is the fraction of profiles for a cell type assigned to its dominant K-means cluster. Higher values indicate concentration in one cluster; they do not measure tumour purity or prove deconvolution accuracy.';
 if(d.id.includes('confusion_crosstab'))return 'Counts show how each cell type is distributed across K-means clusters. Cluster numbers have meaning only within a compartment; several cell types can share a cluster.';
 if(d.id.includes('separation_metrics'))return 'Silhouette measures separation of expression profiles; adjusted Rand index measures agreement between clusters and cell-type labels. Raw and row-normalized expression describe different geometries.';
 if(d.id.includes('patient_spread'))return 'Within-cell-type standard deviations along PC1, PC2 and PC3 describe variation between profiles. These values are spreads, not PCA coordinates.';
 if(d.id.includes('pca_component'))return 'Squared-axis contributions show which cell types contribute most to each principal component. Contributions are neither abundance fractions nor BMI effects.';
 if(d.id.includes('pca_explained'))return 'Explained variance quantifies how much expression variation each principal component captures. QC subsets are reported separately.';
 if(d.id.includes('pca_pc_bmi'))return 'Associations between principal components and BMI describe exploratory patterns in expression geometry. They are distinct from the primary Bayesian BMI effects.';
 if(d.id.includes('pca_tiered'))return 'This comparison assesses how PCA geometry changes under different component-quality thresholds. Non-immune results show greater sensitivity to the QC choice.';
 if(d.id.includes('pca_input'))return 'Numbers of profiles, features and cell types contributing to each PCA analysis.';
 if(d.id.includes('theta_bmi'))return 'Theta represents the estimated contribution of a cell type to bulk expression. Associations with BMI are reported with their corresponding multiple-testing adjustment families.';
 if(/inclusion|reliability_all/.test(d.id))return 'Component QC describes sparse expression, low estimated abundance and variability in deconvolution. Inclusion tiers apply different quality thresholds.';
 if(/dependence/.test(d.id))return 'Correlations describe dependence between programme scores. Normalization and residual adjustment change the contribution of shared expression patterns.';
 if(/raw_vs|cap3_vs/.test(d.id))return 'This sensitivity comparison shows how score correlations and BMI slopes change with expression scaling or clipping thresholds.';
 if(d.kind==='annotation')return 'Reference-based annotations assign single cells to cell-type labels. The reference panels are alternative annotations of the same data and should not be summed.';
 if(/concordance|agreement_with_chance/.test(d.id))return 'Agreement statistics compare effect patterns between data layers. Confidence intervals describe uncertainty; agreement does not establish cell-specific protein expression.';
 if(d.id.includes('celltype_mean_slopes'))return 'Mean associations are grouped by cell type and measurement layer. These are aggregate estimates with confidence intervals, distinct from per-feature posterior HDIs.';
 if(d.kind==='coverage')return 'Feature counts describe each stage from measurement to eligibility and adequate numerical sampling.';
 if(d.section==='ptm')return 'Post-translational modification summaries describe measured sites or glycoforms and their associations with BMI. Conditional effects account for protein abundance but are not direct occupancy measurements. Site identity, assay and uncertainty should be considered together.';
 return 'This table provides the corresponding analysis quantities and uncertainty measures. Select one comparison and measurement scale before interpreting differences.';
}

const collectionLabel=d=>{const part=d.id.split('/')[1];return ({effects:'Effect estimates',run_diagnostics:'Sampling overview',diagnostic_families:'R-hat & ESS by parameter family',diagnostic_family_summary:'R-hat & ESS by parameter family',ppc_domain_checks:'Predictive checks · domains',ppc_v2_domain_checks:'Predictive checks · domains',ppc_feature_checks:'Predictive checks · features',ppc_v2_feature_checks:'Predictive checks · features',fit_status:'Sampling checks · all fits',run_status:'Sampling checks · all fits',predictive_checks:'Posterior predictive checks',all_predictive_checks:'Posterior predictive checks'})[part]||d.title;};
function syncRunControls(d){if(d?.section==='rna'&&d.run_key!=='overview'){const categorical=d.run_key.startsWith('categorical_');$('bmi-model').value=categorical?'categorical':'continuous';$('compartment').value=d.run_key.replace('categorical_','');}}
function availableCollections(section){return state.catalog.filter(d=>d.section===section).filter(d=>{if(state.mode==='advanced')return d.title.toLowerCase().includes($('collection-search').value.trim().toLowerCase());if(state.mode==='biology'&&d.kind!=='effects')return false;if(state.mode==='reliability'&&d.kind!=='diagnostics')return false;if(section==='rna'&&d.id!=='rna/compartment_diagnostics'){const run=($('bmi-model').value==='categorical'?'categorical_':'')+$('compartment').value;return d.run_key===run;}return true;});}
function changeSection(section,preferred){
 fade($('results')); 
 if(['manuscript','simulation','reference','signatures','single_cell'].includes(section)&&state.mode!=='advanced'){state.mode='advanced';updateModeUI();}
 state.section=section;for(const b of $('sections').children)b.setAttribute('aria-selected',b.dataset.section===section);
 $('rna-navigation').hidden=state.mode==='advanced'||section!=='rna';
 const choices=availableCollections(section);const defaultId=section==='ptm'&&state.mode!=='reliability'&&choices.some(d=>d.id==='ptm/effects')?'ptm/effects':section==='reference'?'reference/pca_original':state.mode==='reliability'?(choices.find(d=>d.id==='rna/compartment_diagnostics')?.id||choices.find(d=>/run_diagnostics|fit_status|run_status/.test(d.id))?.id):choices.find(d=>d.kind==='effects')?.id;
 options($('dataset'),choices.map(d=>[d.id,state.mode==='advanced'?d.title:collectionLabel(d)]),choices.some(d=>d.id===preferred)?preferred:choices.some(d=>d.id===defaultId)?defaultId:choices[0]?.id);
 $('dataset').hidden=state.mode==='biology'&&choices.length===1;$('dataset-label').hidden=$('dataset').hidden;
 if(!choices.length){state.version++;state.dataset=null;$('search').disabled=true;$('reset').disabled=true;$('prev').disabled=true;$('next').disabled=true;state.rows=[];state.filtered=[];clearFigure();$('dataset-title').textContent='No matching collections';$('dataset-description').textContent='Clear the collection search to see the full archive.';$('download').disabled=true;$('download-json').hidden=true;$('diagnostic-panel').hidden=true;$('tbody').replaceChildren();$('facets').replaceChildren();$('extra-facets').replaceChildren();$('result-count').textContent='';$('draw').disabled=true;status('');return;}
 loadDataset();
}
function updateModeUI(){for(const b of document.querySelectorAll('[data-mode]'))b.setAttribute('aria-selected',b.dataset.mode===state.mode);$('collection-search-label').hidden=state.mode!=='advanced';$('extra-options').open=state.mode==='advanced';$('mode-help').textContent=state.mode==='biology'?'Start with a data layer and compartment. Primary models are selected by default; sensitivity analyses remain under More options.':state.mode==='reliability'?'Inspect sampling quality and predictive checks. Failed fits remain visible here. Summary diagnostics do not replace chain traces.':'All '+state.catalog.length+' collections: effect estimates, diagnostics, simulations, reference data and sensitivity analyses.';}
function setMode(mode){if(mode!=='advanced'&&['manuscript','reference','signatures','single_cell','simulation'].includes(state.section))state.section='rna';state.mode=mode;$('collection-search').value='';syncRunControls(state.dataset);updateModeUI();changeSection(state.section,state.dataset?.id);}
async function loadDataset(){
 const version=++state.version;state.dataset=state.catalog.find(d=>d.id===$('dataset').value);if(!state.dataset)return;
 const d=state.dataset;$('pca-panel').hidden=true;$('signature-detail').hidden=true;document.querySelector('.filter-card').hidden=false;document.querySelector('.viewbar').hidden=false;$('table-panel').hidden=false;$('analysis-explanation').textContent=explain(d);$('download-json').hidden=d.kind!=='signatures';$('download-json').disabled=true;$('search').disabled=false;$('reset').disabled=false;state.rows=[];state.filtered=[];state.facets={};state.page=0;clearFigure();$('search').value='';$('include-failed').checked=false;$('download').disabled=true;$('facets').replaceChildren();$('extra-facets').replaceChildren();$('diagnostic-panel').hidden=true;$('draw').disabled=true;$('tbody').replaceChildren();$('thead').replaceChildren();
 $('dataset-title').textContent=state.mode==='advanced'||d.run_key==='overview'?d.title:(d.section==='rna'?runLabel[d.run_key.replace('categorical_','')]:sections.find(s=>s[0]===d.section)[1])+' · '+collectionLabel(d);$('dataset-description').textContent=brief(d);$('section-label').textContent=sections.find(s=>s[0]===d.section)[1];$('notice').hidden=!d.provisional;$('notice').textContent=d.provisional?'Provisional analysis: interpret these fine-immune categorical results with the diagnostic evidence.':'';
 if(d.kind==='pca_original'){showPCA();return;}
 status('Loading '+number.format(d.row_count)+' result rows…');history.replaceState(null,'','#dataset='+encodeURIComponent(d.id)+'&mode='+state.mode);
 try{
  let rows=cache.get(d.id);
  if(!rows&&d.id==='rna/compartment_diagnostics'){rows=[];for(const run of ['immune_coarse','immune_fine','nonimmune','categorical_immune_coarse','categorical_immune_fine','categorical_nonimmune']){const ds=await api('atlas_web_rows',{select:'data',dataset_id:'eq.'+run+'/run_diagnostics',order:'row_number.asc',limit:10});if(version!==state.version)return;rows.push(...ds.map(x=>({...x.data,compartment:run.replace('categorical_',''),bmi_model:run.startsWith('categorical_')?'categorical':'continuous'})));}cache.set(d.id,rows);}
  if(!rows){rows=[];for(let offset=0;offset<d.row_count;offset+=1000){const page=await api('atlas_web_rows',{select:'row_number,data',dataset_id:'eq.'+d.id,order:'row_number.asc',limit:1000,offset});if(version!==state.version)return;if(!page.length)throw Error('Incomplete data response. Please retry.');rows.push(...page.map(x=>x.data));status('Loading results: '+number.format(rows.length)+' / '+number.format(d.row_count));}if(rows.length!==d.row_count)throw Error('The collection changed during loading. Please reload.');cache.set(d.id,rows);}
  if(version!==state.version)return;state.rows=rows;$('include-failed').checked=d.kind==='diagnostics';buildFacets();applyFilters();view(state.mode==='biology'&&plotTypes().length?'plot':'table');$('view-plot').hidden=!plotTypes().length;status('');
 }catch(e){if(version===state.version){status(e.message,true);$('result-count').textContent='';}}
}
const facetFields=['template','quantity','numerical_pass','bmi_model','entity_level','variant','component','contrast','units','cell_type','assay','geometry','kind','rule','scenario','check','domain','parameter_family','compartment','analysis_set','PC','CellType','reference','resolution','scope','endpoint','cap','mode','layer','metric','comparison','family','stratum'];
function buildFacets(preserve=false){
 const previous=preserve?{...state.facets}:{};state.facets={};
 $('facets').replaceChildren();$('extra-facets').replaceChildren();for(const key of facetFields){const eligible=state.rows.filter(r=>Object.entries(state.facets).every(([k,v])=>!v||String(r[k])===v));const values=[...new Set(eligible.map(r=>r[key]).filter(v=>v!==null&&v!==undefined&&v!==''))].sort();if(values.length<2||values.length>150)continue;
 const mandatory=!(state.dataset.id==='ptm/extended_sensitivity_effects'&&key==='variant')&&((state.dataset.id==='manuscript/ptm_calibration_templates'&&key==='template')||(state.dataset.id==='manuscript/sensitivity_comparisons'&&key==='comparison')||(state.dataset.kind==='signatures'&&key==='cell_type')||(state.dataset.id==='rna/compartment_diagnostics'&&key==='bmi_model')||(state.dataset.section==='simulation'&&['kind','geometry','contrast','rule'].includes(key))||(state.dataset.kind==='purity'&&key==='compartment')||(state.dataset.kind==='scree'&&['compartment','analysis_set'].includes(key))||(state.dataset.kind==='contributions'&&['compartment','analysis_set','PC'].includes(key))||(state.dataset.kind==='annotation'&&['reference','resolution'].includes(key))||(state.dataset.kind==='effects'&&['entity_level','variant','component','contrast','units'].includes(key))||(state.dataset.kind==='diagnostics'&&['variant','check','domain'].includes(key))||(/celltype_mean_slopes|agreement_with_chance|concordance_95CI/.test(state.dataset.id)&&['metric','comparison','layer','scope','stratum'].includes(key)));
 let initial=previous[key]&&values.map(String).includes(previous[key])?previous[key]:mandatory?(key==='component'&&values.includes('PTM_given_Protein')?'PTM_given_Protein':values.includes('continuous')?'continuous':values.includes('primary')?'primary':values.includes('feature')?'feature':values[0]):'';
 const fieldNames={variant:'Analysis version',component:'Quantity to compare',cell_type:'Cell type',entity_level:'Feature level',contrast:'BMI contrast',units:'Measurement scale',assay:'Assay'};const label=el('label',fieldNames[key]||human(key));const select=el('select',undefined,{'aria-label':fieldNames[key]||human(key)});options(select,[...(!mandatory?[['','All']]:[]),...values.map(v=>[String(v),human(v)])],String(initial));state.facets[key]=String(initial);select.dataset.key=key;select.onchange=()=>{state.facets[key]=select.value;state.page=0;buildFacets(true);applyFilters();};label.append(select);const extra=state.mode!=='advanced'&&['variant','entity_level','units','parameter_family'].includes(key);$(extra?'extra-facets':'facets').append(label);
 }
 $('include-failed').closest('label').hidden=state.dataset.kind==='diagnostics'||!state.rows.some(r=>r.numerical_checks_pass!==undefined);$('extra-options').hidden=!$('extra-facets').children.length;
}
function applyFilters(){
 const term=$('search').value.trim().toLowerCase();state.filtered=state.rows.filter(r=>Object.entries(state.facets).every(([k,v])=>!v||String(r[k])===v)&&(state.dataset.kind==='diagnostics'||$('include-failed').checked||![false,'False','false'].includes(r.numerical_checks_pass))&&(!term||searchable(r).some(v=>v!=null&&String(v).toLowerCase().includes(term))));
 state.page=0;clearFigure();$('download').disabled=!state.filtered.length;$('download-json').disabled=!state.filtered.length;$('result-count').textContent=number.format(state.filtered.length)+' / '+number.format(state.rows.length)+' rows';drawTable();setupPlot();renderDiagnostics();renderSignatureDetail();$('view-plot').hidden=!plotTypes().length;$('selection-summary').textContent=[state.dataset.kind==='effects'?'Pointwise 95% intervals':'Saved analysis summaries',...['variant','component','contrast','units'].map(k=>state.facets[k]?human(state.facets[k]):'').filter(Boolean)].join(' · ');
}
function drawTable(){
 const rows=state.filtered;const keys=[...new Set(rows.flatMap(Object.keys))];const first=state.dataset.kind==='signatures'?['signature','gene_count','genes','cell_type']:state.dataset.kind==='effects'?['label','cell_type','gene','estimate','hdi_lower','hdi_upper','units','contrast','component','variant','numerical_checks_pass']:[];const cols=state.mode==='biology'&&state.dataset.kind==='effects'?['label','cell_type','gene','estimate','hdi_lower','hdi_upper','units'].filter(k=>keys.includes(k)):[...first.filter(k=>keys.includes(k)),...keys.filter(k=>!first.includes(k))];
 const tr=el('tr');for(const c of cols)tr.append(el('th',human(c),{scope:'col'}));$('thead').replaceChildren(tr);$('tbody').replaceChildren();
 for(const row of rows.slice(state.page*25,(state.page+1)*25)){const tr=el('tr');for(const c of cols){const td=el('td',safeValue(c==='genes'&&Array.isArray(row[c])?row[c].join(', '):c==='label'&&state.mode!=='advanced'?programmeLabel(row):row[c]));if(c==='genes')td.className='genes';if(c==='numerical_checks_pass'&&[false,'False','false'].includes(row[c]))td.className='flag';tr.append(td);}$('tbody').append(tr);}
 if(!rows.length){const tr=el('tr');tr.append(el('td','No results match these filters. Try resetting the filters.'));$('tbody').append(tr);}
 const pages=Math.max(1,Math.ceil(rows.length/25));$('page-label').textContent=`Page ${state.page+1} of ${pages}`;$('prev').disabled=state.page===0;$('next').disabled=state.page>=pages-1;
}
function plotTypes(){
 const d=state.dataset;if(!d)return [];
 if(/ptm_heatmap_lookup|phospho_labels|glyco_labels/.test(d.id))return [['ptm_primary','PTM heatmap'],['ptm_colourbar','Shared colour scale']];
 if(/rna_heatmap_values|rna_protein_heatmap_values/.test(d.id))return [['manuscript_heatmap','Heatmap']];
 if(/ptm_calibration_summary|ptm_calibration_templates/.test(d.id))return [['calibration_coverage','Interval coverage'],['calibration_bias','Bias'],['calibration_rmse','RMSE']];
 if(d.id==='manuscript/sensitivity_comparisons')return [['sensitivity_comparison','Matched effect intervals']];
 if(d.id==='manuscript/sensitivity_classifications')return [['sensitivity_classification','Interval classifications']];
 if(d.id==='manuscript/ptm_sbc_bins')return [['sbc_ranks','SBC rank histograms']];
 if(d.id==='manuscript/ptm_calibration_completion')return [['calibration_completion','Calibration sample counts']];
 if(d.id==='ptm/extended_sensitivity_effects')return [['sensitivity_grid','Sensitivity heatmap'],['forest','Effect intervals'],['heatmap','Heatmap']];
 if(d.kind==='effects')return [['forest','Effect intervals'],['heatmap','Heatmap']];
 if(d.kind==='signatures')return [['signatures','Genes per signature']];
 if(d.kind==='coverage')return [['coverage','Feature coverage']];
 if(d.kind==='purity')return [['purity','Cell-type clustering purity']];
 if(d.kind==='annotation')return [['annotation','Reference annotation counts']];
 if(d.section==='simulation')return [['simulation','Simulation metric by scenario']];
 if(/celltype_mean_slopes|agreement_with_chance|concordance_95CI/.test(d.id))return [['intervals','Estimates and confidence intervals']];
 return [];
}
function setupPlot(){
 const types=plotTypes(),old=$('plot-type').value;options($('plot-type'),types,types.some(t=>t[0]===old)?old:types[0]?.[0]);
 const metrics=['coverage','coverage_mean','rmse','bias','mean_absolute_error','mean_interval_width','mean_hdi_width','fdp','fdp_mean','power','power_mean','seconds','discoveries','discoveries_mean','power_mean_BHM','power_mean_OLS','fdp_mean_BHM','fdp_mean_OLS'];
 options($('plot-y'),metrics.filter(k=>state.filtered.some(r=>r[k]!==''&&r[k]!=null&&Number.isFinite(Number(r[k])))).map(k=>[k,human(k)]));plotControls();
}
function plotControls(){const t=$('plot-type').value;$('x-wrap').hidden=true;$('y-wrap').hidden=t!=='simulation';$('plot-y').parentElement.firstChild.textContent='Metric';$('plot-limit').parentElement.hidden=!['forest','heatmap'].includes(t);$('draw').disabled=!state.filtered.length||!plotTypes().length;}
function fade(node){node.classList.remove('fade-in');void node.offsetWidth;node.classList.add('fade-in');}
function view(which){if(which==='plot'&&!plotTypes().length)which='table';state.view=which;$('table-panel').hidden=which!=='table';$('plot-panel').hidden=which!=='plot';$('view-table').setAttribute('aria-selected',which==='table');$('view-plot').setAttribute('aria-selected',which==='plot');fade($(which==='plot'?'plot-panel':'table-panel'));}
let busy=false,activePlotId=null,pendingPlot=null;
function dispatchPlot(job){busy=true;activePlotId=job.id;state.worker.postMessage(job);}
function startWorker(){
 if(state.worker)return;
 state.worker=new Worker('./plot-worker.js?v=manuscript-python-2',{type:'module'});
 state.worker.onmessage=({data})=>{
  if(data.status){if(data.id===state.plotId)status(data.status);return;}
  if(data.id===activePlotId){busy=false;activePlotId=null;}
  if(data.id===state.plotId){
   $('draw').disabled=!state.filtered.length||!plotTypes().length;$('cancel-plot').hidden=true;
   if(data.error)status(data.error.split('\n').filter(Boolean).slice(-1)[0],true);
   else{state.figure=data.result;const img=el('img',undefined,{alt:state.dataset.title+' — '+$('plot-type').selectedOptions[0].text});const url=URL.createObjectURL(new Blob([data.result.svg],{type:'image/svg+xml'}));img.onload=()=>URL.revokeObjectURL(url);img.src=url;$('figure-container').replaceChildren(img);$('plot-caption').textContent=data.result.caption;$('svg').disabled=$('png').disabled=false;status('Figure ready.');}
  }
  if(pendingPlot){const next=pendingPlot;pendingPlot=null;if(next.id===state.plotId)dispatchPlot(next);}
 };
 state.worker.onerror=()=>{busy=false;activePlotId=null;pendingPlot=null;$('draw').disabled=!state.filtered.length;$('cancel-plot').hidden=true;status('Plotting stopped. Choose Draw figure to restart it; the data remain available.',true);state.worker.terminate();state.worker=null;};
}
function drawPlot(){
 if(!plotTypes().length||!state.filtered.length)return;
 const id=++state.plotId;startWorker();$('draw').disabled=true;$('cancel-plot').hidden=false;
 const job={id,request:{scale_rows:state.dataset.id==='manuscript/ptm_heatmap_lookup'?state.rows:undefined,rows:state.filtered.map(r=>r.feature&&PROGRAMME_LABELS[r.feature]?{...r,label:programmeLabel(r)}:r),type:$('plot-type').value,y:$('plot-y').value,limit:Number($('plot-limit').value),title:state.dataset.title}};
 if(busy){pendingPlot=job;status('Figure queued. You can continue browsing or return Home.');}else{status('Preparing the figure. Navigation remains available.');dispatchPlot(job);}
}
$('cancel-plot').onclick=()=>{clearFigure();$('draw').disabled=!state.filtered.length;status('Figure discarded. The plotting engine remains loaded.');};
$('manuscript-shortcut').onclick=()=>{state.mode='advanced';$('collection-search').value='';updateModeUI();changeSection('manuscript','manuscript/ptm_heatmap_lookup');};
function showPCA(){
 state.rows=[];state.filtered=[];clearFigure();$('pca-panel').hidden=false;$('diagnostic-panel').hidden=true;document.querySelector('.filter-card').hidden=true;document.querySelector('.viewbar').hidden=true;$('table-panel').hidden=true;$('plot-panel').hidden=true;$('download').disabled=true;status('');history.replaceState(null,'','#dataset=reference%2Fpca_original&mode=advanced');loadPCAImage();
}
function loadPCAImage(){const c=$('pca-compartment').value,m=$('pca-scaling').value;const url='./pca/'+c+'_'+m+'.png';$('pca-image').src=url;$('pca-image').alt=runLabel[c]+' PCA, '+(m==='raw'?'raw expression':'row-normalized expression');fade($('pca-panel'));}
$('pca-compartment').onchange=$('pca-scaling').onchange=loadPCAImage;
function renderSignatureDetail(){const box=$('signature-detail');box.hidden=state.dataset?.kind!=='signatures';if(box.hidden)return;box.replaceChildren(el('h3','Signature composition'),el('p','Select a cell type, then search by signature name or gene. The table lists the genes in each signature; the chart compares gene counts. Counts describe the gene sets, not expression or enrichment.'));}
function home(){state.mode='biology';$('collection-search').value='';$('search').value='';$('compartment').value='immune_coarse';$('bmi-model').value='continuous';clearFigure();updateModeUI();changeSection('rna','immune_coarse/effects');window.scrollTo({top:0,behavior:'auto'});}
$('home').onclick=home;

const diagnosticSpecs=[
 {id:'rhat',label:'R-hat · maximum',keys:['max_all_parameter_rhat','maximum_rhat','max_rhat'],higherBad:true,reference:1.01,note:'Each point is the reported maximum R-hat for one parameter family or fit. Values near 1 indicate agreement between chains. The dashed line is 1.01; this is a diagnostic reference, not an overall model-validity verdict.'},
 {id:'bulk',label:'Bulk ESS · minimum',keys:['min_all_parameter_ess_bulk','minimum_bulk_ess','min_ess_bulk'],higherBad:false,note:'Each point is the reported minimum bulk effective sample size for one family or fit, not a chain or patient count. It describes information available for central posterior estimates; larger values generally mean better Monte Carlo precision.'},
 {id:'tail',label:'Tail ESS · minimum',keys:['min_all_parameter_ess_tail','minimum_tail_ess','min_ess_tail'],higherBad:false,note:'Each point is the reported minimum tail effective sample size for one family or fit. Tail ESS concerns precision in the posterior tails and should be reviewed separately from bulk ESS.'},
 {id:'bfmi',label:'BFMI · minimum across chains',keys:['min_bfmi'],higherBad:false,note:'Minimum reported BFMI across chains. This summarizes energy exploration; it does not show individual chain trajectories.'},
 {id:'divergences',label:'Divergences · count',keys:['divergences'],higherBad:true,note:'Recorded divergence counts for each saved fit. Zero is desirable. No percentage is inferred when the required draw totals are unavailable.'},
 {id:'ppc',label:'Observed versus predictive interval',keys:['observed_draw_median','observed_count','observed_median'],higherBad:true,note:'Observed summaries (magenta points) compared with the saved predictive interval (blue lines). Different predictive statistics must be viewed separately. These checks assess model fit, not convergence or calibrated discovery probabilities.'},
 {id:'tail-rate',label:'Predictive tail outside rate',keys:['tail_outside_rate'],higherBad:true,note:'Saved predictive tail outside rate per feature. This is a model-check summary, not an MCMC diagnostic or a discovery probability.'},
 {id:'coverage',label:'Predictive 90% interval coverage',keys:['predictive_90pct_coverage'],higherBad:false,reference:.9,note:'Saved predictive coverage per feature, with a 0.90 reference. Coverage and MCMC convergence answer different questions.'}
];
const diagnosticNumber=(r,keys)=>{for(const k of keys){const v=r[k];if(v!==null&&v!==undefined&&v!==''&&Number.isFinite(Number(v)))return Number(v);}return null;};
const isFalse=v=>v===false||v==='False'||v==='false';
const isTrue=v=>v===true||v==='True'||v==='true';
let activeDiagnosticSpecs=[],diagnosticSVG='';
function renderDiagnostics(){
 const visible=state.dataset?.kind==='diagnostics';$('diagnostic-panel').hidden=!visible;if(!visible)return;
 const rows=state.filtered;$('diagnostic-cards').replaceChildren();
 const card=(label,value,review=false)=>{const box=el('div',undefined,{class:'metric-card'+(review?' review':'')});box.append(el('strong',value),el('span',label));$('diagnostic-cards').append(box);};
 card('Saved summary rows',number.format(rows.length));
 if(!rows.length){$('diagnostic-explanation').textContent='No diagnostic rows match the filters.';$('diagnostic-chart').replaceChildren();$('diagnostic-metric-label').hidden=true;$('diagnostic-svg').hidden=true;return;}
 const flags=rows.filter(r=>r.numerical_checks_pass!==undefined||r.gate_pass!==undefined||r.publication_convergence_pass!==undefined);
 if(flags.length){const failed=flags.filter(r=>isFalse(r.numerical_checks_pass)||isFalse(r.gate_pass)||isFalse(r.publication_convergence_pass)).length;card('Fits failing a recorded numerical gate',failed+' / '+flags.length,failed>0);}
 if(rows.length===1){const r=rows[0];const chains=diagnosticNumber(r,['chains','observed_chains']);const draws=diagnosticNumber(r,['draws_per_chain','observed_draws_per_chain']);if(chains!==null)card('Recorded chains',String(chains));if(draws!==null)card('Draws per chain',number.format(draws));}
 for(const spec of diagnosticSpecs.slice(0,5)){const values=rows.map(r=>diagnosticNumber(r,spec.keys)).filter(v=>v!==null);if(!values.length)continue;const value=spec.higherBad?Math.max(...values):Math.min(...values);card(spec.id==='divergences'&&rows.length>1?'Maximum divergences in one fit':spec.label,number.format(value),spec.id==='rhat'?value>1.01:spec.id==='divergences'?value>0:false);}
 const predictiveFlags=rows.filter(r=>isTrue(r.flag)||isTrue(r.residual_dispersion_outside)||isTrue(r.studentized_tail_outside)||isTrue(r.flag_upper)||isTrue(r.review_flag_p_le_0_025)).length;
 if(rows.some(r=>['flag','residual_dispersion_outside','flag_upper','review_flag_p_le_0_025'].some(k=>k in r)))card('Rows with a recorded predictive flag',String(predictiveFlags),predictiveFlags>0);
 activeDiagnosticSpecs=diagnosticSpecs.filter(s=>rows.some(r=>diagnosticNumber(r,s.keys)!==null));const old=$('diagnostic-metric').value;
 const combined=activeDiagnosticSpecs.some(s=>['rhat','bulk','tail'].includes(s.id));options($('diagnostic-metric'),[...(combined?[['overview','R-hat, bulk ESS and tail ESS']]:[]),...activeDiagnosticSpecs.map(s=>[s.id,s.label])],old==='overview'&&combined?'overview':activeDiagnosticSpecs.some(s=>s.id===old)?old:combined?'overview':activeDiagnosticSpecs[0]?.id);
 $('diagnostic-metric-label').hidden=!activeDiagnosticSpecs.length;$('chain-note').hidden=!activeDiagnosticSpecs.some(s=>['rhat','bulk','tail','bfmi','divergences'].includes(s.id));
 renderDiagnosticChart();
}
function renderDiagnosticChart(){
 if($('diagnostic-metric').value==='overview'){renderMCMCOverview();return;}
 const spec=activeDiagnosticSpecs.find(s=>s.id===$('diagnostic-metric').value);const root=$('diagnostic-chart');root.replaceChildren();diagnosticSVG='';$('diagnostic-svg').hidden=true;
 if(!spec){$('diagnostic-explanation').textContent='Review the recorded flags and full table below. No compatible quantitative diagnostic is available for this selection.';return;}
 let points=state.filtered.map((r,i)=>({r,i,value:diagnosticNumber(r,spec.keys)})).filter(p=>p.value!==null);
 if(spec.id==='ppc'){
  const checks=new Set(points.map(p=>p.r.check||p.r.domain||''));if(checks.size>1){$('diagnostic-explanation').textContent='Select one predictive check or domain above to compare observed values with the corresponding predictive interval.';return;}
  for(const p of points){p.low=diagnosticNumber(p.r,['replicated_95_lower','replicated_q025']);p.high=diagnosticNumber(p.r,['replicated_95_upper','replicated_q975']);p.outside=p.low!==null&&p.high!==null&&(p.value<p.low||p.value>p.high);}
  points.sort((a,b)=>Number(b.outside)-Number(a.outside));
 }else points.sort((a,b)=>spec.higherBad?b.value-a.value:a.value-b.value);
 const total=points.length;if(!total)return;if(spec.id==='ppc'&&total>80){renderPredictiveScatter(points,spec);return;}if(total>80&&spec.id!=='ppc'){renderDiagnosticDistribution(spec,points);return;}if(spec.id==='bfmi'&&total===1){$('diagnostic-explanation').textContent=spec.note+' This scalar is shown in the summary above. Energy-density plots require draw-level energy values.';return;}
 $('diagnostic-explanation').textContent=spec.note+' '+'All '+total+' summary rows shown.';
 const width=860,left=330,right=40,top=42,step=29,height=top+points.length*step+65;
 const vals=points.flatMap(p=>[p.value,p.low,p.high]).filter(v=>v!==null&&v!==undefined);if(spec.reference!==undefined)vals.push(spec.reference);
 let min=spec.id==='rhat'?Math.min(.99,...vals):Math.min(0,...vals),max=Math.max(...vals);if(max===min)max=min+1;const pad=(max-min)*.04;max+=pad;if(min<0)min-=pad;
 const x=v=>left+(v-min)/(max-min)*(width-left-right);
 const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.setAttribute('role','img');svg.setAttribute('aria-label',spec.label+' by saved family or fit');svg.setAttribute('xmlns',ns);
 const add=(name,attrs,text)=>{const n=document.createElementNS(ns,name);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;svg.append(n);return n;};
 add('rect',{width,height,fill:'white'});add('text',{x:14,y:20,'font-size':13,'font-family':'Arial',fill:'#1F1F1F'},spec.label);
 for(let i=0;i<=4;i++){const value=min+(max-min)*i/4,xx=x(value);add('line',{x1:xx,x2:xx,y1:top-10,y2:height-45,stroke:'#e5e8e8'});add('text',{x:xx,y:height-23,'text-anchor':'middle','font-family':'Arial','font-size':11,fill:'#6F7378'},number.format(value));}
 if(spec.reference!==undefined)add('line',{x1:x(spec.reference),x2:x(spec.reference),y1:top-12,y2:height-45,stroke:'#6F7378','stroke-dasharray':'5 4'});
 points.forEach((p,i)=>{const yy=top+i*step,label=human(p.r.compartment||p.r.parameter_family||programmeLabel(p.r)||p.r.key||p.r.global_feature||p.r.global_celltype||p.r.celltype||p.r.global_celltype||'Run summary');const t=add('text',{x:left-12,y:yy+4,'text-anchor':'end','font-size':11,'font-family':'Arial',fill:'#1F1F1F'},label.length>46?label.slice(0,43)+'…':label);const title=document.createElementNS(ns,'title');title.textContent=label+' = '+p.value;t.append(title);if(spec.id==='ppc'&&p.low!==null&&p.high!==null)add('line',{x1:x(p.low),x2:x(p.high),y1:yy,y2:yy,stroke:'#2C5AA0','stroke-width':3});const circle=add('circle',{cx:x(p.value),cy:yy,r:4,fill:spec.id==='ppc'||(spec.id==='rhat'&&p.value>1.01)||(spec.id==='divergences'&&p.value>0)?'#D81B60':'#0E7C7B'});const ct=document.createElementNS(ns,'title');ct.textContent=String(p.value);circle.append(ct);});
 root.append(svg);diagnosticSVG=new XMLSerializer().serializeToString(svg);$('diagnostic-svg').hidden=false;
}

function svgRoot(width,height,title){
 const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.setAttribute('xmlns',ns);svg.setAttribute('role','img');svg.setAttribute('aria-label',title);
 const add=(name,attrs,text)=>{const n=document.createElementNS(ns,name);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);if(text!==undefined)n.textContent=text;svg.append(n);return n;};
 add('rect',{width,height,fill:'white'});return {svg,add};
}
function showDiagnosticSVG(svg,note){$('diagnostic-chart').replaceChildren(svg);diagnosticSVG=new XMLSerializer().serializeToString(svg);$('diagnostic-svg').hidden=false;$('diagnostic-explanation').textContent=note;}
function renderMCMCOverview(){
 const rows=state.filtered, specs=diagnosticSpecs.slice(0,3),large=rows.length>80;
 const W=1000,H=large?350:Math.max(280,rows.length*29+130),{svg,add}=svgRoot(W,H,'R-hat and effective sample sizes across the selection');
 const text=(x,y,t,anchor='start')=>add('text',{x,y,'font-family':'Arial','font-size':12,'text-anchor':anchor,fill:'#1F1F1F'},t);
 const L=large?50:250,panel=(W-L-35)/3;
 for(let j=0;j<3;j++){
  const spec=specs[j],x0=L+j*panel+15,x1=L+(j+1)*panel-25;let points=rows.map((r,i)=>({i,r,v:diagnosticNumber(r,spec.keys)})).filter(p=>p.v!==null);if(!points.length)continue;
  const lo=spec.id==='rhat'?Math.min(.99,...points.map(p=>p.v)):0,hi=Math.max(...points.map(p=>p.v),spec.id==='rhat'?1.01:1),max=hi+(hi-lo)*.07;const x=v=>x0+(v-lo)/(max-lo)*(x1-x0);
  text(x0,24,['Maximum R-hat','Minimum bulk ESS','Minimum tail ESS'][j]);
  if(large){const bins=20,counts=Array(bins).fill(0);for(const p of points)counts[Math.min(bins-1,Math.floor((p.v-lo)/(max-lo)*bins))]++;const m=Math.max(...counts,1);counts.forEach((n,i)=>add('rect',{x:x0+i*(x1-x0)/bins,y:270-n/m*205,width:(x1-x0)/bins-1,height:n/m*205,fill:j===0?'#0E7C7B':j===1?'#2C5AA0':'#D81B60'}));text(x0,48,points.length+' fits');text(x0,325,'Distribution of fit summaries');}
  else{
   points.forEach(p=>{const yy=70+p.i*29;add('circle',{cx:x(p.v),cy:yy,r:4,fill:j===0&&p.v>1.01?'#D81B60':j===2?'#D81B60':'#2C5AA0'});});
   if(j===0)rows.forEach((r,i)=>text(L-8,74+i*29,human(r.compartment||r.parameter_family||r.programme||r.key||'Entire model').slice(0,37),'end'));
  }
  const bottom=large?270:H-55;add('line',{x1:x0,x2:x1,y1:bottom,y2:bottom,stroke:'#555'});
  for(let k=0;k<=3;k++){const v=lo+(max-lo)*k/3;text(x(v),bottom+22,number.format(v),'middle');}
  if(spec.id==='rhat')add('line',{x1:x(1.01),x2:x(1.01),y1:45,y2:bottom,stroke:'#777','stroke-dasharray':'4 4'});
 }
 showDiagnosticSVG(svg,'All '+rows.length+' selected '+(large?'fits are included in these distributions.':'parameter-family or compartment summaries are shown together.')+' R-hat compares chains (reference 1.01); bulk and tail ESS combine information across chains. The values are recorded maxima/minima, not individual parameter distributions or single-chain ESS. BFMI is a fit-level energy diagnostic and is reported separately.');
}
function renderPredictiveScatter(points,spec){
 const vals=points.map(p=>({...p,median:diagnosticNumber(p.r,['replicated_median'])})).filter(p=>p.median!==null);
 if(!vals.length){$('diagnostic-explanation').textContent='The complete predictive-check table is shown below. An observed-versus-predicted plot requires saved predictive medians.';return;}
 const W=780,H=490,{svg,add}=svgRoot(W,H,'Observed versus predictive median');const all=vals.flatMap(p=>[p.value,p.median]);let lo=Math.min(...all),hi=Math.max(...all);if(lo===hi)hi=lo+1;const pad=(hi-lo)*.04;lo-=pad;hi+=pad;const x=v=>70+(v-lo)/(hi-lo)*650,y=v=>420-(v-lo)/(hi-lo)*360;
 add('line',{x1:x(lo),x2:x(hi),y1:y(lo),y2:y(hi),stroke:'#777','stroke-dasharray':'4 4'});
 vals.forEach(p=>add('circle',{cx:x(p.median),cy:y(p.value),r:2.3,fill:p.outside?'#D81B60':'#2C5AA0','fill-opacity':.4}));
 for(let i=0;i<=4;i++){const v=lo+(hi-lo)*i/4;add('text',{x:x(v),y:445,'text-anchor':'middle','font-family':'Arial','font-size':12},number.format(v));add('text',{x:60,y:y(v),'text-anchor':'end','font-family':'Arial','font-size':12},number.format(v));}
 add('text',{x:370,y:478,'text-anchor':'middle','font-family':'Arial','font-size':13},'Predictive median');add('text',{x:70,y:22,'font-family':'Arial','font-size':13},'Observed summary');
 showDiagnosticSVG(svg,'All '+vals.length+' comparable checks are shown. The diagonal indicates equality of observed and predictive medians. Magenta identifies observed values outside the saved predictive interval; blue indicates values inside it. This is a posterior predictive comparison, not a convergence diagnostic.');
}


function renderDiagnosticDistribution(spec,points){
 const W=850,H=350,{svg,add}=svgRoot(W,H,spec.label+' across all fits');const lo=Math.min(0,...points.map(p=>p.value)),hi=Math.max(...points.map(p=>p.value),1),bins=30,counts=Array(bins).fill(0);points.forEach(p=>counts[Math.min(bins-1,Math.floor((p.value-lo)/(hi-lo)*bins))]++);const max=Math.max(...counts,1),x=v=>65+(v-lo)/(hi-lo)*735;
 add('text',{x:65,y:24,'font-family':'Arial','font-size':14},spec.label+' across '+points.length+' fits');counts.forEach((n,i)=>add('rect',{x:65+i*735/bins,y:285-230*n/max,width:735/bins-1,height:230*n/max,fill:'#0E7C7B'}));
 for(let i=0;i<=5;i++){const v=lo+(hi-lo)*i/5;add('text',{x:x(v),y:310,'text-anchor':'middle','font-size':12,'font-family':'Arial'},number.format(v));}
 add('text',{x:15,y:45,'font-size':12,'font-family':'Arial'},'Fits');add('text',{x:20,y:65,'font-size':12,'font-family':'Arial'},String(max));
 showDiagnosticSVG(svg,spec.note+' Distribution across all '+points.length+' selected fits. These values are fit summaries, not posterior draws.');
}

$('diagnostic-svg').onclick=()=>save(diagnosticSVG,'atlas_'+$('diagnostic-metric').value+'.svg','image/svg+xml');

$('dataset').onchange=loadDataset;document.querySelectorAll('[data-mode]').forEach(b=>b.onclick=()=>setMode(b.dataset.mode));$('compartment').onchange=$('bmi-model').onchange=()=>changeSection(state.section);$('collection-search').oninput=()=>changeSection(state.section);$('diagnostic-metric').onchange=renderDiagnosticChart;$('search').oninput=applyFilters;$('include-failed').onchange=applyFilters;$('reset').onclick=()=>{$('search').value='';$('include-failed').checked=state.dataset?.kind==='diagnostics';buildFacets();applyFilters();};$('prev').onclick=()=>{state.page--;drawTable();};$('next').onclick=()=>{state.page++;drawTable();};$('view-table').onclick=()=>view('table');$('view-plot').onclick=()=>view('plot');$('plot-type').onchange=()=>{clearFigure();plotControls();};$('draw').onclick=drawPlot;
for(const id of ['download-json','download','svg','png','pca-download','manifest-download','diagnostic-svg']) $(id).onclick=requestDownload;
$('download-contact-close').onclick=()=> $('download-contact').close();
for(const [i,[id,name]]of sections.entries()){const b=el('button',undefined,{role:'tab','aria-selected':id===state.section});b.dataset.section=id;b.append(el('span',String(i+1).padStart(2,'0'),{class:'num'}),el('span',name));b.onclick=()=>changeSection(id);$('sections').append(b);}
(async()=>{try{state.catalog=await api('atlas_web_datasets',{select:'id,section,title,kind,description,run_key,row_count,provisional',order:'id.asc',limit:200});for(const d of state.catalog){if(d.kind==='signatures')d.title='Gene signature catalogue';}if(state.catalog.length)state.catalog.push({id:'reference/pca_original',section:'reference',title:'PCA by cell type',kind:'pca_original',row_count:0,run_key:'reference'},{id:'rna/compartment_diagnostics',section:'rna',title:'All RNA compartments: sampling overview',kind:'diagnostics',row_count:6,run_key:'overview'});if(!state.catalog.length){$('connection').textContent='Website ready · data release pending';$('dataset-title').textContent='Collections awaiting publication';status('The results have been prepared. Public access will open after the lab approves the data release.');$('search').disabled=true;$('reset').disabled=true;return;}$('connection').textContent=state.catalog.length+' collections · connected to the research data';const p=new URLSearchParams(location.hash.slice(1));const d=state.catalog.find(d=>d.id===p.get('dataset'));state.mode=['biology','reliability','advanced'].includes(p.get('mode'))?p.get('mode'):(d&&d.kind!=='effects'?'advanced':'biology');syncRunControls(d);updateModeUI();changeSection(d?.section||'rna',d?.id||'immune_coarse/effects');}catch(e){$('connection').textContent='Connection unavailable';$('dataset-title').textContent='Results temporarily unavailable';status(e.message,true);}})();
