import {CONFIG} from './config.js';
const $=id=>document.getElementById(id);
const sections=[['rna','RNA programmes'],['rna_protein','RNA–protein'],['ptm','Protein & PTM'],['simulation','Simulations']];
const state={catalog:[],section:'rna',dataset:null,rows:[],filtered:[],facets:{},page:0,view:'table',version:0,figure:null,worker:null,plotId:0};
const cache=new Map();
const number=new Intl.NumberFormat('en',{maximumSignificantDigits:5});
const human=s=>String(s).replaceAll('_',' ');
function el(tag,text,attrs={}){const e=document.createElement(tag);if(text!==undefined)e.textContent=text;for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);return e;}
function status(message,error=false){$('status').textContent=message;$('status').classList.toggle('error',error);}
async function api(table,query){const url=new URL(CONFIG.url+'/rest/v1/'+table);for(const[k,v]of Object.entries(query))url.searchParams.set(k,v);const res=await fetch(url,{headers:{apikey:CONFIG.key},signal:AbortSignal.timeout(30000)});if(!res.ok)throw Error('Results could not be loaded ('+res.status+'). Please try again.');return res.json();}
function options(select,values,selected){select.replaceChildren();for(const [value,label]of values){const o=el('option',label,{value});select.append(o);}if(selected!==undefined)select.value=selected;}
function safeValue(v){if(v===null||v===undefined||v==='')return '—';if(typeof v==='number')return number.format(v);if(typeof v==='object')return JSON.stringify(v);return String(v);}
function save(content,name,type){const url=URL.createObjectURL(content instanceof Blob?content:new Blob([content],{type}));const a=el('a',undefined,{href:url,download:name});document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
function csvCell(v){let s=v==null?'':typeof v==='object'?JSON.stringify(v):String(v);if(/^[=+@\t\r]/.test(s)||(/^-.+/.test(s)&&!Number.isFinite(Number(s))))s="'"+s;return '"'+s.replaceAll('"','""')+'"';}
function filteredCSV(rows){const cols=[...new Set(rows.flatMap(Object.keys))];return [cols.map(csvCell).join(','),...rows.map(r=>cols.map(c=>csvCell(r[c])).join(','))].join('\r\n');}
function clearFigure(){state.plotId++;state.figure=null;$('svg').disabled=$('png').disabled=true;$('plot-caption').textContent='';$('figure-container').replaceChildren(el('div','Choose a figure type and draw the current selection.',{class:'plot-empty'}));}
function changeSection(section,preferred){state.section=section;for(const b of $('sections').children)b.setAttribute('aria-selected',b.dataset.section===section);const choices=state.catalog.filter(d=>d.section===section);options($('dataset'),choices.map(d=>[d.id,d.title]),preferred&&choices.some(d=>d.id===preferred)?preferred:choices[0]?.id);loadDataset();}
async function loadDataset(){
 const version=++state.version;state.dataset=state.catalog.find(d=>d.id===$('dataset').value);if(!state.dataset)return;
 const d=state.dataset;state.rows=[];state.filtered=[];state.facets={};state.page=0;clearFigure();$('search').value='';$('include-failed').checked=false;$('download').disabled=true;$('facets').replaceChildren();$('tbody').replaceChildren();$('thead').replaceChildren();
 $('dataset-title').textContent=d.title;$('dataset-description').textContent=d.description;$('section-label').textContent=sections.find(s=>s[0]===d.section)[1].toUpperCase();$('notice').hidden=!d.provisional;$('notice').textContent=d.provisional?'Provisional analysis: interpret these fine-immune categorical results with the diagnostic evidence.':'';
 status('Loading '+number.format(d.row_count)+' result rows…');history.replaceState(null,'','#dataset='+encodeURIComponent(d.id));
 try{
  let rows=cache.get(d.id);
  if(!rows){rows=[];for(let offset=0;offset<d.row_count;offset+=1000){const page=await api('atlas_web_rows',{select:'row_number,data',dataset_id:'eq.'+d.id,order:'row_number.asc',limit:1000,offset});if(version!==state.version)return;if(!page.length)throw Error('Incomplete data response. Please retry.');rows.push(...page.map(x=>x.data));status('Loading results: '+number.format(rows.length)+' / '+number.format(d.row_count));}if(rows.length!==d.row_count)throw Error('The collection changed during loading. Please reload.');cache.set(d.id,rows);}
  if(version!==state.version)return;state.rows=rows;buildFacets();applyFilters();status('');
 }catch(e){if(version===state.version){status(e.message,true);$('result-count').textContent='';}}
}
const facetFields=['entity_level','variant','component','contrast','units','cell_type','assay','geometry','kind','rule','scenario','check','parameter_family'];
function buildFacets(preserve=false){
 const previous=preserve?{...state.facets}:{};state.facets={};
 $('facets').replaceChildren();for(const key of facetFields){const eligible=state.rows.filter(r=>Object.entries(state.facets).every(([k,v])=>!v||String(r[k])===v));const values=[...new Set(eligible.map(r=>r[key]).filter(v=>v!==null&&v!==undefined&&v!==''))].sort();if(values.length<2||values.length>150)continue;
 const mandatory=state.dataset.kind==='effects'&&['entity_level','variant','component','contrast','units'].includes(key);
 let initial=previous[key]&&values.map(String).includes(previous[key])?previous[key]:mandatory?(values.includes('primary')?'primary':values.includes('feature')?'feature':values[0]):'';
 const label=el('label',human(key));const select=el('select',undefined,{'aria-label':human(key)});options(select,[...(!mandatory?[['','All']]:[]),...values.map(v=>[String(v),human(v)])],String(initial));state.facets[key]=String(initial);select.dataset.key=key;select.onchange=()=>{state.facets[key]=select.value;state.page=0;buildFacets(true);applyFilters();};label.append(select);$('facets').append(label);
 }
 $('include-failed').closest('label').hidden=!state.rows.some(r=>r.numerical_checks_pass!==undefined);
}
function applyFilters(){
 const term=$('search').value.trim().toLowerCase();state.filtered=state.rows.filter(r=>Object.entries(state.facets).every(([k,v])=>!v||String(r[k])===v)&&($('include-failed').checked||![false,'False','false'].includes(r.numerical_checks_pass))&&(!term||Object.values(r).some(v=>v!=null&&String(v).toLowerCase().includes(term))));
 state.page=0;clearFigure();$('download').disabled=!state.filtered.length;$('result-count').textContent=number.format(state.filtered.length)+' / '+number.format(state.rows.length)+' rows';drawTable();setupPlot();
}
function drawTable(){
 const rows=state.filtered;const keys=[...new Set(rows.flatMap(Object.keys))];const first=state.dataset.kind==='effects'?['label','cell_type','gene','estimate','hdi_lower','hdi_upper','units','contrast','component','variant','numerical_checks_pass']:[];const cols=[...first.filter(k=>keys.includes(k)),...keys.filter(k=>!first.includes(k))];
 const tr=el('tr');for(const c of cols)tr.append(el('th',human(c),{scope:'col'}));$('thead').replaceChildren(tr);$('tbody').replaceChildren();
 for(const row of rows.slice(state.page*25,(state.page+1)*25)){const tr=el('tr');for(const c of cols){const td=el('td',safeValue(row[c]));if(c==='numerical_checks_pass'&&[false,'False','false'].includes(row[c]))td.className='flag';tr.append(td);}$('tbody').append(tr);}
 if(!rows.length){const tr=el('tr');tr.append(el('td','No results match these filters. Try resetting the filters.'));$('tbody').append(tr);}
 const pages=Math.max(1,Math.ceil(rows.length/25));$('page-label').textContent=`Page ${state.page+1} of ${pages}`;$('prev').disabled=state.page===0;$('next').disabled=state.page>=pages-1;
}
function setupPlot(){
 const kind=state.dataset.kind;const old=$('plot-type').value;const types=kind==='effects'?[['forest','Effect intervals'],['heatmap','Original-style heatmap']]:kind==='coverage'?[['coverage','Coverage counts']]:[['histogram','Distribution'],['scatter','Scatterplot']];options($('plot-type'),types,types.some(t=>t[0]===old)?old:types[0][0]);
 const cols=[...new Set(state.filtered.flatMap(Object.keys))].filter(k=>state.filtered.some(r=>r[k]!==''&&r[k]!=null&&typeof r[k]!=='boolean'&&Number.isFinite(Number(r[k]))));options($('plot-x'),cols.map(c=>[c,human(c)]));options($('plot-y'),cols.map(c=>[c,human(c)]),cols[1]||cols[0]);plotControls();
}
function plotControls(){const t=$('plot-type').value;$('x-wrap').hidden=!['scatter','histogram'].includes(t);$('y-wrap').hidden=t!=='scatter';$('plot-limit').parentElement.hidden=!['forest','heatmap'].includes(t);$('draw').disabled=!state.filtered.length;}
function view(which){state.view=which;$('table-panel').hidden=which!=='table';$('plot-panel').hidden=which!=='plot';$('view-table').setAttribute('aria-selected',which==='table');$('view-plot').setAttribute('aria-selected',which==='plot');}
let busy=false;
async function drawPlot(){
 if(busy)return;busy=true;$('draw').disabled=true;const id=++state.plotId;
 status('Preparing the figure…');
 if(!state.worker){state.worker=new Worker('./plot-worker.js');state.worker.onmessage=({data})=>{if(data.status){if(busy)status(data.status);return;}busy=false;$('draw').disabled=!state.filtered.length;if(data.id!==state.plotId)return;if(data.error){status(data.error.split('\n').filter(Boolean).slice(-1)[0],true);return;}state.figure=data.result;const img=el('img',undefined,{alt:state.dataset.title+' — '+$('plot-type').selectedOptions[0].text});const url=URL.createObjectURL(new Blob([data.result.svg],{type:'image/svg+xml'}));img.onload=()=>URL.revokeObjectURL(url);img.src=url;$('figure-container').replaceChildren(img);$('plot-caption').textContent=data.result.caption;$('svg').disabled=$('png').disabled=false;status('Figure ready.');};state.worker.onerror=()=>{busy=false;$('draw').disabled=false;status('The plotting worker could not start. Reload the page and try again.',true);state.worker.terminate();state.worker=null;};}
 state.worker.postMessage({id,request:{rows:state.filtered,type:$('plot-type').value,x:$('plot-x').value,y:$('plot-y').value,limit:Number($('plot-limit').value),title:state.dataset.title}});
}
$('dataset').onchange=loadDataset;$('search').oninput=applyFilters;$('include-failed').onchange=applyFilters;$('reset').onclick=()=>{$('search').value='';$('include-failed').checked=false;buildFacets();applyFilters();};$('prev').onclick=()=>{state.page--;drawTable();};$('next').onclick=()=>{state.page++;drawTable();};$('view-table').onclick=()=>view('table');$('view-plot').onclick=()=>view('plot');$('plot-type').onchange=()=>{clearFigure();plotControls();};$('draw').onclick=drawPlot;
$('download').onclick=()=>save(filteredCSV(state.filtered),state.dataset.id.replaceAll('/','_')+'.csv','text/csv;charset=utf-8');$('svg').onclick=()=>save(state.figure.svg,'atlas_figure.svg','image/svg+xml');$('png').onclick=()=>{const b=Uint8Array.from(atob(state.figure.png),c=>c.charCodeAt(0));save(new Blob([b],{type:'image/png'}),'atlas_figure.png');};
for(const [i,[id,name]]of sections.entries()){const b=el('button',undefined,{role:'tab','aria-selected':id===state.section});b.dataset.section=id;b.append(el('span',String(i+1).padStart(2,'0'),{class:'num'}),el('span',name));b.onclick=()=>changeSection(id);$('sections').append(b);}
(async()=>{try{state.catalog=await api('atlas_web_datasets',{select:'id,section,title,kind,description,run_key,row_count,provisional',order:'id.asc',limit:200});if(!state.catalog.length){$('connection').textContent='Website ready · data release pending';$('dataset-title').textContent='Collections awaiting publication';status('The results have been prepared. Public access will open after the lab approves the data release.');$('search').disabled=true;$('reset').disabled=true;return;}$('connection').textContent=state.catalog.length+' collections · connected to the research data';const p=new URLSearchParams(location.hash.slice(1));const d=state.catalog.find(d=>d.id===p.get('dataset'));changeSection(d?.section||'rna',d?.id||'immune_coarse/effects');}catch(e){$('connection').textContent='Connection unavailable';$('dataset-title').textContent='Results temporarily unavailable';status(e.message,true);}})();
