const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const source=fs.readFileSync(__dirname+'/../web/app.js','utf8');
const clear=source.slice(source.indexOf('function clearFigure()'),source.indexOf('const runLabel='));
const flow=source.slice(source.indexOf('function schedulePlot()'),source.indexOf('function showPCA()'));
const nodes=new Map(),timers=new Map(),workers=[];let nextTimer=0;
function node(id){if(!nodes.has(id))nodes.set(id,{hidden:false,disabled:false,value:id==='plot-limit'?'20':'forest',selectedOptions:[{text:'Effect intervals'}],replaceChildren(x){this.content=x;}});return nodes.get(id);}
const context=vm.createContext({$:node,state:{view:'plot',dataset:{id:'test/effects',title:'Test'},filtered:[{estimate:1}],rows:[],plotId:0,worker:null},autoPlotTimer:null,PROGRAMME_LABELS:{},programmeLabel:r=>r.feature,plotTypes:()=>[['forest','Effect intervals']],status:()=>{},el:(tag,text)=>({tag,text}),URL:{createObjectURL:()=> 'blob:test',revokeObjectURL:()=>{}},Blob:class {},setTimeout:fn=>{timers.set(++nextTimer,fn);return nextTimer;},clearTimeout:id=>timers.delete(id),Worker:class{constructor(){this.jobs=[];workers.push(this);}postMessage(job){this.jobs.push(job);}terminate(){}},updateModeUI:()=>{},changeSection:()=>{}});
vm.runInContext(clear+flow,context);
function flush(){for(const [id,fn] of [...timers]){timers.delete(id);fn();}}
context.schedulePlot();context.schedulePlot();context.schedulePlot();assert.equal(timers.size,1);flush();assert.equal(workers[0].jobs.length,1);
const first=workers[0].jobs[0];
context.state.filtered=[{estimate:2}];context.schedulePlot();flush();
context.state.filtered=[{estimate:3}];context.schedulePlot();flush();
assert.equal(workers[0].jobs.length,1,'only one active render');
workers[0].onmessage({data:{id:first.id,result:{svg:'old',caption:'old'}}});
assert.equal(context.state.figure,null,'stale result must not display');assert.equal(workers[0].jobs.length,2);
const latest=workers[0].jobs[1];assert.equal(latest.request.rows[0].estimate,3);
workers[0].onmessage({data:{id:latest.id,result:{svg:'latest',caption:'latest'}}});assert.equal(context.state.figure.caption,'latest');
context.schedulePlot();context.state.view='table';context.schedulePlot();flush();assert.equal(workers[0].jobs.length,2,'hidden plot must not render');
context.state.view='plot';context.state.filtered=[];context.schedulePlot();flush();assert.equal(workers[0].jobs.length,2,'empty selection must not render');
context.state.filtered=[{estimate:4}];context.schedulePlot();flush();const failed=workers[0].jobs.at(-1);workers[0].onmessage({data:{id:failed.id,error:'test failure'}});assert.equal(node('draw').hidden,false,'retry shown after rendering failure');
console.log('PASS: debounce, latest-only queue, stale results, table/empty guards and retry');
