// Checks the cellular-model picker: which collections each model and mode offers,
// and that a dataset id sets the controls back correctly.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const source=fs.readFileSync(__dirname+'/../web/app.js','utf8');
const slice=(a,b)=>source.slice(source.indexOf(a),source.indexOf(b));
const code=slice('const runLabel=','function brief(d)')+slice('function syncRunControls','function changeSection');
const controls={'cellular-model':{value:'st'},'bmi-model':{value:'continuous'},compartment:{value:'immune_fine'},'collection-search':{value:''}};
const ds=(id,kind,run_key,cellular_model)=>({id,section:'rna',title:id,kind,run_key,cellular_model});
const catalog=[
 ds('immune_fine/effects','effects','immune_fine','orthogonal_patient_effect'),
 ds('immune_fine/run_diagnostics','diagnostics','immune_fine','orthogonal_patient_effect'),
 ds('st_immune_fine/effects','effects','st_immune_fine','student_t_intercept'),
 ds('st_immune_fine/run_diagnostics','diagnostics','st_immune_fine','student_t_intercept'),
 ds('st_immune_fine/sensitivity_summary','table','st_immune_fine','student_t_intercept'),
 ds('st_categorical_immune_fine/effects','effects','st_categorical_immune_fine','student_t_intercept'),
 ds('compare_immune_fine/features','comparison','compare_immune_fine','comparison'),
 ds('compare_immune_fine/celltypes','comparison','compare_immune_fine','comparison'),
 ds('compare_rna/summary','comparison','compare_rna','comparison'),
 ds('rna/compartment_diagnostics','diagnostics','overview',null)];
const context=vm.createContext({$:id=>controls[id],state:{catalog,mode:'biology'}});
vm.runInContext(code+';this.api={availableCollections,syncRunControls,currentRun,runParts,preferredRNA};',context);
const {availableCollections,syncRunControls,currentRun,runParts,preferredRNA}=context.api;
const ids=()=>availableCollections('rna').map(d=>d.id).sort();

assert.deepEqual(ids(),['st_immune_fine/effects'],'Student-t biology view');
controls['cellular-model'].value='orth';
assert.deepEqual(ids(),['immune_fine/effects'],'orthogonal biology view');
controls['cellular-model'].value='compare';
assert.deepEqual(ids(),['compare_immune_fine/celltypes','compare_immune_fine/features','compare_rna/summary'],'comparison view');
context.state.mode='reliability';
assert.deepEqual(ids(),['immune_fine/run_diagnostics','rna/compartment_diagnostics','st_immune_fine/run_diagnostics'],'comparison reliability shows both models');
controls['cellular-model'].value='st';
assert.deepEqual(ids(),['rna/compartment_diagnostics','st_immune_fine/run_diagnostics'],'Student-t reliability');
controls['bmi-model'].value='categorical';context.state.mode='biology';
assert.equal(currentRun(),'st_categorical_immune_fine');
assert.deepEqual(ids(),['st_categorical_immune_fine/effects']);

for(const [id,model,bmi,comp] of [['immune_fine/effects','orth','continuous','immune_fine'],['st_categorical_immune_fine/effects','st','categorical','immune_fine'],['compare_immune_fine/features','compare','continuous','immune_fine']]){
 syncRunControls(catalog.find(d=>d.id===id));
 assert.deepEqual([controls['cellular-model'].value,controls['bmi-model'].value,controls.compartment.value],[model,bmi,comp],id);
}
syncRunControls(catalog.find(d=>d.id==='compare_rna/summary'));assert.equal(controls['cellular-model'].value,'compare');
assert.deepEqual({...runParts('st_categorical_nonimmune')},{model:'st',base:'categorical_nonimmune',categorical:true,compartment:'nonimmune'});
assert.equal(preferredRNA(),'immune_coarse/effects','falls back when the Student-t collection is not published');
console.log('PASS: model picker collections, reliability view, control sync and fallback');
