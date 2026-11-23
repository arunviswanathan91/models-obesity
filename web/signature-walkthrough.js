// In-page replay of the signature-curation workflow. Scripted: no model is called.
// Rules mirror analysis/signatures/curate_signatures.py: overlap coefficient on
// uppercased symbols, pairs above 0.50 reviewed, proposals above 0.20 rejected.
const DEDUP=0.5, NEW_SIG=0.2;
const MENU={1:'Keep the first only',2:'Keep the second only',3:'Keep both',4:'Merge into the first',5:'Merge into the second'};
const up=g=>new Set(g.map(x=>String(x).toUpperCase()));
const overlap=(a,b)=>a.size&&b.size?[...a].filter(x=>b.has(x)).length/Math.min(a.size,b.size):0;
const h=(tag,cls,text)=>{const e=document.createElement(tag);if(cls)e.className=cls;if(text!==undefined)e.textContent=text;return e;};
const chips=(genes,cls='')=>{const w=h('div','sw-chips');for(const g of genes)w.append(h('span','sw-chip '+cls,g));if(!genes.length)w.append(h('span','sw-none','none, it is a subset'));return w;};

export function buildSteps(run){
  const steps=[];
  for(const [i,r] of run.dedup_replies.entries()){
    const sigs=run.catalogue[r.cell_type];const a=up(sigs[r.sig1]),b=up(sigs[r.sig2]);
    const key=`${r.cell_type}|${r.sig1}|${r.sig2}`;
    steps.push({kind:'pair',title:`${r.cell_type}: ${r.sig1} and ${r.sig2}`,ov:overlap(a,b),limit:DEDUP,
      shared:[...a].filter(x=>b.has(x)).sort(),only1:[...a].filter(x=>!b.has(x)).sort(),only2:[...b].filter(x=>!a.has(x)).sort(),
      reply:r.reasoning,suggested:r.choice,decision:run.dedup_decisions[key]});
  }
  for(const [cell,list] of Object.entries(run.suggestions)){
    const existing=new Set(Object.values(run.catalogue[cell]).flat().map(x=>x.toUpperCase()));
    for(const s of list){
      const g=up(s.genes);const inter=[...g].filter(x=>existing.has(x)).length;const ov=g.size?inter/g.size:0;
      steps.push({kind:'new',title:`${cell}: proposed ${s.name}`,ov,limit:NEW_SIG,genes:s.genes,reply:s.rationale,
        rejected:ov>NEW_SIG,short:s.genes.length<8,decision:run.suggestion_decisions[`${cell}|${s.name}`]});
    }
  }
  return steps;
}

function outcome(s){
  if(s.kind==='pair'){const d=s.decision;return {ok:true,text:`Person chose: ${MENU[d]}.`};}
  if(s.rejected)return {ok:false,text:`Rejected automatically. Overlap ${s.ov.toFixed(2)} is above ${s.limit.toFixed(2)}, so a person never reviews it.`};
  if(s.short)return {ok:false,text:'Proposal is shorter than the minimum gene count, so the person declined it.'};
  return {ok:s.decision==='y',text:s.decision==='y'?'Person accepted it. Added to the catalogue.':'Person declined it.'};
}

export async function mountWalkthrough(root,{onBack,url='./signature-toy-run.json'}){
  root.replaceChildren();
  const back=h('button','secondary sw-back','← Back to gene signatures');back.type='button';back.onclick=onBack;
  const head=h('div','sw-head');head.append(h('h3',undefined,'How the gene signatures were curated'));
  const note=h('p','sw-banner');note.append(h('strong',undefined,'Scripted replay. '),document.createTextNode('A small toy catalogue and fixed replies written in advance. Nothing here calls a model or needs a key, and the gene lists are illustrative, not the manuscript catalogue.'));
  const intro=h('p','help','The real workflow asks a language model to flag near-duplicate signatures and to propose missing ones. Two rules sit around it: pairs that overlap by more than 0.50 are reviewed by a person, and proposed signatures that overlap the existing catalogue by more than 0.20 are rejected before anyone sees them. Every change needs a person’s decision.');
  const stage=h('div','sw-stage');stage.setAttribute('aria-live','polite');
  const nav=h('div','sw-nav');
  root.append(back,head,note,intro,stage,nav);
  let steps;try{steps=buildSteps(await (await fetch(url)).json());}catch{stage.textContent='The demonstration data could not be loaded.';return;}
  let i=0;
  const prev=h('button','secondary','Previous'),next=h('button','primary','Next step'),count=h('span','help sw-count');
  prev.type=next.type='button';nav.append(prev,next,count);
  const show=()=>{
    const s=steps[i];stage.replaceChildren();
    stage.append(h('h4',undefined,s.title));
    const m=h('div','sw-meter');const bar=h('i');bar.style.width=Math.min(100,s.ov*100)+'%';bar.className=s.ov>s.limit?'over':'under';
    const tick=h('u');tick.style.left=s.limit*100+'%';m.append(bar,tick);stage.append(m);
    stage.append(h('p','help',`Overlap coefficient ${s.ov.toFixed(2)}. Limit ${s.limit.toFixed(2)}. `+(s.kind==='pair'?'Pairs above the limit go to review.':'Proposals above the limit are rejected. Exactly 0.20 passes.')));
    if(s.kind==='pair'){
      stage.append(h('p','sw-lab',`Shared (${s.shared.length})`),chips(s.shared,'shared'),h('p','sw-lab','Only in the first'),chips(s.only1),h('p','sw-lab','Only in the second'),chips(s.only2));
    }else stage.append(h('p','sw-lab',`Proposed genes (${s.genes.length})`),chips(s.genes));
    const say=h('div','sw-say');say.append(h('em',undefined,'Scripted reply, no model called'),h('span',undefined,s.reply));stage.append(say);
    const o=outcome(s);stage.append(h('p',o.ok?'sw-ok':'sw-no',o.text));
    prev.disabled=i===0;next.disabled=i===steps.length-1;count.textContent=`Step ${i+1} of ${steps.length}`;
  };
  prev.onclick=()=>{i--;show();};next.onclick=()=>{i++;show();};show();
}
