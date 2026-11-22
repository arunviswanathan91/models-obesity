const groupNames = ['Normal weight', 'Overweight', 'Obesity'];
const fmt = v => Number(v.toPrecision(3)).toString();
const element = (tag, text) => {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  return node;
};
function chart(rows, continuous, upper, title) {
  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', '0 0 450 260');
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label', title);
  const add = (tag, attrs, text) => {
    const item = document.createElementNS(ns, tag);
    Object.entries(attrs).forEach(([key, value]) => item.setAttribute(key, value));
    if (text !== undefined) item.textContent = text;
    svg.append(item);
    return item;
  };
  const y = v => 194 - 151 * v / upper;
  const left = 68, right = 424;
  const x = (r, i) => continuous
    ? left + 8 + (right-left-16) * (r.bmi-rows[0].bmi)/(rows.at(-1).bmi-rows[0].bmi)
    : left + (right-left)*(i+0.5)/rows.length;
  const label = (px, py, value, anchor='middle') =>
    add('text', {x:px,y:py,'text-anchor':anchor,'font-size':12,fill:'#374151'}, value);
  label(left, 21, 'Estimated RNA fraction (%)', 'start');
  for (let i=0;i<=4;i++) {
    const v=upper*i/4;
    add('line',{x1:left,x2:right,y1:y(v),y2:y(v),stroke:'#dde1e4'});
    label(left-7,y(v)+4,fmt(100*v),'end');
  }
  rows.forEach((r,i)=>{
    const px=x(r,i);
    add('line',{x1:px,x2:px,y1:y(r.low),y2:y(r.high),stroke:'#196e73','stroke-width':3});
    for(const v of [r.low,r.high]) add('line',{x1:px-5,x2:px+5,y1:y(v),y2:y(v),stroke:'#196e73'});
    add('circle',{cx:px,cy:y(r.median),r:5,fill:'#16595d'});
    if(!continuous){
      label(px,216,groupNames[i]);
      label(px,235,'n = '+r.n);
    }
  });
  if(continuous){
    const start=rows[0].bmi, end=rows.at(-1).bmi;
    for(let tick=Math.ceil(start/2)*2;tick<=end;tick+=2) label(x({bmi:tick}),216,tick);
    label((left+right)/2,242,'BMI (kg/m²)');
  }
  return svg;
}
export function renderTheta(container, rows) {
  const select=element('select');
  const picker=element('label','Reference and component ');
  picker.append(select);
  for(const family of [...new Set(rows.map(r=>r.family))]){
    const group=element('optgroup');
    group.label=family;
    for(const row of rows.filter(r=>r.family===family)){
      const option=element('option',row.component.replaceAll('_',' '));
      option.value=String(rows.indexOf(row));
      group.append(option);
    }
    select.append(group);
  }
  select.value=String(Math.max(0,rows.findIndex(r=>r.component==='TAM')));
  const panels=element('div'); panels.className='theta-panels';
  const details=element('p');details.className='help';
  details.textContent='Both tests use individual tumors and are unadjusted for covariates. BH-adjusted P values account for the 47 components separately for each test. Fractions are relative to the selected deconvolution reference.';
  const draw=()=>{
    const r=rows[Number(select.value)];
    panels.replaceChildren();
    const upper=Math.max(...[...r.trend,...r.groups].map(b=>b.high))*1.08;
    for(const continuous of [true,false]){
      const panel=element('section');
      const title=continuous?'Continuous BMI':'BMI groups';
      const stats=continuous
        ? 'Spearman ρ = '+fmt(r.rho)+' · P = '+fmt(r.continuous_p)+' · BH P = '+fmt(r.continuous_q)
        : 'Kruskal–Wallis H(2) = '+fmt(r.h)+' · P = '+fmt(r.group_p)+' · BH P = '+fmt(r.group_q);
      panel.append(element('h3',title),element('p',stats),
        chart(continuous?r.trend:r.groups,continuous,upper,title+' and '+r.component+' RNA fraction'),
        element('p',continuous?'Median and interquartile range within eight BMI quantile bins, positioned at median BMI. Spearman correlation uses all '+r.n+' tumors.':
          'Median and interquartile range by BMI group (18.5–<25, 25–<30 and ≥30 kg/m²).'));
      panels.append(panel);
    }
  };
  select.onchange=draw;
  container.replaceChildren(picker,panels,details);
  draw();
}
