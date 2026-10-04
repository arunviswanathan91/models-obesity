# Heatmap renderer copied verbatim from the original BiB figure report.
# Browser adapter supplies the same columns; no replacement heatmap algorithm.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns
import io, json, base64
BLUE, MAGENTA, TEAL = "#2C5AA0", "#D81B60", "#0E7C7B"
DIVERGING = mcolors.LinearSegmentedColormap.from_list("bmi_diverging", [BLUE, "#FFFFFF", MAGENTA])
EMPTY_CELL = "#FAF9F6"
CELL_INCH = 0.14
HM_TICK, HM_LABEL, HM_TITLE, HM_SPINE = 8, 10, 11.5, 1.0
SHARED_LABEL = "Shared across cell types"
_MEASURE = {}
_CAPTURED = None
plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],
 "font.size":9.5,"axes.spines.top":False,"axes.spines.right":False,"axes.edgecolor":"#1F1F1F",
 "grid.color":"#DCDDE0","grid.linewidth":0.45,"svg.fonttype":"none","figure.facecolor":"white","savefig.facecolor":"white"})
def export_figure(fig, stem, section):
    global _CAPTURED
    _CAPTURED = fig

def text_width(strings, size, weight="normal"):
    """Rendered width in inches of the widest string."""
    if "fig" not in _MEASURE:
        _MEASURE["fig"] = plt.figure(figsize=(1, 1))
        _MEASURE["renderer"] = _MEASURE["fig"].canvas.get_renderer()
    fig, renderer = _MEASURE["fig"], _MEASURE["renderer"]
    widest = 0.0
    for text in strings:
        artist = fig.text(0, 0, str(text), fontsize=size, fontweight=weight)
        widest = max(widest, artist.get_window_extent(renderer).width / fig.dpi)
        artist.remove()
    return widest

def draw_heatmap(features, row_names, columns, palette, home, vmax, stem, section, title, value_col):
    subset = features[features["signature_label"].isin(row_names)]
    means = subset.pivot_table(index="signature_label", columns="celltype_label", values=value_col, aggfunc="first")
    flags = subset.pivot_table(index="signature_label", columns="celltype_label", values="hdi_95_excludes_zero",
                               aggfunc="first")
    means = means.reindex(index=row_names, columns=columns)
    flags = flags.reindex_like(means).astype(float).fillna(0).astype(bool)
    n_rows, n_cols = means.shape
    left = text_width(row_names, HM_TICK, "bold") + 0.42
    bottom = text_width(columns, HM_TICK) * 0.71 + 0.5
    top = 0.42
    axes_w, axes_h = CELL_INCH * n_cols, CELL_INCH * n_rows
    width, height = left + axes_w + 0.1, bottom + axes_h + top
    fig = plt.figure(figsize=(max(width, 2.4), height))
    ax = fig.add_axes([left / max(width, 2.4), bottom / height, axes_w / max(width, 2.4), axes_h / height])
    sns.heatmap(means, cmap=DIVERGING, center=0, vmin=-vmax, vmax=vmax, ax=ax, cbar=False,
                linewidths=0.1, linecolor="gray", square=False, mask=means.isna())
    ax.set_facecolor(EMPTY_CELL)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_xticklabels(columns, rotation=45, ha="right", rotation_mode="anchor", fontsize=HM_TICK)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_yticklabels(row_names, rotation=0, fontsize=HM_TICK, fontweight="bold")
    for label in ax.get_yticklabels():
        label.set_color(palette.get(home.get(label.get_text(), SHARED_LABEL), "black"))
    ax.set_xticks(np.arange(n_cols + 1), minor=True)
    ax.set_yticks(np.arange(n_rows + 1), minor=True)
    ax.grid(which="minor", color="gray", linewidth=0.3)
    ax.tick_params(which="both", length=0, pad=2)
    values = means.to_numpy(float)
    for r in range(n_rows):
        for c in range(n_cols):
            if np.isfinite(values[r, c]) and flags.iat[r, c]:
                ax.add_patch(mpatches.Circle((c + 0.5, r + 0.5), 0.2, facecolor="white", edgecolor="black",
                                             linewidth=0.6, zorder=10))
    ax.set_xlabel("Cell type", fontsize=HM_LABEL, fontweight="bold", labelpad=4)
    ax.set_ylabel("Signature", fontsize=HM_LABEL, fontweight="bold", labelpad=4)
    ax.set_title(title, fontsize=HM_TITLE, fontweight="bold", pad=6, loc="right")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(HM_SPINE)
        spine.set_edgecolor("black")
    export_figure(fig, stem, section)

# Adapters for the supplied manuscript Python plotting code.
# Figures are computed from result rows; no pre-rendered assets are fetched.
import matplotlib as mpl
from matplotlib.lines import Line2D
import textwrap, re
BLUE, LIGHT_BLUE, MAGENTA, LIGHT_MAGENTA = "#2C5AA0", "#A9BEDF", "#D81B60", "#F2B6CB"
TEAL, LIGHT_TEAL, MID_TEAL = "#0E7C7B", "#BFE3E2", "#5FAFAE"
INK, MID_GREY, GRID_GREY = "#1F1F1F", "#6F7378", "#DCDDE0"
DIVERGING = mpl.colors.LinearSegmentedColormap.from_list("d", [BLUE, "white", MAGENTA])
SIZE = {"base": 9.5, "tick": 9, "ytick": 10.5, "label": 9.5, "title": 11.5, "legend": 10, "value": 8.5,
        "axis_line": 0.9, "grid_line": 0.45, "data_line": 1.6, "interval_line": 1.9, "zero_line": 1.0, "marker": 32}
MANUSCRIPT_STYLE = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": SIZE["base"], "axes.titlesize": SIZE["title"], "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.titlepad": 6, "axes.labelsize": SIZE["label"], "axes.labelpad": 4,
    "axes.edgecolor": INK, "axes.linewidth": SIZE["axis_line"], "axes.spines.top": False,
    "axes.spines.right": False, "axes.facecolor": "white", "xtick.labelsize": SIZE["tick"],
    "ytick.labelsize": SIZE["ytick"], "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.color": INK, "ytick.color": INK,
    "legend.fontsize": SIZE["legend"], "legend.frameon": False, "lines.linewidth": SIZE["data_line"],
    "grid.color": GRID_GREY, "grid.linewidth": SIZE["grid_line"], "figure.facecolor": "white",
    "savefig.facecolor": "white", "svg.fonttype": "none"}


def panel(w, h, left, right=0.14, bottom=0.62, top=0.36):
    fig = plt.figure(figsize=(left + w + right, bottom + h + top))
    ax = fig.add_axes([left / (left + w + right), bottom / (bottom + h + top),
                       w / (left + w + right), h / (bottom + h + top)])
    return fig, ax


def style_axis(ax, grid="x"):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(SIZE["axis_line"])
    if grid:
        ax.grid(axis=grid, zorder=0)
        ax.set_axisbelow(True)


def set_title(ax, text_):
    fig = ax.figure
    fig_w = fig.get_size_inches()[0]
    pos = ax.get_position()
    left_in = pos.x0 * fig_w
    axes_w = pos.width * fig_w
    title = ax.set_title(text_, fontsize=SIZE["title"], fontweight="bold", loc="left")
    title.set_x(-left_in / axes_w + 0.005)


SCENARIOS=['null','protein_only','conditional_signal','heavy_tails','sbc']
SCENARIO_LABELS={'null':'Null','protein_only':'Protein only','conditional_signal':'Conditional signal','heavy_tails':'Heavy tails','sbc':'SBC (prior predictive)'}
QUANTITIES={'beta_protein':'Protein BMI slope','beta_ptm':'Marginal PTM BMI slope','beta_conditional':'Protein-conditioned BMI slope','rho':'Residual correlation','coupling':'Protein-to-PTM coupling','nu':'Student-t degrees of freedom'}
VARIANTS=['primary','fixed_plex','stage_subset_base','stage_adjusted','wes_subset_base','wes_adjusted','mutation_adjusted','exclude_recorded_weight_loss','exclude_adenosquamous']
VLAB={'primary':'Primary','fixed_plex':'Fixed plex','stage_subset_base':'Stage subset','stage_adjusted':'Stage adjusted','wes_subset_base':'Purity subset','wes_adjusted':'Purity adjusted','mutation_adjusted':'Mutation adjusted','exclude_recorded_weight_loss':'Exclude weight loss','exclude_adenosquamous':'Exclude adenosquamous'}
COMPARISONS=[('primary','fixed_plex'),('primary','stage_subset_base'),('stage_subset_base','stage_adjusted'),('primary','wes_subset_base'),('wes_subset_base','wes_adjusted'),('primary','mutation_adjusted'),('primary','exclude_recorded_weight_loss'),('primary','exclude_adenosquamous')]
CLAB=['Fixed plex','Stage subset restriction','Stage adjustment (matched)','Purity subset restriction','Purity adjustment (matched)','Mutation adjustment','Exclude recorded weight loss','Exclude adenosquamous']


def feature_labels(df, gene='gene', feature='feature_id'):
    # Keep peptide sequence where needed to distinguish equal glycan compositions.
    def label(r):
        f=re.sub(r'\[[^]]*\]', '', str(r[feature]))
        f=re.sub(r'^n(?=[A-Z])', '', f)
        f=re.sub(r'([NHFSG])0(?=[NHFSG]|$)', '', f)
        if '_' in f and re.search(r'_[STY]\d+$',f): f=f.rsplit('_',1)[-1]
        return str(r[gene])+' | '+f
    return df.apply(label,axis=1).tolist()

def manuscript_plot(df, req):
    kind=req['type']
    with mpl.rc_context(MANUSCRIPT_STYLE):
        return _manuscript_plot(df,req,kind)

def _manuscript_plot(df,req,kind):
    if kind in ('ptm_primary','ptm_colourbar'):
        scale=pd.DataFrame(req.get('scale_rows') or req['rows'])
        vmax=max(float(np.abs(scale[['posterior_mean_marginal','posterior_mean_conditional']].to_numpy()).max()),1e-8)
        if kind=='ptm_colourbar':
            fig=plt.figure(figsize=(1.25,3)); ax=fig.add_axes([.12,.05,.16,.9])
            bar=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-vmax,vmax),cmap=DIVERGING),cax=ax)
            bar.outline.set_linewidth(1);bar.ax.tick_params(length=3,width=.8,labelsize=9)
            bar.ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5,symmetric=True))
            bar.set_label('Posterior mean (response SD per\n5 kg/m² higher BMI, within plex)',fontsize=9.5,fontweight='bold',labelpad=6)
        else:
            labels=df.short_label.tolist()
            if len(set(labels))!=len(labels):labels=feature_labels(df)
            if len(set(labels))!=len(labels):raise ValueError('Feature labels are not unique; narrow the selection.')
            rows=[]
            for (_,r),label in zip(df.iterrows(),labels):
                for suffix,col in [('marginal','Marginal PTM'),('conditional','Protein-conditioned PTM')]:
                    rows.append(dict(signature_label=label,celltype_label=col,value=r['posterior_mean_'+suffix],hdi_95_excludes_zero=r['hdi_95_lower_'+suffix]>0 or r['hdi_95_upper_'+suffix]<0))
            draw_heatmap(pd.DataFrame(rows),labels,['Marginal PTM','Protein-conditioned PTM'],{},{},vmax,'browser','extended','Supported BMI effects','value')
            fig=_CAPTURED;fig.axes[0].set_xlabel('BMI association',fontsize=HM_LABEL,fontweight='bold');fig.axes[0].set_ylabel('PTM feature',fontsize=HM_LABEL,fontweight='bold')
        return fig, f'Colour indicates the posterior mean BMI effect (response SD per 5 kg/m²), on a symmetric scale from −{vmax:.3g} to +{vmax:.3g}. White circles indicate pointwise 95% HDIs excluding zero.'
    if kind=='manuscript_heatmap':
        if df.duplicated(['signature_label','celltype_label']).any():raise ValueError('Duplicate heatmap cells; narrow the selection.')
        rows=list(dict.fromkeys(df.signature_label));cols=list(dict.fromkeys(df.celltype_label));vmax=max(float(df.value.abs().max()),1e-8)
        draw_heatmap(df,rows,cols,{},{},vmax,'browser','extended','','value');fig=_CAPTURED
        return fig,f'Colour indicates the posterior mean effect, on a symmetric scale from −{vmax:.3g} to +{vmax:.3g}. White circles indicate pointwise 95% HDIs excluding zero. Blank cells indicate unavailable estimates.'
    if kind=='sensitivity_grid':
        if df.quantity.nunique()!=1:raise ValueError('Select one model quantity.')
        identities=df[['key','gene','feature_id']].drop_duplicates().sort_values(['gene','feature_id']);order=identities.key.tolist();labels=feature_labels(identities)
        variants=[v for v in VARIANTS if v in set(df.variant)]
        if df.duplicated(['key','variant']).any():raise ValueError('Select one quantity and contrast.')
        matrix=df.pivot(index='key',columns='variant',values='posterior_mean').reindex(index=order,columns=variants)
        good=df.pivot(index='key',columns='variant',values='numerical_checks_pass').reindex(index=order,columns=variants)
        resolved=df.pivot(index='key',columns='variant',values='hdi_excludes_zero').reindex(index=order,columns=variants)
        vmax=float(np.nanmax(np.abs(matrix.to_numpy()))) or .01
        f,ax=panel(4.9,.245*len(order),left=text_width(labels,9)+.2,bottom=1.55,top=.42,right=.15)
        ax.imshow(matrix,aspect='auto',cmap=DIVERGING,vmin=-vmax,vmax=vmax)
        for i in range(len(order)):
            for j in range(len(variants)):
                if pd.isna(matrix.iloc[i,j]):continue
                if not good.iloc[i,j]:ax.text(j,i,'×',ha='center',va='center',color=INK)
                elif resolved.iloc[i,j]:ax.plot(j,i,'o',ms=3.5,mfc='none',mec=INK,mew=.65)
        ax.set_xticks(range(len(variants)),[VLAB[v] for v in variants],rotation=50,ha='right',fontsize=9)
        ax.set_yticks(range(len(order)),labels);ax.tick_params(length=0)
        for i in range(1,len(identities)):
            if identities.gene.iloc[i]!=identities.gene.iloc[i-1]:ax.axhline(i-.5,color='white',lw=1.3)
        set_title(ax,'Expanded PTM sensitivities')
        return f,f'Posterior mean effects across model variants. Colour scale ±{vmax:.3g}; '+str(df.effect_units.iloc[0])+'. Open circles: pointwise 95% HDIs exclude zero. Subset restriction and covariate adjustment are separate analyses. No multiplicity control.'
    if kind=='sensitivity_comparison':
        if df.comparison.nunique()!=1:raise ValueError('Select one matched comparison.')
        x=df.sort_values(['gene_base','feature_id_base']).reset_index(drop=True);order=x.key.tolist();labels=feature_labels(x,'gene_base','feature_id_base')
        f,ax=panel(4.65,.26*len(order),left=text_width(labels,9)+.2,bottom=.65,top=.62,right=.45)
        for i,r in x.iterrows():
            if i%2==0:ax.axhspan(i-.5,i+.5,color='#F6F6F6',zorder=0)
            for suffix,offset,color in [('base',-.14,BLUE),('alt',.14,MAGENTA)]:
                lo,hi,mu=r['hdi_95_lower_'+suffix],r['hdi_95_upper_'+suffix],r['posterior_mean_'+suffix]
                adequate=bool(r['numerical_checks_pass_'+suffix]);color=color if adequate else MID_GREY
                ax.hlines(i+offset,lo,hi,color=color,lw=SIZE['interval_line'],zorder=2)
                ax.plot(mu,i+offset,marker='o' if adequate else 'x',ms=4,mec=color,mfc=color if (lo>0 or hi<0) else 'white',ls='',zorder=3)
        ax.axvline(0,color=MID_GREY,lw=1,ls='--');ax.set_yticks(range(len(order)),labels);ax.invert_yaxis();ax.set_xlabel('Protein-conditioned BMI slope (SD per 5 kg/m²)');style_axis(ax)
        set_title(ax,str(x.comparison.iloc[0])+'\n'+VLAB[x.variant_base.iloc[0]]+' versus '+VLAB[x.variant_alt.iloc[0]])
        return f,'Blue: reference; pink: sensitivity. Filled/open points exclude/include zero in pointwise 95% HDIs. Stage and purity adjustment use matched subset baselines. Classification differences are descriptive, not posterior contrast tests.'
    if kind=='sensitivity_classification':
        t=df;ids=range(int(t.total.max()))
        f,ax=panel(4.6,2.55,left=2.65,bottom=.62,top=.42,right=.2);y=np.arange(len(t))
        ax.barh(y,t.same_class,color=BLUE,label='Same interval class');ax.barh(y,t.changed_class,left=t.same_class,color=MAGENTA,label='Changed interval class')
        ax.barh(y,t.numerical_failures,left=t.adequate,color=GRID_GREY,label='Numerical failure')
        ax.set_yticks(y,t.comparison,fontsize=9);ax.invert_yaxis();ax.set_xlabel('Selected modified features');ax.set_xlim(0,len(ids));style_axis(ax);set_title(ax,'Interval classification across sensitivity comparisons')
        return f,'Same (blue) versus changed (pink) interval classification within the selected features; gray denotes numerical failures. Descriptive counts, not discovery rates.'
    if kind.startswith('calibration_') and kind!='calibration_completion':
        metric=kind.removeprefix('calibration_');summary=df;scope_colors={'All completed':BLUE,'Numerically adequate':MAGENTA}
        if 'template' in df and df.template.nunique()!=1:raise ValueError('Select one calibration template.')
        if df.duplicated(['quantity','scope','scenario']).any():raise ValueError('Repeated calibration summaries; narrow the selection.')
        f,axs=plt.subplots(3,2,figsize=(7.4,8.6));f.subplots_adjust(left=.22,right=.98,bottom=.065,top=.94,wspace=.30,hspace=.60)
        for ax,(q,title) in zip(axs.flat,QUANTITIES.items()):
            for scope,offset in [('All completed',-.12),('Numerically adequate',.12)]:
                g=summary[(summary.scope==scope)&(summary.quantity==q)].set_index('scenario').reindex(SCENARIOS);xx=g[metric].to_numpy();yy=np.arange(5)+offset
                if metric=='coverage':err=np.vstack([xx-g.coverage_low,g.coverage_high-xx])
                elif metric=='bias':err=1.96*g.bias_mcse.to_numpy()
                else:err=None
                ax.errorbar(xx,yy,xerr=err,fmt='o',ms=4,color=scope_colors[scope],lw=1.4,capsize=2)
            ax.set_yticks(range(5),[SCENARIO_LABELS[s] for s in SCENARIOS] if list(QUANTITIES).index(q)%2==0 else ['']*5,fontsize=9);ax.invert_yaxis();ax.set_title(textwrap.fill(title,28),fontsize=10,loc='left');style_axis(ax)
            if metric=='coverage':ax.axvline(.95,color=MID_GREY,ls='--',lw=1);ax.set_xlim(0,1.03);ax.set_xlabel('95% HDI coverage')
            elif metric=='bias':ax.axvline(0,color=MID_GREY,ls='--',lw=1);ax.set_xlabel('Bias (parameter units)')
            else:ax.set_xlim(left=0);ax.set_xlabel('RMSE (parameter units)')
        return f,'Blue: all completed; pink: numerically adequate. Coverage bars: 95% Wilson intervals; bias bars: ±1.96 Monte Carlo SE; RMSE has no uncertainty bars. Dashed lines: nominal 0.95 coverage or zero bias. SBC is prior predictive. No screen-wide false-discovery guarantee.'
    if kind=='calibration_completion':
        shown=[s for s in SCENARIOS if s in set(df.scenario)];counts=df.set_index('scenario').loc[shown]
        f,ax=panel(4.2,1.85,left=1.75,bottom=.65,top=.4,right=.2);y=np.arange(len(shown))
        ax.barh(y,counts.adequate,color=BLUE);ax.barh(y,counts.failed,left=counts.adequate,color=MAGENTA);ax.set_yticks(y,[SCENARIO_LABELS[s] for s in shown]);ax.invert_yaxis();ax.set_xlabel('Completed calibration replicates');style_axis(ax);set_title(ax,'PTM calibration numerical checks')
        return f,'Completed simulated datasets; blue passes numerical checks, pink fails. Completion and numerical adequacy are distinct.'
    if kind=='sbc_ranks':
        f,axs=plt.subplots(2,3,figsize=(7.6,5.5));f.subplots_adjust(left=.08,right=.98,bottom=.09,top=.94,wspace=.3,hspace=.5)
        for ax,(q,title) in zip(axs.flat,QUANTITIES.items()):
            g=df[df.quantity.eq(q)].sort_values('bin_left')
            if g.empty:ax.set_visible(False);continue
            edges=np.r_[g.bin_left.to_numpy(),g.bin_right.iloc[-1]];centres=(g.bin_left.to_numpy()+g.bin_right.to_numpy())/2
            counts=g['count'].to_numpy();bc=g.failed.to_numpy();expected=g.expected.to_numpy();low=g.reference_low.to_numpy();high=g.reference_high.to_numpy()
            ax.bar(centres,counts,width=.112,color=BLUE,edgecolor='white',zorder=2)
            ax.bar(centres,bc,bottom=counts-bc,width=.112,color=MAGENTA,edgecolor='white',zorder=3)
            ax.plot(centres,expected,color=INK,lw=1,ls='--');ax.step(edges,np.r_[high,high[-1]],where='post',color=MID_GREY,lw=.8);ax.step(edges,np.r_[low,low[-1]],where='post',color=MID_GREY,lw=.8)
            ax.set_title(textwrap.fill(title,24),fontsize=10,loc='left');ax.set_xlabel('Rank / (draws + 1)');ax.set_ylabel('Replicates' if list(QUANTITIES).index(q)%3==0 else '');ax.set_xlim(0,1);ax.set_ylim(0,max(float(counts.max()),float(high.max()))*1.22);ax.text(.98,.94,f'n={counts.sum()}; failed={bc.sum()}',ha='right',va='top',transform=ax.transAxes,fontsize=8);style_axis(ax,grid='y')
        return f,'Prior-predictive SBC only. Blue: adequate; pink: failures retained. Dashed expected counts and gray pointwise 95% binomial limits use discrete uniform ranks. Not simultaneous bands or a formal calibration test.'
    raise ValueError('Unknown manuscript plot type')

def render(request_json):
    global _CAPTURED
    req = json.loads(request_json)
    df = pd.DataFrame(req['rows'])
    if df.empty:
        raise ValueError('No rows match these filters.')
    kind = req['type']
    limit = min(60, max(1, int(req.get('limit',20))))
    title = req.get('title','Extended results')
    note = ''
    if kind in ('ptm_primary','ptm_colourbar','manuscript_heatmap','sensitivity_grid','sensitivity_comparison','sensitivity_classification','calibration_coverage','calibration_bias','calibration_rmse','calibration_completion','sbc_ranks'):
        fig,note=manuscript_plot(df,req)
    elif kind in ('forest','heatmap'):
        for c in ['units','contrast','component','variant','entity_level']:
            if c in df and df[c].fillna('').nunique() > 1:
                raise ValueError('Choose one '+c.replace('_',' ')+' before drawing. Different quantities must not be pooled.')
        for c in ['estimate','hdi_lower','hdi_upper']:
            df[c] = pd.to_numeric(df[c], errors='coerce')
        df = df[np.isfinite(df[['estimate','hdi_lower','hdi_upper']]).all(axis=1)].copy()
        if df.empty:
            raise ValueError('No finite effect intervals in this selection.')
        units = str(df['units'].iloc[0])
        df = df.sort_values('estimate', key=lambda x:x.abs(), ascending=False)
        if kind == 'forest':
            df = df.head(limit).iloc[::-1]
            fig, ax = plt.subplots(figsize=(8,max(3,len(df)*.25+1.3)))
            for i,(_,r) in enumerate(df.iterrows()):
                color = MAGENTA if r['estimate'] >= 0 else BLUE
                ax.plot([r.hdi_lower,r.hdi_upper],[i,i],color=color,lw=1.8)
                ax.scatter([r.estimate],[i],s=26,edgecolor=color,facecolor=color if r.interval_excludes_zero else 'white',zorder=3)
            labels=[str(r.get('label') or r['feature']) + (' · '+str(r['cell_type']) if r.get('cell_type') else '') for r in df.to_dict('records')]
            ax.set_yticks(range(len(df)),labels=labels,fontsize=8)
            ax.axvline(0,color='#6F7378',lw=.8,ls='--');ax.set_xlabel(units)
            ax.set_title(title,loc='left',fontsize=11);ax.grid(axis='x',alpha=.5)
            fig.tight_layout()
            note=f'{len(df)} rows with largest absolute posterior means in the filtered selection. Bars: pointwise 95% HDIs; filled dots: HDI excludes zero. No discovery claim or multiplicity correction.'
        else:
            # Unique feature strings distinguish multiple PTM sites within a gene.
            if 'cell_type' not in df:df['cell_type']=None
            if df['cell_type'].notna().any():
                df['signature_label']=df['label'].astype(str)
                df['celltype_label']=df['cell_type'].fillna('Other')
            else:
                df['signature_label']=df.apply(lambda r: str(r.get('label') or '')+' · '+str(r['feature']) if r.get('label')!=r['feature'] else str(r['feature']),axis=1)
                df['celltype_label']=df['component'].astype(str)
            row_names=list(dict.fromkeys(df.signature_label))[:limit]
            columns=list(dict.fromkeys(df.celltype_label))
            df=df[df.signature_label.isin(row_names)].copy()
            if df.duplicated(['signature_label','celltype_label']).any():
                raise ValueError('More than one estimate per heatmap cell. Narrow the filters.')
            df['hdi_95_excludes_zero']=df.interval_excludes_zero
            vmax=max(float(df.estimate.abs().max()),.01)
            draw_heatmap(df,row_names,columns,{}, {},vmax,'browser','extended','', 'estimate')
            fig=_CAPTURED
            fig.axes[0].set_xlabel('Cell type' if df.cell_type.notna().any() else 'Effect component')
            fig.axes[0].set_ylabel('Signature / feature')
            note=f'{len(row_names)} rows ranked by maximum absolute posterior mean. Scale: −{vmax:.3g} (blue) to +{vmax:.3g} (magenta), white = 0; units: {units}. Circle: pointwise 95% HDI excludes zero. Empty cells are unavailable.'
    elif kind == 'signatures':
        if df.cell_type.nunique()!=1:
            raise ValueError('Select one cell type to compare its signatures.')
        df=df.sort_values('signature').copy()
        df['gene_count'] = df.genes.map(len)
        labels=df.signature.str.replace('_Signature','',regex=False).str.replace('_',' ',regex=False)
        fig,ax=plt.subplots(figsize=(max(10,len(df)*.25),6))
        ax.bar(np.arange(len(df)),df.gene_count,color=TEAL,width=.75)
        ax.set_xticks(np.arange(len(df)),labels=labels,rotation=65,ha='right',fontsize=8)
        ax.set_ylabel('Number of genes');ax.set_xlabel('Signature');ax.set_title(str(df.cell_type.iloc[0]).replace('_',' ')+' — signature composition',loc='left')
        ax.yaxis.get_major_locator().set_params(integer=True);fig.tight_layout()
        note=f'All {len(df)} signatures in the filtered cell type are shown. Gene counts describe set composition, not expression, enrichment or statistical significance. Gene membership is listed in the table.'
    elif kind == 'simulation':
        for field in ['kind','geometry','contrast','rule']:
            if field in df and df[field].fillna('').nunique()>1:
                raise ValueError('Select one '+field+' before comparing scenarios.')
        metric=req.get('y')
        if metric not in df:raise ValueError('Choose an available simulation metric.')
        df[metric]=pd.to_numeric(df[metric],errors='coerce');df=df.dropna(subset=[metric])
        if df.empty:raise ValueError('No finite values for this metric.')
        fig,ax=plt.subplots(figsize=(max(8,df.scenario.nunique()*.7),5))
        if 'replicate' in df:
            sns.boxplot(data=df,x='scenario',y=metric,color=TEAL,ax=ax,fliersize=2)
            note='Distribution across simulation replicates within each scenario. Boxes show the median and interquartile range; whiskers extend to 1.5 times that range.'
        else:
            if df.duplicated('scenario').any():raise ValueError('Narrow the filters to one summary per scenario.')
            ax.scatter(df.scenario,df[metric],color=TEAL)
            mcse=metric.replace('_mean','_mcse') if '_mean' in metric else ''
            if mcse in df:
                err=pd.to_numeric(df[mcse],errors='coerce')
                ax.errorbar(df.scenario,df[metric],yerr=err,fmt='none',color=TEAL,capsize=3)
                note='Scenario-level estimates with ±1 Monte Carlo standard error. These error bars are not 95% confidence intervals.'
            else:note='Scenario-level summaries. No uncertainty bars are inferred from unavailable values.'
        ax.set_xlabel('Scenario');ax.set_ylabel(metric.replace('_',' '));ax.tick_params(axis='x',labelrotation=45);ax.set_title(title,loc='left');fig.tight_layout()
    elif kind == 'intervals':
        for field in ['metric','comparison','layer','scope','stratum']:
            if field in df and df[field].fillna('').nunique()>1:raise ValueError('Select one '+field+' before comparing intervals.')
        for k in ['estimate','ci_95_lower','ci_95_upper']:df[k]=pd.to_numeric(df[k],errors='coerce')
        df=df.dropna(subset=['estimate','ci_95_lower','ci_95_upper'])
        if df.empty:raise ValueError('No complete confidence intervals in this selection.')
        if len(df)>120:raise ValueError('Select a cell type or family to display at most 120 intervals.')
        df=df.sort_values('estimate')
        fig,ax=plt.subplots(figsize=(8,max(4,len(df)*.27)))
        labels=[' · '.join(str(r[k]).replace('_',' ') for k in ['family','celltype_label','CellType','comparison'] if k in r and str(r[k])) for r in df.to_dict('records')]
        for i,(_,r) in enumerate(df.iterrows()):
            ax.plot([r.ci_95_lower,r.ci_95_upper],[i,i],color=BLUE,lw=2);ax.scatter(r.estimate,i,color=TEAL,s=25)
        ax.set_yticks(range(len(df)),labels=labels);ax.set_xlabel(str(df['metric'].iloc[0]).replace('_',' ') if 'metric' in df else 'Estimate');ax.set_title(title,loc='left');fig.tight_layout()
        note='Point estimates and reported 95% confidence intervals. These confidence intervals are distinct from Bayesian highest-density intervals.'
    elif kind == 'purity':
        if df.compartment.nunique()!=1:
            raise ValueError('Select one compartment before drawing.')
        df=df.sort_values('max_cluster_purity')
        fig,ax=plt.subplots(figsize=(8,max(4,len(df)*.25)))
        ax.barh(df.CellType,df.max_cluster_purity.astype(float),color=TEAL)
        ax.set_xlim(0,1);ax.set_xlabel('Fraction assigned to dominant cluster');ax.set_ylabel('Cell type / state')
        ax.set_title(title,loc='left',fontsize=11);fig.tight_layout()
        note='Within-type cluster concentration, not tumour purity or deconvolution accuracy. A value of 1 does not rule out two cell types sharing a cluster. The saved purity table does not identify preprocessing mode.'
    elif kind in ('scree','contributions','annotation'):
        fields = ['compartment','analysis_set'] if kind == 'scree' else ['compartment','analysis_set','PC'] if kind == 'contributions' else ['reference','resolution']
        for field in fields:
            if field in df and df[field].nunique() > 1:
                raise ValueError('Select one '+field+' before drawing.')
        fig,ax=plt.subplots(figsize=(8, max(4, len(df)*.22) if kind != 'scree' else 4))
        if kind == 'scree':
            df=df.assign(pc_number=df.PC.str.extract(r'(\d+)',expand=False).astype(int)).sort_values('pc_number')
            ax.bar(df.PC,df.explained_variance_ratio.astype(float)*100,color=TEAL)
            ax.set_xlabel('Principal component');ax.set_ylabel('Explained variance (%)')
            note='Saved PCA of deconvolved bulk-expression scores, not single-cell PCA. No PCA is refitted. Compare tiered and strict QC separately; non-immune geometry has a saved review flag.'
        elif kind == 'contributions':
            df=df.sort_values('percent_axis_contribution')
            ax.barh(df.CellType,df.percent_axis_contribution.astype(float),color=BLUE)
            ax.set_xlabel('Contribution to PC axis (%)');ax.set_ylabel('Cell type')
            note='Saved squared-axis contributions aggregated by cell type. These are PCA contributions, not abundance fractions or BMI effect estimates.'
        else:
            df=df.sort_values('cell_count')
            ax.barh(df.annotation,df.cell_count.astype(int),color=TEAL)
            ax.set_xlabel('Annotated cells');ax.set_ylabel('Reference label')
            note='Saved reference-based annotation counts. Reference panels are alternative annotations of cells and must not be added together. Counts do not establish annotation accuracy or BMI association.'
        ax.set_title(title,loc='left',fontsize=11);fig.tight_layout()
    elif kind == 'coverage':
        fig,ax=plt.subplots(figsize=(8,4))
        stages=sorted(df.stage_order.astype(int).unique());width=.32
        for k,(assay,color) in enumerate([('phospho',TEAL),('glyco',MAGENTA)]):
            part=df[df.assay==assay].copy();part['stage_order']=part.stage_order.astype(int)
            vals=[int(part.loc[part.stage_order==s,'count'].iloc[0]) if s in part.stage_order.values else 0 for s in stages]
            ax.bar(np.arange(len(stages))+(k-.5)*width,vals,width,color=color,label=assay)
        ax.set_yscale('log');ax.set_ylabel('Features (log scale)')
        names=[str(df.loc[df.stage_order.astype(int)==s,'stage'].iloc[0]).replace('_','\n') for s in stages]
        ax.set_xticks(range(len(stages)),labels=names,fontsize=8);ax.legend();fig.tight_layout()
        note='Counts use the original coverage-figure filters. The numerical-check stage includes primary fits only.'
    else:
        x,y=req.get('x'),req.get('y')
        if x not in df:raise ValueError('Choose an available numeric variable.')
        df[x]=pd.to_numeric(df[x],errors='coerce')
        if kind=='scatter':
            if y not in df:raise ValueError('Choose a Y variable.')
            df[y]=pd.to_numeric(df[y],errors='coerce');df=df.dropna(subset=[x,y])
            fig,ax=plt.subplots(figsize=(7,4.5));sns.scatterplot(data=df,x=x,y=y,color=TEAL,s=26,alpha=.65,ax=ax)
            note=f'{len(df)} finite rows from the filtered table. Descriptive scatterplot; no regression or inferential test.'
        else:
            df=df.dropna(subset=[x]);fig,ax=plt.subplots(figsize=(7,4.5));sns.histplot(data=df,x=x,bins=min(30,max(1,int(np.sqrt(len(df))))),color=TEAL,ax=ax)
            note=f'{len(df)} finite rows from the filtered table. A distribution of saved values, not posterior draws. No pooling-based inference.'
        if df.empty:raise ValueError('No finite numeric values for these axes.')
        ax.set_title(title,loc='left',fontsize=11);fig.tight_layout()
    # Seaborn can leave axis-label clipping enabled on narrow heatmaps.
    for ax in fig.axes:
        ax.xaxis.label.set_clip_on(False)
        ax.yaxis.label.set_clip_on(False)
    svg=io.StringIO();png=io.BytesIO()
    fig.savefig(svg,format='svg',bbox_inches='tight',pad_inches=.12,bbox_extra_artists=[label for ax in fig.axes for label in (ax.xaxis.label,ax.yaxis.label)])
    fig.savefig(png,format='png',dpi=180,bbox_inches='tight',pad_inches=.12,bbox_extra_artists=[label for ax in fig.axes for label in (ax.xaxis.label,ax.yaxis.label)])
    plt.close(fig)
    return json.dumps({'svg':svg.getvalue(),'png':base64.b64encode(png.getvalue()).decode(),'caption':note})
