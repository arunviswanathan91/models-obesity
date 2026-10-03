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
    if kind in ('forest','heatmap'):
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
            note=f'{len(row_names)} rows ranked by maximum absolute posterior mean. Original study heatmap renderer. Scale: −{vmax:.3g} (blue) to +{vmax:.3g} (magenta), white = 0; units: {units}. Circle: pointwise 95% HDI excludes zero. Empty cells are unavailable.'
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
    svg=io.StringIO();png=io.BytesIO()
    fig.savefig(svg,format='svg',bbox_inches='tight',pad_inches=.12)
    fig.savefig(png,format='png',dpi=180,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    return json.dumps({'svg':svg.getvalue(),'png':base64.b64encode(png.getvalue()).decode(),'caption':note})
