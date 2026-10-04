import importlib.util, json, unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location("plotting",Path(__file__).resolve().parents[1]/"web/plotting.py")
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

class ManuscriptRenderer(unittest.TestCase):
    def test_primary_heatmap_preserves_cells_and_flags(self):
        rows=[dict(short_label="GENE | S1",gene="GENE",feature_id="S1",posterior_mean_marginal=.2,posterior_mean_conditional=-.3,hdi_95_lower_marginal=-.1,hdi_95_upper_marginal=.4,hdi_95_lower_conditional=-.5,hdi_95_upper_conditional=-.1)]
        fig,note=p.manuscript_plot(p.pd.DataFrame(rows),dict(type="ptm_primary",rows=rows))
        ax=fig.axes[0]
        self.assertEqual(len(ax.patches),1)
        self.assertEqual([x.get_text() for x in ax.get_xticklabels()],["Marginal PTM","Protein-conditioned PTM"])
        self.assertEqual(ax.patches[0].center,(1.5,.5))
        self.assertAlmostEqual(ax.get_position().width*fig.get_figwidth(),2*p.CELL_INCH)
        p.plt.close(fig)
    def test_subset_adjustment_keeps_reference(self):
        row=dict(key="k",gene_base="GENE",feature_id_base="S1",comparison="Stage adjustment (matched)",variant_base="stage_subset_base",variant_alt="stage_adjusted")
        for suffix,mean in [("base",.1),("alt",.2)]:
            row.update({"posterior_mean_"+suffix:mean,"hdi_95_lower_"+suffix:-.1,"hdi_95_upper_"+suffix:.4,"numerical_checks_pass_"+suffix:True})
        fig,_=p.manuscript_plot(p.pd.DataFrame([row]),dict(type="sensitivity_comparison"))
        self.assertIn("Stage subset versus Stage adjusted",fig.axes[0].get_title(loc="left"))
        p.plt.close(fig)
    def test_literal_null_calibration_scenario(self):
        rows=[dict(scope=scope,scenario="null",quantity=q,n=50,coverage=.96,coverage_low=.85,coverage_high=.99,bias=.01,bias_mcse=.01,rmse=.1) for q in p.QUANTITIES for scope in ["All completed","Numerically adequate"]]
        result=json.loads(p.render(json.dumps(dict(type="calibration_coverage",rows=rows))))
        self.assertIn("Null",result["svg"])
        self.assertIn("Wilson",result["caption"])

if __name__=="__main__":unittest.main()
