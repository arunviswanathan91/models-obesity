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

    def test_coverage_at_one_preserves_estimate_and_interval(self):
        # Actual template summaries contain 1.0 with a Wilson upper bound one
        # floating-point step below 1.0. This previously raised negative xerr.
        rows=[dict(scope=scope,scenario="null",quantity=q,coverage=1.0,
                   coverage_low=.8668077490609515,
                   coverage_high=p.np.nextafter(1.0,0.0))
              for q in p.QUANTITIES for scope in ["All completed","Numerically adequate"]]
        fig,_=p.manuscript_plot(p.pd.DataFrame(rows),dict(type="calibration_coverage"))
        try:
            fig.canvas.draw()
            for ax in fig.axes:
                for bars in ax.containers:
                    self.assertEqual(bars.lines[0].get_xdata()[0],1.0)
                    segment=bars.lines[2][0].get_segments()[0]
                    self.assertEqual(segment[0,0],.8668077490609515)
                    self.assertEqual(segment[1,0],1.0)
        finally:p.plt.close(fig)

    def test_roundoff_only_is_clamped(self):
        errors=p.np.array([-.5*p.np.finfo(float).eps,0.,.13,p.np.nan])
        actual=p.nonnegative_error_lengths(errors,'Test')
        p.np.testing.assert_equal(actual,[0.,0.,.13,p.np.nan])
        self.assertLess(errors[0],0.,'source values must not be mutated')
        for invalid in [[-.001],[p.np.inf],[-p.np.inf]]:
            with self.assertRaisesRegex(ValueError,'source interval or standard error'):
                p.nonnegative_error_lengths(invalid,'Test')

    def test_zero_coverage_lower_bound_roundoff(self):
        rows=[dict(scope='All completed',scenario='null',quantity='beta_protein',
                   coverage=0.,coverage_low=p.np.finfo(float).eps,coverage_high=.14)]
        result=json.loads(p.render(json.dumps(dict(type='calibration_coverage',rows=rows))))
        self.assertIn('<svg',result['svg'])

    def test_simulation_standard_errors_share_roundoff_guard(self):
        rows=[dict(scenario='null',coverage_mean=.95,coverage_mcse=-1e-16),
              dict(scenario='signal',coverage_mean=.96,coverage_mcse=.02)]
        result=json.loads(p.render(json.dumps(dict(type='simulation',y='coverage_mean',rows=rows))))
        self.assertIn('Monte Carlo',result['caption'])
        rows[0]['coverage_mcse']=-.02
        try:
            with self.assertRaisesRegex(ValueError,'source interval or standard error'):
                p.render(json.dumps(dict(type='simulation',y='coverage_mean',rows=rows)))
        finally:p.plt.close('all')

if __name__=="__main__":unittest.main()

