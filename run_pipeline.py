import os

from analysis_pipeline.lc_fitting import fit_csp, fit_des, fit_sdss

file_dir = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.join(file_dir, 'pipeline_outputs')
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

fit_csp(out_dir)
fit_des(out_dir)
fit_sdss(out_dir)
