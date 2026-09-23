# Author-controlled isolated reproduction runtime

The completed reproduction ran on Windows in a new author-controlled working directory on a separate drive. It was isolated from the original result directory, but it was not an independent third-party reproduction.

## Verified Python environment

- Python 3.12.13
- NumPy 2.3.5
- h5py 3.15.1
- Numba 0.63.1
- llvmlite 0.46.0
- threadpoolctl 3.6.0
- openpyxl 3.1.5 (workbook readback only)
- floating-point type: float64

Independent-path continuous-value comparisons used `atol=1e-12` and `rtol=1e-10`. Membership, directions, identifiers, counts, seeds, null extreme counts, decisions and labels required exact equality. A sign disagreement or a significance-boundary disagreement was a hard failure regardless of numerical tolerance.

The R environment recorded for the GSE85358 analysis is separately specified in `environment.yml`; it was not used for the Phase 171B GSE70138 scoring/null computation.
