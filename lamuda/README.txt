# VSG synchronization benchmark

This folder contains the reproducible MATLAB code used to compute the four critical loading factors for six benchmark networks.

## Main entry point

Run the following commands in MATLAB:

```matlab
[resultsTable, caseResults, cfg] = run_vsg_two_proposed_conditions_benchmark;
```

The main script calls `run_static_ieee_benchmark.m`. The numerical routines
are in the `+normbench` package, and the six bundled network data files are in
`data_matpower81`.

## Conditions and systems

For each network, the program computes:

- `lambda_1`: existing 2-norm condition;
- `lambda_2`: existing infinity-norm condition;
- `lambda_3`: proposed 2-norm condition;
- `lambda_4`: proposed infinity-norm condition.

The tested systems are WSCC 9-bus, IEEE 14-bus, IEEE 30-bus, New England
39-bus, IEEE 57-bus, and IEEE 118-bus.

The proposed tests implemented in this release are

```text
||u||_2   < g(sqrt((2||P||_2^2 + 1)/3))
||u||_inf < g((2||P||_inf + 1)/3)
```

The benchmark is self-contained. It does not require a separate MATPOWER
installation, CVX, or SDPT3.

## Generated outputs

After a successful run, the following files are written to `results`:

- `static_benchmark_results.csv`: complete numerical results;
- `critical_loading_factors.tex`: LaTeX table for the four loading factors;
- `critical_ratio_comparison.tex`: compatibility copy of the table;
- `certified_point_verification.csv`: certificate checks;
- `static_benchmark_workspace.mat`: MATLAB workspace containing all results;
- `converted_case_data/`: converted network data used by the benchmark.

The script also prints the normalized loading ratios and the independent
phase-cohesive boundary for each network in the MATLAB command window.

## Notes

The code uses the bundled case files and intentionally disables optional AC
power-flow calls by default. This keeps the table reproducible without
external solver packages. To rerun the study from a clean state, delete the
`results` folder and execute the main entry point again.
