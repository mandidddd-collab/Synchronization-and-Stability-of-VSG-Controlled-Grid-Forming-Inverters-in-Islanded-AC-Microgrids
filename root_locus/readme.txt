# Dominant-Pole Root Locus

This folder contains the stand-alone MATLAB implementation for the dominant
eigenvalue root locus in the IEEE 39-bus case study.

## Purpose

The script evaluates the equilibrium equation

```text
B*A*sin(B' * theta) = lambda * varpi
```

along an increasing loading path.  At each equilibrium, it eliminates the
algebraic buses using a Schur complement and forms the reduced VSG state
matrix.  The rotational zero eigenvalue caused by uniform phase shifts is
removed, and the dominant non-rotational pole is tracked.

The plotted color represents `lambda/lambda_4`.  The critical point is the
first phase-angle boundary at which the dominant non-rotational pole reaches
the imaginary axis.  For the supplied IEEE 39-bus data,

```text
lambda_4 = 4.153119151666
lambda_cr = 5.581250408027
lambda_cr/lambda_4 = 1.343869560253
```

The script stops at this critical equilibrium.  It does not invent an
unstable eigenvalue branch beyond the phase-angle boundary.

## Run

Open MATLAB and run:

```matlab
cd('path/to/root_locus')
draw_dominant_root_locus
```

The script is self-contained and does not require MATPOWER, CVX, or SDPT3.
It contains the IEEE 39-bus voltage magnitudes, active-power data, branch
reactances, transformer ratios, VSG parameters, equilibrium continuation,
boundary refinement, and plotting code in one file.

## Outputs

The script writes the following files into this folder:

- `dominant_root_locus.pdf`: vector figure for the paper.
- `dominant_root_locus.png`: raster preview.
- `dominant_root_locus_data.csv`: loading factors, dominant poles, phase
  angles, and residuals.
- `dominant_root_locus_data.mat`: MATLAB data for further analysis.

The main plot includes the full dominant-pole locus and a local view of the
critical real pole near the origin.  The limiting branch and boundary
residual are printed in the MATLAB Command Window.

## Numerical interpretation

The locus describes the linearized reduced-VSG model at synchronous
equilibria.  It is therefore a local small-signal calculation.  The
time-domain nonconvergence at a tested loading level should be interpreted
together with the existence of the synchronous equilibrium and the
phase-angle boundary.
