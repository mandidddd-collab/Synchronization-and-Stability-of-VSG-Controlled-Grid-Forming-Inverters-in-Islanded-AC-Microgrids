# IEEE 39-Bus Simulink Case Study

This folder contains the Simulink model and MATLAB plotting script used for the  frequency-response case study.

## Model

```text
IEEE39_29VSG.slx
```

The model represents the IEEE 39-bus system with 29 VSG buses and applies the sequential loading experiment described in the paper.  
The frequency logging block records the nine source-VSG channels associated with buses 30--38 in the variable `f_vsg30_38_Hz`.

The model uses Simulink and Simscape Electrical Specialized Power Systems.
The supplied files were checked with MATLAB/Simulink R2025a.

## Loading schedule

The default sequential loading schedule is

```text
0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.35 times lambda_4
```

For the bundled model, the loading stages start at approximately

```text
0, 0.05, 3, 6, 9, 12, 15 seconds
```

with

```text
lambda_4 = 4.153119151666
```

The `1.35 lambda_4` stage is the first tested nonconverged case in the
provided simulation setup.

## Run the simulation

1. Open `IEEE39_29VSG.slx` in Simulink.
2. Run the model with the default parameters.
3. Confirm that the model produces the logged variable
   `f_vsg30_38_Hz`.
4. From this folder, run:

```matlab
plot_VSG_frequency_paper_noinit
```

The plotting function does not re-simulate the model and does not modify the
Simulink file.  It reads the logged frequency data and exports the
publication-style figure to `VSG_paper_figures/`.

## Replot the supplied reference data

The repository includes the previously generated frequency data.  If a new
simulation is not needed, the figure can be regenerated from the saved data:

```matlab
S = load(fullfile('VSG_paper_figures', ...
    'VSG_frequency_paper_data.mat'));
data = timeseries(S.f_Hz,S.t);
plot_VSG_frequency_paper_noinit(data);
```

The plotting function also accepts a `SimulationOutput` object or a
time-domain log structure with `time` and `signals` fields.

## Outputs

The `VSG_paper_figures/` directory contains:

- `VSG_frequency_paper.pdf`: vector figure for the paper.
- `VSG_frequency_paper.png`: high-resolution raster figure.
- `VSG_frequency_paper.fig`: editable MATLAB figure.
- `VSG_frequency_paper_data.mat`: original time and frequency data,
  plotting settings, and loading-stage metadata.

The plotting script preserves the original samples.  It does not smooth, rescale, or truncate the frequency response.
