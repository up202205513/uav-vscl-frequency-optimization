# Free Vibration Optimisation of Composite Laminates for a UAV Wing

MATLAB source code for the project of the same name: a genetic algorithm maximises the
fundamental natural frequency (Mode 1) of a two-zone, fibre-steered Variable Stiffness
Composite Laminate (VSCL) panel. Each candidate is evaluated with the CalculiX FE solver,
and an analytical Automated Fibre Placement (AFP) curvature filter enforces manufacturability.

> **Full methodology, equations, results, and discussion are in the project report.**
> This repository contains only the source code and the baseline mesh needed to run it.

## Requirements
- **MATLAB** R2020b or newer.
- **CalculiX `ccx_dynamic`** — required, not included here. The code expects it at
  `tools/PrePoMax v2.5.1 dev/PrePoMax v2.5.1 dev/Solver/ccx_dynamic.exe` (e.g. from a
  PrePoMax install) and prepends that folder to `PATH` at runtime.
- **OS:** Windows (the solver is invoked as `ccx_dynamic.exe`).

## Running
From the MATLAB command window (scripts self-locate, so the current folder doesn't matter):
- **Main multi-start optimisation:** `run('src/optimization/ga_optimization_core.m')`
- **Conventional fixed-angle baseline:** `run('src/postprocessing/conventional_baseline_study.m')`
- **Visualise a steering field:** `addpath('src/preprocessing'); build_vscl_input_deck`

## Layout
- `src/optimization/` — GA driver, single seeded run, fitness evaluator, AFP curvature filter, grid sweep
- `src/preprocessing/` — steered input-deck generator, mesh parser
- `src/postprocessing/` — Mode-1 frequency extractor, conventional baseline sweep
- `fea/mesh/Baseline_Trial.inp` — baseline unsteered CalculiX mesh (S8R shells)

## License
MIT — see [LICENSE](LICENSE).
