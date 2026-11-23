# gauge-arbitrage-theory
Vázquez and Farinelli (V&F, arXiv 0908.3043) propose a gauge-invariant measure of arbitrage ($`\alpha^A`$) and a detector for it. We want to find out what that detector responds to. Specifically, we want to know what it means mathematically (what statistic of the return process it computes), what it means financially (what market behavior makes it positive, negative, or zero), and whether it is equivalent to, or a disguised version of, something already known, such as return autocorrelation, lead-lag structure, or cointegration. We are currently in the simulation phase, and we are finding out which kinds of market behavior make V&F's arbitrage detector read positive, negative, or zero using Monte Carlo simulations.

## Objectives
- Implement the detector from V&F and reproduce the paper's 21 asset simulation
- Run the detector on simulated market data (autocorrelation, lead-lag, cointegration)
- Derive the detector's outputs for each model and check them using Monte Carlo

## Layout
```
gauge_arbitrage_theory/
├── README.md           
├── walkthrough.ipynb                 Project walkthrough
├── docs/
│   ├── figures/                      Output plots
│   ├── DEFINITIONS.md                Shared definitions
│   ├── STATUS.md                     Progress tracker
│   ├── project_overview.md           Reference doc
│   ├── simulation_phase_overview.md  Reference for simulation phase
├── notebooks/                
├── results/                  Simulation phase results
├── src/gauge_arbitrage
│   ├── time_series.py        Creates time series using NumPy
│   └── vf/
│   │  ├── detector.py        V&F's algorithm (section 4) used in the simulations
│   │  ├── predict.py         Prediction calculator
│   │  └── vf_series.py       V&F's 21 asset simulation random walks (Fig. 1 in V&F)
```

## Results
