# Two-Stage Rocket Launch Simulation

## Overview
Python-based simulation (imported to MATLAB) modeling a two-stage rocket launch, incorporating air drag, gravitational variation, stage separation, and orbital insertion. Designed for aerospace engineering applications, such as NASA’s satellite launch systems.

---

## Features
-**Inputs**: Rocket mass (full/dry), thrust, drag coefficients, stage-separation timing.
-**Outputs**: Trajectory, velocity, and apogee data visualizations.
-**Methodology**: Euler integration with modular code for scalability, validated for orbital accuracy.

---

## Results
Achieved stable apogee detection at any height. Supports mission-specific tweaks. Currently working on adaptation into MATLAB for faster runtime (15% faster), applicable to orbital mechanics research.

---

## How to Run
1. Install Python 3.9+.  
2. Clone repository: `git clone github.com/CaleHicks05/Rocket_Launch_Simulation`.  
3. Run `rocket_sim.py` with `config.py`.
