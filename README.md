# FlashGrinder3D / FlashGrinderX

**Thermo-Mechanical and Kinematic Simulation Suite for Laser-Textured cBN Grinding Wheels**

---

## Overview

**FlashGrinderX** is a MATLAB-based simulation and modeling suite developed for analyzing, simulating, and optimizing precision grinding processes using laser-textured cubic Boron Nitride (cBN) grinding wheels. 


---

## Core Capabilities

- **Kinematic Grit Trajectory Simulation:** Computes 3D discrete trajectories for individual abrasive grits considering wheel speed, feed rate, machine-tool axis offset (eccentricity/runout), and initial wheel orientation angle (`GTracker4.m`).
- **3D Surface Topography & Stochastic Fitting:** Processes confocal microscopy surface data to statistically model grit height distributions, rake angles, clearance angles, nose radii, and widths of cut (`GSpecifier.m`, `GParameters.m`).
- **In-Situ Grit-Workpiece Engagement:** Dynamically updates workpiece surface profile and calculates instantaneous uncut chip thickness distributions (`WIntersections.m`, `WSurface2.m`, `intersections.m`).
- **Thermo-Mechanical Force Prediction:** Calculates cutting forces, flank-face ploughing forces (plastic sticking and elastic sliding zones), and dead metal zone (DMZ) forces using oblique cutting mechanics and Johnson-Cook material parameters for AISI 4140 steel (`GForce1.m`, `GForce2.m`, `OneGrit.m`).
- **Laser Texture Geometry Generator:** Simulates various laser-patterned wheel topographies (stripe width, orientation angle, grooved surface area ratio) and evaluates passive grit activation (`GPatterner.m`, `GShaper.m`).
- **Thermal & Surface Roughness Modeling:** Computes contact zone finite-difference thermal fields and predicts ground surface roughness ($R_a$) in feed and cross-feed directions (`GTempFD.m`, `GTempFD3D3.m`, `SRoughness.m`).

---

## Repository Structure

```
FlashGrinderX/
├── script.m                      # Main simulation execution pipeline
├── GTracker4.m                   # Kinematic grit trajectory generator
├── GForce1.m & GForce2.m         # Multi-zone thermo-mechanical force solvers
├── GPatterner.m & GShaper.m      # Texture pattern & slot geometry generators
├── GParameters.m & GSpecifier.m  # Grit statistics & geometrical parameter solvers
├── GTempFD.m & GTempFD3D3.m      # 2D/3D Finite Difference thermal solvers
├── OneGrit.m                     # Single-grit interaction mechanics
├── SRoughness.m                  # Surface roughness (Ra) estimation module
├── WIntersections.m & WSurface2.m# Dynamic workpiece surface height map update
├── Deflector.m & ArcPlotter.m    # Dynamic deflection & trajectory visualization
├── intersections.m               # Kinematic intersection detection routine
├── mtrlz.txt                     # Workpiece material constitutive properties (AISI 4140)
├── grtz.txt                      # Grit topography statistical distribution input
├── clntz.txt                     # Coolant / boundary condition settings
├── vars2.txt                     # Simulation process state parameters
└── vars-Just-For-Explanation.txt # Variable naming key and documentation
```

---

## Requirements & Quick Start

### Prerequisites
- **MATLAB** R2020b or later (Image Processing & Signal Processing Toolboxes recommended).

### Running the Simulation
1. Clone the repository:
   ```bash
   git clone https://github.com/v1v4v/FlashGrinder3D.git
   ```
2. Open MATLAB, set the current working directory to `FlashGrinderX`.
3. Execute the main script:
   ```matlab
   script
   ```

---

## Citation & Contact

If you use this codebase or simulation methodology in your research, please cite:

```bibtex
@article{moussavi2025texturing,
  title={On texturing cBN grinding wheels and its effect on reduction of ploughing share and increased productivity},
  author={Moussavi, Vahid and Behrouzbaraghi, Suzan and Araghizad, Arash Ebrahimi and Budak, Erhan},
  journal={Journal of Manufacturing Processes},
  volume={154},
  pages={386--401},
  year={2025},
  publisher={Elsevier}
}
```

*For inquiries, contact Vahid Moussavi at Sabanci University Manufacturing Research Laboratory.*

