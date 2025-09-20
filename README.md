# Gravitational_Lensing_and_Photon_Rings

This repository contains Python scripts used to solve and analyze gravitational lensing problems in the Ernst Spacetime. The project is part of my master's research on Gravitational Lensing and Photon rings in magnetized black hole environments.

## 📄 Project Report
For a full explanation of the methodology, results, and conclusions, please see the Report (./Lensing and Photon Rings in a Magnetized Black Hole Spacetime.pdf).

## 🗂 Contents

- `Functions.py` – Defines all the functions used in multiple files.
- `Angular_radii(r_O).py` – Use case of the impact parameter determination and calculation of the shadow angular radii for various observer radii.
- `Circular_Timelike_Geodesics.py` – Solves for the case of circular time-like observers in the Equatorial plane and plots the potential and its second derivative for the given magnetic field strength B values.
- `Critical_observer_radius.py` – Solves for the case of  critical observer radius beyond which the horizontal shadow component formulas are invalid for an equatorial observer.
- `Effective_Potential.py` – Solves for and plots the equatorial potential for lightlike geodesics.
- `Gap_paramter.py` – Solves for and plots the gap parameter for list of B-values for a given n (=2 by default).
- `Timelike_Geodesics.py` – Solves equations for the case of circular timelike geodesics and plots the L^2 (z-component of the angular momentum squared), proving that circular timelike geodesics are only possible between r_1 and r_2 (the radii for the circular lightlike geodesics in the equatorial plane).
- `r_roots_calculation.py` – Calculates and plots the roots of the cubic equation in r for a series of B-values.
- `shadow_plots.py` – Plots the horizontal and vertical angular radii of the black hole shadow for a given observer radius r_O.

```bash
pip install numpy matplotlib scipy math

