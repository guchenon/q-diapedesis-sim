# Q-DIAPEDESIS
**Physiology-Guided Quantum Floquet Control of Biomimetic Microrobotic Swarms for Cerebral Microvascular Recanalization**

[![OSF Preregistration](https://img.shields.io/badge/OSF_Registry-osf.io%2Fhtswk-blue)](https://osf.io/htswk)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0008--1357--9003-brightgreen)](https://orcid.org/0009-0008-1357-9003)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview
`q-diapedesis-sim` is the official 100% in-silico multi-physics and open-system Floquet–Lindblad simulation framework for the **Q-DIAPEDESIS** protocol. It investigates whether untethered magnetic microrobotic swarms can recanalize occluded cerebral penetrating arterioles (15–50 µm) to resolve the post-thrombectomy no-reflow phenomenon while preserving neurovascular physiology.

The framework integrates four computational layers:
1. **In-Vivo Microvascular Rheology (`src/env_hemodynamics/`):** Pulsatile non-Newtonian blood flow governed by the Pries–Secomb Fåhræus–Lindqvist viscosity law and cerebral perfusion pressure ($\text{CPP} = \text{MAP} - \text{ICP}$).
2. **Dual-Kingdom Biomimetic Mechanics (`src/biomimicry_mechanics/`):** Leukocyte-inspired selectin catch-bond boundary margination and decentralized hydro-magnetic phase synchronization.
3. **Floquet–Lindblad Swarm Control (`src/quantum_controller/`):** Bimodal Floquet Hamiltonian engineering ($80\text{ Hz}$ thrombolytic drilling and $5\text{ Hz}$ distal embolic trapping) paired with Lindblad dissipative tensor networks and sub-voxel state estimation.
4. **Neurocritical Control Barrier Functions (`src/neurocritical_cbf/`):** Real-time safety filters enforcing endothelial glycocalyx preservation, shear-induced NETosis prevention ($\tau_{\text{wall}} \le 5.5\text{ Pa}$), pulsed acidotic buffering, and magnetothermal safety ($\Delta T < 0.3\ ^\circ\text{C}$).

Reproducibility & Preregistration
The mathematical formulation, falsifiable hypotheses, and physiological safety thresholds of this study are preregistered on the Open Science Framework (OSF): https://osf.io/htswk.

Author
Chanon Singson — ORCID: 0009-0008-1357-9003

License
This project is licensed under the MIT License — see the LICENSE file for details.

## Repository Structure
```text
q-diapedesis-sim/
├── README.md                          # Project overview and architecture
├── LICENSE                            # MIT License
├── CITATION.cff                       # Bibliographic citation metadata (ORCID & OSF linked)
├── environment.yml                    # Reproducible Python 3.11 environment specification
├── configs/                           # Centralized YAML simulation parameters
│   ├── vessel_mra_config.yaml
│   ├── quantum_floquet_config.yaml
│   └── neurocritical_cbf_config.yaml
├── src/                               # Core multi-physics & control modules
│   ├── env_hemodynamics/
│   ├── biomimicry_mechanics/
│   ├── quantum_controller/
│   ├── neurocritical_cbf/
│   └── observer_retrieval/
├── notebooks/                         # Reproducible Jupyter execution workflows
├── experiments/                       # Benchmark & ablation study logs
└── docs/                              # Preregistered mathematical specifications
