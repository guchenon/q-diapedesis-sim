# Q-DIAPEDESIS
**Physiology-Guided Quantum Floquet Control of Biomimetic Microrobotic Swarms for Cerebral Microvascular Recanalization**

[![OSF Registration](https://img.shields.io/badge/OSF_Registration-osf.io%2F8mysz-blue)](https://osf.io/8mysz)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0008--1357--9003-brightgreen)](https://orcid.org/0009-0008-1357-9003)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23599019.svg)](https://doi.org/10.5281/zenodo.23599019)

## Overview
`q-diapedesis-sim` is the official 100% in-silico multi-physics and open-system Floquet–Lindblad simulation framework for the **Q-DIAPEDESIS** protocol. It investigates whether untethered magnetic microrobotic swarms can recanalize occluded cerebral penetrating arterioles (15–50 µm) to resolve the post-thrombectomy no-reflow phenomenon while preserving neurovascular physiology.

The framework integrates four computational layers:
1. **In-Vivo Microvascular Rheology (`src/env_hemodynamics/`):** Pulsatile non-Newtonian blood flow governed by the Pries–Secomb Fåhræus–Lindqvist viscosity law and cerebral perfusion pressure ($\text{CPP} = \text{MAP} - \text{ICP}$).
2. **Dual-Kingdom Biomimetic Mechanics (`src/biomimicry_mechanics/`):** Leukocyte-inspired selectin catch-bond boundary margination and decentralized hydro-magnetic phase synchronization.
3. **Floquet–Lindblad Swarm Control (`src/quantum_controller/`):** Bimodal Floquet Hamiltonian engineering ($80\text{ Hz}$ thrombolytic drilling and $5\text{ Hz}$ distal embolic trapping) paired with Lindblad dissipative tensor networks and sub-voxel state estimation.
4. **Neurocritical Control Barrier Functions (`src/neurocritical_cbf/`):** Real-time safety filters enforcing endothelial glycocalyx preservation, shear-induced NETosis prevention ($\tau_{\text{wall}} \le 5.5\text{ Pa}$), pulsed acidotic buffering, and magnetothermal safety ($\Delta T < 0.3\ ^\circ\text{C}$).

Reproducibility & Preregistration
The mathematical formulation, falsifiable hypotheses, and physiological safety thresholds of this study are registered on the Open Science Framework (OSF): registration [osf.io/8mysz](https://osf.io/8mysz), associated project [osf.io/hnstf](https://osf.io/hnstf).

**Verification status (2026-10-24).** The registration GUID resolves on the OSF *registrations* API endpoint, which is what distinguishes a registration from an ordinary project, but it currently answers `401 Unauthorized` — it is private or under embargo and so **cannot yet be checked by a reader**. The parameter lock that *is* independently verifiable is the Zenodo snapshot [10.5281/zenodo.23599020](https://doi.org/10.5281/zenodo.23599020) and the git tag `v0.1.0-alpha` (`b5899d7`), which `experiments/check_continuity.py` re-checks on every run. See ELN #005 and Appendix C of the evidence matrix.

Author
Chanon Singson — ORCID: 0009-0008-1357-9003

License
This project is licensed under the MIT License — see the LICENSE file for details.

## Method naming — what the modules actually implement

Two module names in this repository promise more than the code delivers, and both are recorded
here rather than quietly left in place (D11 / ELN #003):

* **`PINOSubVoxelObserver` (Module B2)** is a **Physics-Informed Spectral MAP Observer**: a
  five-mode regularized Fourier estimator of the cerebral haemodynamic waveform (mean plus the
  first and second cardiac harmonics in sine and cosine), solved as a maximum-a-posteriori fit
  whose weighting follows from the measurement and prior uncertainties. **There is no neural
  network in it and no operator is learned.** `PhysicsInformedSpectralObserver` is provided as
  the accurate alias; the original name is kept so imports do not break. The neural-operator
  literature cited alongside it (Li et al., ICLR 2021) describes a method this code does not use.
  Its measured accuracy advantage comes mostly from the physics prior, not from the estimator:
  giving a Kalman filter the same cardiac-harmonic model closes about three quarters of the gap
  (see the harmonic-aware Kalman baseline it reports), and the remaining margin exists partly
  because the batch fit is a *smoother* that uses the whole window rather than a causal filter.
* **The two-level `LindbladEquivariantTTN` density matrix** is a reduced-order surrogate
  (Status B), not a derivation from the N-body dipolar Hamiltonian, and the magnetothermal term
  is a specific heating rate, not a Pennes bioheat solution.

## Repository Structure
```text
q-diapedesis-sim/
├── README.md                          # Project overview and architecture
├── LICENSE                            # MIT License
├── CITATION.cff                       # Bibliographic citation metadata (ORCID & OSF linked)
├── environment.yml                    # Reproducible Python 3.11 environment specification
├── pytest.ini                         # pythonpath / testpaths (R1-R2, ELN #002)
├── .github/workflows/sanity_ci.yml    # CI: installs pyyaml pytest numpy scipy, runs tests/
├── configs/                           # Centralized YAML simulation parameters
│   ├── vessel_mra_config.yaml
│   ├── quantum_floquet_config.yaml
│   └── neurocritical_cbf_config.yaml
├── src/                               # Core multi-physics & control modules
│   ├── env_hemodynamics/
│   ├── biomimicry_mechanics/
│   ├── quantum_controller/
│   ├── neurocritical_cbf/
│   ├── observer_retrieval/
│   └── physics/                       # First-principles parameter derivations (D7-D8)
├── tests/                             # Preregistration sanity & hypothesis tests (run in CI)
├── experiments/                       # Electronic lab notebook, ablation logs, DOI audit
│   ├── ELN_Logbook.md
│   ├── check_continuity.py            # config -> code -> docs continuity audit
│   ├── check_dois.py
│   └── doi_audit.tsv
└── docs/                              # Preregistered mathematical specifications
```
