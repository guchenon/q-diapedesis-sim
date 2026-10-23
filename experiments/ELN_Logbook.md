# Electronic Lab Notebook (ELN) — Q-DIAPEDESIS
**Project:** Physiology-Guided Quantum Floquet Control of Biomimetic Microrobotic Swarms for Cerebral Microvascular Recanalization  
**OSF Preregistration:** [https://osf.io/htswk](https://osf.io/htswk)  
**Lead Investigator:** Chanon Singson ([ORCID: 0009-0008-1357-9003](https://orcid.org/0009-0008-1357-9003))

---

## Log Entry #001 — Preregistration & Baseline Configuration Lock
* **Date:** 2026-10-23
* **Stage / Step:** Stage 0 (Steps 01–04)
* **Objective:** Lock mathematical protocol on OSF Registries and initialize centralized YAML parameter configurations for 100% in-silico reproducibility.
* **Locked Configuration Files:**
  1. `configs/vessel_mra_config.yaml`: Penetrating arteriole diameter $d \in [15.0, 50.0]\ \mu\text{m}$ (default $30.0\ \mu\text{m}$), baseline $\text{CPP} = 70.0\text{ mmHg}$ ($\text{MAP} = 85.0$, $\text{ICP} = 15.0\text{ mmHg}$), initial post-occlusive $\text{pH} = 6.65$, selectin catch-bond critical force $F^* = 15.0\text{ pN}$, washout ceiling $< 2.0\%$.
  2. `configs/quantum_floquet_config.yaml`: Super-particle representation ($N_{\text{super}} = 200 \equiv 10,000$ particles), bimodal Floquet drive ($\omega_\alpha = 80.0\text{ Hz}$ proximal drill, $\omega_\beta = 5.0\text{ Hz}$ distal mesh), TTN max bond dimension $\chi = 16$, target dual-task fidelity $\ge 0.92$.
  3. `configs/neurocritical_cbf_config.yaml`: Total wall shear ceiling $\tau_{\text{wall}} \le 5.50\text{ Pa}$, endothelial glycocalyx integrity $\ge 0.92$, brain temperature elevation $\Delta T < 0.30\ ^\circ\text{C}$, 4-stage CPP recanalization stepping.
* **Verification Status:** Automated unit test suite (`tests/test_config_sanity.py`) initialized.
