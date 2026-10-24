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

---

## Log Entry #002 — CI Repairs (R1–R3) and Code/Hypothesis Audit at `eed0414`
* **Date:** 2026-10-24
* **Stage / Step:** Stages 1–3 (Steps 06–20), audit before Stage 4
* **Prepared by:** Claude (Anthropic), reviewed and approved by Chanon Singson and Argon (decisions D1–D6)

### Resolved issues
| ID | Symptom | Root cause | Fix | Commit |
| :--- | :--- | :--- | :--- | :--- |
| R1 | CI `ModuleNotFoundError: No module named 'src'` after the CI command changed to `pytest tests/ -v` | Repository root not on `sys.path` | Added `pytest.ini` (`pythonpath = .`, `testpaths = tests`) | `9048d82` |
| R2 | pytest exit code 4 | Leading whitespace in `pytest.ini` | Removed indentation | `cb045d3` |
| R3 | `test_stage1_coupled_hemodynamics_and_cbf_protection`: 5.504865 Pa > 5.5001 Pa | One QP relaxation slack δ shared by shear, thermal and CPP constraints; δ (driven by the CPP trade-off) relaxed the shear barrier by δ/κ = 0.0243/5 ≈ 0.0049 Pa | Shear barrier made hard (no slack) plus a post-solve projection; preregistered 5.50 Pa ceiling unchanged; peak shear over the 40-step test 5.44 Pa | `4d52b73` |

### Open findings at `eed0414` (full evidence in the audit report shared with Argon)
* **C1 — Shear ceiling infeasible after recanalization.** Hydrodynamic wall shear alone (no robots) reached 11.22 Pa at 90% and 13.23 Pa at 100% recanalization. The full closed loop failed its own criteria at the configured 10 s horizon (peak shear 13.50 Pa, glycocalyx 0.230, NETosis 0.275). Tests only ran 2–3 s.
* **C2 — Two Stage 3 tests failing.** `step17`: gradient variance 9.5×10⁻⁴ < 1×10⁻³ (and 6.6 < 10 orders of magnitude). `step18`: PINO p95 lost to dead reckoning at 50 ms (1.21 vs 0.69 µm) and 100 ms (1.65 vs 1.47 µm).
* **C3 — Outcomes fixed by model construction.** Stage 1 unfiltered run also passed the final assertions (glycocalyx 0.935, NETosis 0.009, ΔT 0.082 °C); NETosis and mechanical stripping were zero below 5.50 Pa by formula; Kuramoto r(0) = 0.94 > 0.85 target; PINO estimate = ground truth + noise; Kalman given an artificial +3.2 µm; Lindblad state blended toward the target every step (overlap 0.96 vs 0.54 for pure GKLS); TTN gradient variance from a closed-form expression; Floquet fidelity from fixed weights recorded before CBF scaling.
* **C4 — Unfair ablation.** CBF-off baseline multiplied B by 1.95 and rate by 1.8, disabled pH buffering, and reached 100% recanalization in 3 s versus 39% for the full system.
* **H1 — Viscosity law.** Code used the in-vitro η₀.₄₅ term `220·exp(−1.3D)` and `(D/(D+1.1))²`; the Pries in-vivo law uses `6·exp(−0.085D)` and `(D/(D−1.1))²`. Underestimate: −36% at 30 µm, −82% to −92% at 5–10 µm.
* **H2 — CPP rate cap.** The shared slack still relaxed the CPP barrier: at distal pH ≥ 7.10 the filtered rate was 0.142–0.146 s⁻¹ versus the 0.139 s⁻¹ cap (surge 2.62 vs 2.50 mmHg/s).
* **H3 — Unregistered constants** embedded in code (CBF couplings, kinetic rates, step-out frequencies, fidelity weights, retrieval constants, 1200 ms degaussing; 800 ms gives retrieval 0.989 < 0.99).
* **M1/M2, L1/L2** — coupling gaps, borderline washout (1.96% at 50 µm), minor code and documentation items.

### Decisions approved (D1–D6)
D1-c (keep 5.50 Pa; fix H1; add distal bed resistance calibrated to physiological arteriolar shear), D2 a→b (xfail with reason, then real models against unchanged thresholds), D3 b+a (falsifiable mechanistic models, labelled reduced-order Status B), D4 b+a (no multipliers, separate pH switch, equal-time and iso-recanalization comparisons), D5 a+b (constants into `model_implementation_constants`, documented here), D6 (execute work order 1–7, one issue per commit, before/after numbers logged).

### Integrity rules adopted for the work that follows
1. No preregistered threshold is lowered or changed.
2. Model parameters introduced by D1/D3 are fixed from cited sources or stated a-priori rationale **before** outcomes are inspected, and the resulting outcomes are reported as they come out, including negative results.
3. Any change to how a preregistered comparison is evaluated is logged below as a **deviation** and must also be filed on OSF.
