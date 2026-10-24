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

---

## Log Entry #003 — Implementation of Decisions D1–D6
* **Date:** 2026-10-24
* **Branch:** `claude/argon-d1-d6` (one issue per commit)

| Step | ID | Change | Before | After | Commit message prefix |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 2 | C2 (D2-a) | `step17`, `step18` marked `xfail(strict=True)` with reasons; thresholds unchanged | CI: 10 passed, 2 failed | CI: 10 passed, 2 xfailed | `D2a:` |
| 3 | H2 | CPP barrier made hard (slack removed), rate bound set to the staged/pH-gated cap, post-solve clip; new test `tests/test_cbf_cpp_rate_cap.py` (18 cases) written first | 16/18 cases over cap. Rate/cap at rf=0.5, pH 6.65: **0.1765 / 0.0208 (8.5× over; pH gate ineffective)**; at pH 7.15: 0.1455 / 0.1389 | 18/18 within cap. rf=0.5, pH 6.65: 0.0208 / 0.0208; pH 7.15: 0.1389 / 0.1389; rf=0, pH 6.65: 0.1111 (was 0.1160) | `Fix H2:` |
| 4a | H1 | Viscosity law corrected to Pries 1994 in-vivo form (`6·exp(−0.085D)`, `z = (D/max(D−1.1, 1.5))²`); test `tests/test_pries_in_vivo_viscosity.py` written first (5 failed before) | η at 5/10/30/50 µm: 1.26 / 1.26 / 1.92 / 2.30 mPa·s. Hydro peak τ at rf = 0 / 0.9 / 1.0 (d = 30): 0.89 / 11.22 / 13.23 Pa. Closed loop 10 s full: rf 0.55, PASS; d = 15 µm 10 s: peak 5.86 Pa, FAIL | η: 15.14 / 7.05 / 3.03 / 2.74 mPa·s. Hydro peak τ: 1.07 / 11.53 / 13.23 Pa (open-vessel τ = ΔP·r/2L is viscosity-independent, so C1 persists). Closed loop 10 s full: rf 0.41 (slower pH buffering), PASS; d = 15 µm: 5.50 Pa, PASS | `Fix H1:` |
| 4b | C1 (D1-c) | Added series distal-bed resistance `R_distal_bed` = 2.144×10¹⁴ Pa·s/m³, calibrated **once** at 30 µm to the cited open-vessel peak WSS 3.9 Pa (Rovainen et al., JCBFM 1993: 39 ± 14 dyn/cm² in 21–71 µm surface arterioles; the authors note several-fold higher shear is likely in the smallest terminal arterioles) and then held fixed for all diameters. CPP pressure-drop fraction (0.08) moved to `implementation_constants`. Test `tests/test_distal_bed_calibration.py` written first (1 failed before) | Hydro peak τ at rf = 0.9 / 1.0 (d = 30): 11.53 / 13.23 Pa; rf = 1 at d = 15 / 50: 6.39 / 22.35 Pa. CBF-off 10 s: peak 18.75 Pa, glycocalyx 0.000 | Hydro peak τ: 4.19 / 3.89 Pa; d = 15 / 50: **5.93** / 0.97 Pa. CBF-off 10 s: peak 9.42 Pa, glycocalyx 0.103, NETosis 0.301 (FAIL). Full system 10 s: rf only 0.41 (pH-gated cap), PASS | `D1:` |

**Findings from step 4 (reported, not tuned):**
* At 15 µm the open-vessel hydrodynamic shear alone (5.93 Pa) exceeds the 5.50 Pa ceiling, so the preregistered ceiling cannot hold for the narrowest vessel in the sweep once recanalization completes.
* With the CPP/pH gate now enforced (step 3), the full controller recanalizes only ~41% of the clot in the configured 10 s horizon; the open-vessel regime is therefore not reached within 10 s.
| 5 | C4 (D4), M1 (part) | Removed ×1.95 / ×1.8 multipliers (CBF-off applies the same nominal command); new `enable_pulsed_pH_buffering` switch separate from `enable_cbf_filter`; `stop_at_recanalized_fraction` for iso-recanalization checkpoints; default horizon from `simulation_meta.total_time_seconds` (10 s); peak CPP surge added to `all_hypotheses_passed` (preregistered 2.50 mmHg/s); `mode_label` names every disabled module; `step20` now runs 10 s. Tests `tests/test_ablation_fairness.py` written first (4 failed before) | CBF-off (3 s): field ×1.95, rate ×1.8, pH buffering also off; rf 1.00 vs 0.39 for the full system | See table below | `D4:` |

**Step 5 ablation results, d = 30 µm (preregistered thresholds; PASS = all gated hypotheses):**

| Configuration | Equal time 10 s: rf / peak τ / glycocalyx / NETosis / peak CPP surge / PASS | Iso-recanalization (rf = 0.41): time / peak τ / glycocalyx / peak CPP surge / PASS |
| :--- | :--- | :--- |
| Full | 0.41 / 3.58 Pa / 0.993 / 0.000 / 2.02 / PASS | 10.0 s / 3.58 / 0.993 / 2.02 / PASS |
| No CBF | 1.00 / 6.73 Pa / 0.846 / 0.013 / 5.76 / FAIL | 1.3 s / 3.57 / 0.999 / 5.76 / FAIL (CPP only) |
| No bimodal Floquet | as full; fidelity 0.523 / FAIL | as full; FAIL (fidelity) |
| No catch-bond | as full; washout 4.24%, retrieval 0.9884 / FAIL | same / FAIL |
| No degaussing | as full; retrieval 0.8679 / FAIL | same / FAIL |
| No pH buffering | 0.39 / 3.58 / 0.991 / 0.000 / 2.00 / PASS (final pH 6.67 vs 6.78) | 11.05 s / 3.58 / 0.990 / 2.00 / PASS |

**Findings from step 5:** At equal progress the CBF's measurable benefit is the CPP surge limit; its shear/glycocalyx protection only appears once recanalization goes further than the full controller reaches in 10 s. Pulsed pH buffering changes the 10 s outcome very little (pH 6.78 vs 6.67) because buffering is slow relative to the horizon. NETosis and final pH have no preregistered end-of-episode threshold, so they are reported but not gated.
