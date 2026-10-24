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

| 6a | C3 (D3-b), Stage 1 | Glycocalyx shedding and Piezo1 NETosis made continuous in wall shear: the step at the 5.50 Pa ceiling is replaced by a Hill dose-response `h(τ) = τⁿ/(τ₅₀ⁿ+τⁿ)` in both `glycocalyx_model.py` and `glycocalyx_netosis_dynamics.py`, via the new `src/neurocritical_cbf/shear_dose_response.py`. Oscillatory wear now normalized by τ₅₀ instead of the ceiling; NETosis keeps its −0.01·N clearance. Constants in `configs/neurocritical_cbf_config.yaml → implementation_constants.shear_dose_response`. Test `tests/test_shear_dose_response.py` written first (3 failed before) | Damage exactly **zero** for every τ < 5.50 Pa by construction; discontinuous jump at the ceiling. CBF-off 10 s: glycocalyx 0.846, NETosis 0.013 | Monotone and continuous across the ceiling (5.49 vs 5.51 Pa differ by <2%). CBF-on 10 s: glycocalyx 0.987, NETosis 0.002; CBF-off 10 s: glycocalyx **0.654**, NETosis **0.056**. Full suite 43 passed, 2 xfailed | `D3 Stage 1:` |

**A-priori parameter choices for step 6a (fixed and written down before any outcome was computed):**
* `τ₅₀ = 8.0 Pa` — the highest shear dose tested by Baratchi et al., *Nat. Commun.* 2024 (80 dyn/cm²), assumed half-maximal. Not fitted.
* `n = 6` — chosen by the stated criterion that the Hill fraction at the cited physiological arteriolar shear (3.9 Pa, Rovainen 1993) stays ≤ 1.5 %; (3.9/8)⁶ ≈ 0.0134.
* `k_max` for each process set so the rate at 8 Pa equals the pre-D3 rate at 8 Pa (glycocalyx 0.045·2.5¹·⁵ = 0.178 /s → k_max 0.356 /s; NETosis 0.12·(2.5/5.5)² = 0.0248 /s → k_max 0.0496 /s). This anchors the new law to the old one at one point rather than re-tuning it to an outcome.

**Sensitivity sweep (step 6a), 10 s closed loop at d = 30 µm, k_max rescaled in each cell to hold the 8 Pa rate fixed:**

| τ₅₀ \ n | 2 (glyco / NET, CBF-on → CBF-off) | 4 | 6 |
| :--- | :--- | :--- | :--- |
| 6.0 Pa | 0.587 / 0.076 → 0.302 / 0.149 | 0.880 / 0.021 → 0.398 / 0.118 | 0.967 / 0.006 → 0.451 / 0.102 |
| 8.0 Pa (a priori) | 0.658 / 0.061 → 0.339 / 0.136 | 0.936 / 0.011 → 0.512 / 0.087 | **0.987 / 0.002 → 0.654 / 0.056** |
| 10.0 Pa | 0.699 / 0.052 → 0.365 / 0.128 | 0.954 / 0.008 → 0.576 / 0.073 | 0.991 / 0.001 → 0.735 / 0.041 |

**Findings from step 6a (reported, not tuned):**
* The *direction* of the result is robust: in all nine cells the CBF-filtered run retains more glycocalyx and shows less NETosis than the unfiltered run. The ordering never reverses.
* The *absolute* targets are not robust. At the a-priori parameters the CBF-on run meets the stated goals (glycocalyx 0.987 ≥ 0.92, NETosis 0.002 < 0.03), but at `n = 2` the same CBF-on run gives glycocalyx 0.587–0.699 and NETosis 0.052–0.076, failing both. Any claim phrased as an absolute retention or NETosis number is therefore a claim about `n`, not about the controller. The Hill exponent is unconstrained by the cited data; nothing in Baratchi et al. 2024 fixes it.
* Removing the step makes the no-CBF baseline fail harder than before (glycocalyx 0.846 → 0.654, NETosis 0.013 → 0.056) because sub-ceiling shear now causes damage. This is a consequence of a more honest damage law, not evidence of stronger protection.
* Glycocalyx retention and NETosis are **not** preregistered gated endpoints at end-of-episode, so none of the above changes a hypothesis outcome. They are descriptive.

| 6b | C3 (D3-b), Stage 2 | `soliton_kuramoto_swarm.py`: phases now initialized over the full circle `uniform(−π, π)` instead of `uniform(−0.6, 0.6)`; per-rotor natural frequencies `ω_i ~ N(ω_0g, σ_ω)` with the detuning `ω_i − ω_0g` entering the update; rotational noise `√(2·D_r·dt)·ξ` added. `σ_ω/ω_0 = 0.02` set a priori in `implementation_constants.kuramoto_dispersion`; `D_r` **derived** in code from Stokes–Einstein–Debye `D_r = kT/(8πηa³)` using the existing config radius 0.75 µm, plasma viscosity 1.20 mPa·s and 310 K → **0.336 rad²/s**. Test `tests/test_kuramoto_frequency_dispersion.py` written first (4 failed before) | Initial `r = 0.95` by construction; no detuning and no noise, so the model could not desynchronize. Stage 2 after 2 s: `r_global` 0.98, β coverage 0.976, washout 1.03%. Closed loop 10 s (d = 30): washout 0.72%, PASS | Initial `r = 0.022` (incoherent, ≈1/√200). Stage 2 after 2 s: `r_global` **0.304**, `r_alpha` **0.131**, β coverage **0.793**, washout **2.78%**. Closed loop 10 s (d = 30): washout **2.50%**, FAIL on that component only. Suite 46 passed, 4 xfailed | `D3 Stage 2:` |

**A-priori parameter choices for step 6b (fixed before any outcome was computed):**
* `σ_ω/ω_0 = 2 %` — below the step-out limit a magnetically driven rotor is frequency-locked to the field, so residual natural-frequency dispersion reflects hydrodynamic heterogeneity only. 2 % is deliberately conservative: real magnetic microparticle populations show 5–15 % diameter CV, and magnetic moment scales with volume. Not fitted.
* `D_r` not chosen at all — derived from Stokes–Einstein–Debye with parameters already in the preregistered config. Result 0.336 rad²/s.

**Critical-coupling check (analytic, Sakaguchi–Kuramoto with Gaussian frequencies and noise, `K_c ≈ σ_ω·√(8/π) + 2·D_r`):**

| Sub-swarm | Drive | σ_ω at 2 % CV | K_c | Config `kuramoto_coupling_K` | Synchronizes? |
| :--- | :--- | :--- | :--- | :--- | :--- |
| α (drill) | 80 Hz = 502.7 rad/s | 10.92 rad/s | **18.1 rad/s** | 4.50 rad/s | **no** (K is 4× below K_c) |
| β (mesh) | 5 Hz = 31.4 rad/s | 0.65 rad/s | 1.71 rad/s | 4.50 rad/s | yes |

**σ_ω sweep (mean of seeds 1–3, 2 s, `K` held at the preregistered 4.50 rad/s):**

| frequency CV | σ_ω (α) | r_global | r_alpha | β coverage | washout % | meets r ≥ 0.85 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0 (identical rotors) | 0.00 | 0.795 | 0.826 | 0.872 | 1.41 | **no** |
| 0.1 % | 0.50 | 0.814 | 0.850 | 0.877 | 1.35 | **no** |
| 0.25 % | 1.26 | 0.744 | 0.729 | 0.879 | 1.55 | no |
| 0.5 % | 2.51 | 0.431 | 0.211 | 0.877 | 2.42 | no |
| 1 % | 5.03 | 0.392 | 0.134 | 0.886 | 2.53 | no |
| **2 % (a priori)** | 10.05 | **0.369** | 0.099 | 0.883 | 2.60 | no |
| 5 % | 25.13 | 0.271 | 0.043 | 0.814 | 2.87 | no |

**Findings from step 6b — THREE PREREGISTERED THRESHOLDS NOW FAIL. Reported, not tuned:**
1. **`kuramoto_order_r ≥ 0.85` fails** (0.304 at the a-priori parameters). The previous pass was an artifact of the initial condition: phases started inside a ±0.6 rad window, so `r` began at ≈0.95 and the assertion `r_final > r_initial` was comparing two already-synchronized states. With no detuning and no noise a Kuramoto model is incapable of desynchronizing, so the threshold was true by construction and tested nothing about the coupling.
2. **The failure is not caused by my choice of σ_ω.** At *zero* frequency dispersion — the most favorable honest case, rotational noise only — `r_global` reaches just 0.795 in 2 s, still short of 0.85. No value of σ_ω in the sweep reaches the threshold. The binding constraint is that the α sub-swarm's critical coupling scales with its 80 Hz drive while `kuramoto_coupling_K = 4.50` does not.
3. **`beta_mesh_coverage_ratio ≥ 0.90` fails** (0.79–0.89 across the whole sweep). Coverage is defined as `0.55 + 0.43·r_beta`, so it can only reach 0.90 if `r_beta ≥ 0.814`; it is not an independent measurement.
4. **`transient_washout_percent < 2.0` fails in the closed loop** (2.50 % at d = 30 µm, 2.49 % at 15 µm, 2.83 % at 50 µm). Washout is also an affine function of `r`, so this is the same finding a third time rather than three independent failures. Note the direction: the **no-CBF** run has *lower* washout (0.74 %) purely because it recanalizes faster, which is not a safety benefit.
5. `step20`'s `all_hypotheses_passed` fails on washout **alone**; peak shear 3.58 Pa, glycocalyx 0.987, fidelity 0.993, retrieval 0.9916 and CPP surge 2.02 all still pass.

**Action taken:** no threshold and no preregistered config value was changed. `test_step13_and_step14_soliton_synchronization_and_washout` and `test_step20_closed_loop_pipeline_and_diameter_sweep` are marked `xfail(strict=True)` with the reasons above, following the precedent approved for D2-a. **This needs a decision from the preregistration owner**, and the options are not equivalent:
* Re-derive `kuramoto_coupling_K` from the magnetic dipole–dipole interaction at the modeled inter-particle spacing, so that it is a physical quantity rather than a free parameter. If the derived value exceeds ~18 rad/s the thresholds may hold; if it does not, the Stage 2 hypothesis is falsified as stated. **This is the only option that resolves the question rather than restating it.**
* Keep `K = 4.50` and record Stage 2 as a falsified hypothesis.
* Changing 0.85, 0.90 or 2.0 % to match the new numbers would be fitting the threshold to the result and is **not** recommended; if chosen it must be filed on OSF as a deviation.

| 6c | C3 (D3-b), Stage 3 GKLS | `lindblad_ttn.py`: the step is now a plain RK4 integration of `dρ/dt = −i[H,ρ] + γ(LρL† − ½{L†L,ρ})` and nothing else. Removed (i) the per-step blend `ρ ← (1−a)ρ + a·ρ_attractor` with `a` floored at 0.01, which produced the attractor overlap directly, and (ii) the factor `0.05` multiplying the Hamiltonian commutator, which scaled `H_eff` almost out of the dynamics. Also removed the `clip` on purity and overlap so integration error cannot be hidden. Docstring now states the two-level collective description is a reduced-order surrogate, **Status B**, not derived from the N-body dipolar Hamiltonian. Test `tests/test_lindblad_gkls_integration.py` written first (2 failed before) | Attractor overlap after 25 steps: **0.85+**, produced by the blend. Purity and overlap clipped to their physical ranges. Trace restored by renormalizing each step | Overlap after 25 steps (1.25 s): **0.550**; at steady state (100 s): **0.622**. Trace conserved to < 1e-9 and Hermiticity to < 1e-12 **without any projection**; eigenvalues ≥ −1e-9; purity conserved to < 1e-6 under γ = 0 (unitary limit); convergence verified fourth-order in dt. Suite 52 passed, 4 xfailed | `D3 Stage 3 GKLS:` |

**Findings from step 6c — a FOURTH preregistered threshold fails. Reported, not tuned:**
* **`dissipative_attractor_overlap > 0.80` fails** (0.550 at the tested 25 steps, 0.622 at steady state). The previous pass came entirely from the blend term, which moved the state onto `|+⟩⟨+|` by at least 1 % per step regardless of the physics — including when the dissipation rate was exactly zero, because the blend coefficient was floored at 0.01.
* **Root cause, quantified:** the preregistered dissipation is far too weak relative to the preregistered drive. `γ_eff = 0.19 /s` against `‖H_eff‖ = 11.4 rad/s` — a ratio of 0.017, so coherent Floquet precession dominates the dissipator by about 60×, and the system never settles near the dissipative fixed point.

| `viscous_dissipation_gamma` | γ_eff / ‖H_eff‖ | overlap at 1.25 s | overlap at steady state |
| :--- | :--- | :--- | :--- |
| **0.15 (preregistered)** | 0.017 | **0.550** | **0.622** |
| 1.0 | 0.111 | 0.597 | 0.623 |
| 3.0 | 0.334 | 0.628 | 0.631 |
| 9.0 | 1.002 | 0.741 | 0.690 |
| 15.0 | 1.670 | 0.809 | 0.766 |
| 30.0 | 3.340 | 0.905 | 0.890 |

* Reaching the preregistered 0.80 **at steady state** requires `viscous_dissipation_gamma ≈ 30`, a **200×** increase over the preregistered 0.15. The 25-step value crosses 0.80 at about 15 (100×), but that is a transient, not the attractor.
* The same structural question arises as in step 6b and it is the same shape of question: a preregistered coupling constant was set at a value that does not scale with the drive it is supposed to compete with. **Decision needed from the preregistration owner.** `viscous_dissipation_gamma` should be derived from the rotational drag of the modeled particle in the modeled viscosity rather than chosen — note that the Stokes–Einstein–Debye drag already derived for step 6b gives a rotational relaxation rate of the right order of magnitude to settle this, so this is answerable from first principles rather than by fitting.
* `step17` was already `xfail(strict=True)` for the C2 gradient-variance failure; its reason now records both failures with these numbers.

| 6d | C2 / C3 (D2-b), barren plateau | New `src/quantum_controller/ttn_gradient_variance.py`. The equivariant branch builds a real `ceil(log2 N)`-layer ansatz of `exp(−iθ J_z/j)·exp(−iφ J_x/j)` on the symmetric (Dicke) subspace, dim `N+1`, with the **local** cost `C = ⟨J_z⟩/j` (exactly the single-site Pauli-Z expectation there), and differentiates it **exactly** via `∂_kC = 2 Re⟨∂_kψ|O|ψ⟩`, cross-checked against central finite differences to < 1e-6. The baseline builds a depth-`n` unstructured hardware-efficient circuit on the **full 2ⁿ state vector** with the **global** cost `C = \|⟨0…0\|ψ⟩\|²` and differentiates it by the parameter-shift rule. The 1e-3 floor moved from hard-code to `implementation_constants.barren_plateau` unchanged (D5). Test `tests/test_ttn_gradient_variance.py` written first (1 failed before, after fixing my own adjoint off-by-one) | `Var = Var[sin(a)·cos(b)/√N]` over 64 random angle pairs — not the gradient of any cost, no circuit built. Baseline was the literal constant `0.25·2^−min(N,30)` — never evaluated at all. So the claimed "orders of magnitude" was a ratio of two hand-written numbers | **Measured.** Equivariant at N = 200: **1.654e−6**. Global-cost baseline measured at n = 4…12, decaying as **2^(−1.94 n)** (fit residual ≤ 0.82 in log₂). Extrapolated to N = 200: 7.6e−118, **labelled as an extrapolation** in the returned report. Improvement **111.3** orders | `D2b:` |

**Measured scaling (this is the real content of the barren-plateau claim):**

| N | equivariant local-cost Var[∂C/∂θ] | | n | unstructured global-cost Var[∂C/∂θ] |
| :--- | :--- | :--- | :--- | :--- |
| 8 | 1.449e−2 | | 4 | 1.993e−3 |
| 16 | 3.730e−3 | | 6 | 1.787e−4 |
| 32 | 6.656e−4 | | 8 | 2.130e−5 |
| 64 | 8.032e−5 | | 10 | 8.030e−7 |
| 128 | 7.493e−6 | | 12 | 4.259e−8 |
| 200 | 1.654e−6 | | 200 | 7.6e−118 *(extrapolated)* |

**Findings from step 6d:**
* **The qualitative claim survives, and is now actually measured.** The equivariant local-cost variance decays as **N^(−2.74)** — polynomially, a finite exponent (log–log fit, max residual 0.52 in log₂). The global-cost baseline decays as **2^(−1.94 n)**, i.e. exponentially. That contrast — Ω(1/poly(N)) versus O(2^−N) — is exactly what Barthel & Miao 2025 and Schatzki et al. 2024 predict, and it is now obtained from gradients that were computed rather than written down. The `variance_improvement_orders_of_magnitude ≥ 10` assertion passes with 111 orders.
* **The quantitative threshold fails, and it is the wrong kind of threshold.** `gradient_variance ≥ 1e-3` is false at N = 200 (1.654e−6). But freedom from barren plateaus does not mean the variance exceeds a fixed constant — it means the decay is polynomial rather than exponential. A fixed floor is automatically violated by *any* polynomially decaying quantity at large enough N, including a perfectly plateau-free one; here it is already violated by N ≈ 24. **The preregistered 1e-3 floor tests the wrong property.** The right preregistered quantity is the decay exponent, which `equivariant_scaling_exponent()` now reports. This is a recommendation to the preregistration owner, not a change I have made: the floor is untouched at 1e-3 and `step17` stays `xfail`.
* **Two problems I found in my own first implementation, recorded because they changed the numbers.** (i) The adjoint recursion was off by one unitary, which the finite-difference cross-check caught — the first version disagreed with finite differences by 1.15. This is why the cross-check is a permanent test and not a one-off. (ii) The baseline's rotation axes were assigned deterministically as `i mod 3`, which aliases against the circuit width: whenever `n` is a multiple of 3 every qubit rotates about one fixed axis in every layer, the circuit stops scrambling, and the measured variance jumped by three orders of magnitude at n = 6 and n = 12 (non-monotonic, and stable across seeds, so not sampling noise). Drawing the axes randomly per gate removed it and gave the clean 2^(−1.94n) fit above. Had I not checked monotonicity I would have reported a baseline fit of 2^(−1.03n) — off by a factor of two in the exponent.
* **Honest limitation:** the N = 200 global-cost figure is **not a measurement**; 2²⁰⁰ amplitudes cannot be simulated. It is read off the fitted exponential and the report carries `unstructured_global_grad_var_is_extrapolated: True` together with the fit sizes and residual, so the 111 orders cannot be mistaken for an observation. The equivariant side at N = 200 **is** a direct measurement, because the symmetric subspace has dimension 201.

| 6e | C3 (D3-b), Floquet | `floquet_hamiltonian.py`: added `bimodal_drive_hamiltonian` (the actual time-dependent drive, whose Fourier components are verified numerically to equal the `H_0` and `H_±1` the Magnus expansion uses), `exact_one_period_propagator` (midpoint product integrator, unitary to 1e-10, converged in substeps), and `magnus_propagator_fidelity` — a real gate fidelity `F = \|Tr(U_magnus† U_exact)\|²/4` from state evolution. `closed_loop_pipeline.py`: the recorded fidelity is now evaluated at the **CBF-applied** field, not the nominal one. Test `tests/test_floquet_propagator_fidelity.py` written first | `H_eff` was computed by `compute_floquet_magnus_hamiltonian` and then **discarded** — `dual_task_fidelity` is a weighted score of two step-out ratios that never touches the Hamiltonian. The pipeline recorded fidelity from the nominal command, so the CBF ablation compared **the identical number** to itself (0.992558362332366 in both arms) | Magnus fidelity **0.99999695** at the nominal 8.5/6.0 mT, measured. Pipeline now records 0.993464 (CBF on) vs 0.993391 (CBF off) at d = 15 µm — the arms finally differ. Suite 65 passed, 5 xfailed | `D3 Floquet:` |

**Findings from step 6e:**
* **The Floquet-Magnus engineering itself is sound, and that is now measured rather than assumed.** The truncated `H_eff` reproduces the exact one-period propagator to `F = 0.99999695`. The expansion is well inside its convergence regime (`‖H‖·T ≈ 0.14`), and the fidelity degrades as expected when the drive frequency is lowered at fixed amplitude, which is the behaviour a high-frequency expansion must show. This is a genuine positive result.
* **The `floquet_truncation_order: 2` setting does nothing.** Orders 1 and 2 give the same fidelity to 8 decimal places (0.99999695 both). The second-order term `h_2 = 0.15·(f_β/f_α)²·B_α·B_β·σ_z/ω_α` carries an arbitrary coefficient 0.15 that is not derived from the Magnus series, and it is numerically irrelevant at these parameters. The "2nd-order Floquet-Magnus" description is doing no work; it should either be derived properly or dropped to first order and described as such.
* **A SIXTH preregistered threshold is affected, and this one is a validity problem rather than a numerical one.** The preregistered `dual_task_fidelity ≥ 0.92` has **no term penalizing an inadequate field**, so it is monotonically non-increasing in field amplitude:

| field scale | 0.10 | 0.30 | 0.50 | 0.80 | 0.93 | 1.00 | 1.20 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `dual_task_fidelity` | **0.994371** | 0.994371 | 0.994371 | 0.994371 | 0.993441 | 0.992558 | 0.989673 |
| `alpha_rotation_rate_Hz` | 80.0 | 80.0 | 80.0 | 80.0 | 80.0 | 80.0 | 80.0 |

  A field at **10 % of nominal — which would do no clot drilling whatsoever — scores higher (0.994371) than the full nominal field (0.992558)**. The score is also completely blind to any reduction below 0.85×, because the effective step-out frequency is `f_stepout · clip(scale, 0.85, 1.25)`, and `alpha_rotation_rate_Hz` is pinned at 80 Hz at every amplitude because the step-out law only compares frequencies. **The metric rewards switching the actuator off.** Recorded as `xfail(strict=True)` on a test asserting the property a valid score must have; the definition is left untouched.
* Consequence for the CBF ablation, now visible for the first time: with the applied field recorded, the filtered arm scores *higher* (0.993464 vs 0.993391 at d = 15 µm). That difference is real but it is an artifact of the point above, and it must not be read as the CBF improving actuation fidelity. Any claim that the CBF preserves dual-task fidelity is currently unfalsifiable, because the score cannot go down when the field is cut.
