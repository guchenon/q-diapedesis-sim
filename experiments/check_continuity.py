"""
Data-continuity audit: every number quoted in the docs/ELN must match what the code produces,
and every preregistered threshold must still hold its registered value.

Run: python -I experiments/check_continuity.py
Exit code 0 if the chain is consistent, 1 otherwise.
"""

import re
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

FAILURES: list[str] = []
CHECKS = 0


def check(label: str, ok: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  ok    {label}")
    else:
        print(f"  FAIL  {label}  {detail}")
        FAILURES.append(f"{label} {detail}")


# ---------------------------------------------------------------- 1. preregistered thresholds
PREREGISTERED = {
    ("vessel_mra_config.yaml", ("hemodynamics_pries_secomb", "cpp_baseline_mmHg")): 70.0,
    ("vessel_mra_config.yaml", ("hemodynamics_pries_secomb", "discharge_hematocrit")): 0.45,
    ("vessel_mra_config.yaml", ("thrombus_and_microenvironment", "post_occlusive_pH_initial")): 6.65,
    ("vessel_mra_config.yaml", ("thrombus_and_microenvironment", "physiological_pH_target")): 7.35,
    ("vessel_mra_config.yaml", ("biomimicry_layer", "selectin_catch_bond", "critical_force_pN")): 15.0,
    ("vessel_mra_config.yaml", ("biomimicry_layer", "hydro_magnetic_soliton", "washout_loss_ceiling_percent")): 2.0,
    ("quantum_floquet_config.yaml", ("swarm_scaling", "n_super_particles")): 200,
    ("quantum_floquet_config.yaml", ("bimodal_floquet_drive", "alpha_drill_frequency_Hz")): 80.0,
    ("quantum_floquet_config.yaml", ("bimodal_floquet_drive", "beta_trap_frequency_Hz")): 5.0,
    ("quantum_floquet_config.yaml", ("bimodal_floquet_drive", "target_dual_task_fidelity_min")): 0.92,
    ("quantum_floquet_config.yaml", ("lindblad_and_tensor_network", "ttn_max_bond_dimension_chi")): 16,
    ("neurocritical_cbf_config.yaml", ("cbf_safety_thresholds", "max_wall_shear_stress_Pa")): 5.50,
    ("neurocritical_cbf_config.yaml", ("cbf_safety_thresholds", "min_glycocalyx_integrity_ratio")): 0.92,
    ("neurocritical_cbf_config.yaml", ("cbf_safety_thresholds", "max_brain_temp_elevation_C")): 0.30,
}

print("1. Preregistered thresholds still hold their registered values")
cfgs = {
    name: yaml.safe_load((ROOT / "configs" / name).read_text(encoding="utf-8"))
    for name in ("vessel_mra_config.yaml", "quantum_floquet_config.yaml", "neurocritical_cbf_config.yaml")
}
for (fname, path), expected in PREREGISTERED.items():
    node = cfgs[fname]
    for key in path:
        node = node[key]
    if expected is None:
        check(f"{fname}:{'.'.join(path)} present", node is not None, repr(node))
    else:
        check(f"{fname}:{'.'.join(path)} == {expected}", node == expected, f"found {node!r}")

# ------------------------------------------------- 2. the tagged snapshot has not been rewritten
print("\n2. The v0.1.0-alpha preregistration tag is untouched")
try:
    sha = subprocess.run(
        ["git", "rev-list", "-n", "1", "v0.1.0-alpha"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()
    check("tag v0.1.0-alpha points at b5899d7", sha.startswith("b5899d7"), f"points at {sha[:10]}")
except Exception as exc:                                        # tag not fetched in this clone
    check("tag v0.1.0-alpha resolvable", False, str(exc))

# -------------------------------------------- 3. derived parameters agree with what the ELN says
print("\n3. Derived parameters reproduce the values quoted in ELN #004")
from src.physics import derive_dissipation_gamma, derive_swarm_coupling_K   # noqa: E402

gam = derive_dissipation_gamma(cfgs["vessel_mra_config.yaml"], cfgs["quantum_floquet_config.yaml"])
check("derived moment = 8.935e-16 A m^2", abs(gam["moment_A_m2"] - 8.935e-16) < 2e-19, f"{gam['moment_A_m2']:.4e}")
check("derived gamma (code units) = 6.711", abs(gam["viscous_dissipation_gamma"] - 6.711) < 0.01,
      f"{gam['viscous_dissipation_gamma']:.4f}")
check("mT -> rad/s conversion = 70.22", abs(gam["mT_to_rad_s"] - 70.22) < 0.05, f"{gam['mT_to_rad_s']:.3f}")
for d_um, expected_k in ((15.0, -2.22), (30.0, -1.65), (50.0, -0.55)):
    k = derive_swarm_coupling_K(cfgs["vessel_mra_config.yaml"], cfgs["quantum_floquet_config.yaml"], d_um)
    check(f"derived K at d={d_um:.0f} um = {expected_k}", abs(k - expected_k) < 0.05, f"{k:+.4f}")
    check(f"derived K at d={d_um:.0f} um is negative", k < 0.0, f"{k:+.4f}")

# --------------------------------------- 4. the models actually consume the derived values
print("\n4. The models consume the derived values (no silent fallback to the free numbers)")
from src.biomimicry_mechanics import HydroMagneticSolitonSwarm                # noqa: E402
from src.quantum_controller import LindbladEquivariantTTN                     # noqa: E402

swarm = HydroMagneticSolitonSwarm()
check("swarm coupling is derived", swarm.coupling_is_derived is True)
check("swarm coupling == derived value", swarm.coupling_K == swarm.derived_coupling_K_rad_s,
      f"{swarm.coupling_K} vs {swarm.derived_coupling_K_rad_s}")
check("swarm coupling != the legacy 4.50", swarm.coupling_K != 4.50, f"{swarm.coupling_K}")
ttn = LindbladEquivariantTTN()
check("dissipation is derived", ttn.dissipation_is_derived is True)
check("gamma != the legacy 0.15", ttn.gamma_viscous != 0.15, f"{ttn.gamma_viscous}")

# ------------------------------------------------ 5. every xfail reason cites an ELN entry
print("\n5. Every xfail records a reason that cites an ELN entry")
for path in sorted((ROOT / "tests").glob("test_*.py")):
    text = path.read_text(encoding="utf-8")
    for match in re.finditer(r"@pytest\.mark\.xfail\((.*?)\n\)", text, re.S):
        block = match.group(1)
        check(f"{path.name}: xfail cites an ELN entry", "ELN #" in block, block[:60].replace("\n", " "))
        check(f"{path.name}: xfail is strict", "strict=True" in block)

# ----------------------------------------------- 6. docs reference only files that exist
print("\n6. Paths named in README exist")
readme = (ROOT / "README.md").read_text(encoding="utf-8")
named = set(re.findall(r"`([a-zA-Z0-9_./-]+\.(?:py|yaml|yml|md|cff))`", readme))
# Also the repository-structure tree block, which is where the stale notebooks/ entry lived.
# Plain slicing rather than a regex: the block delimiter is a triple backtick fence and
# escaping it through several layers of quoting is how this check silently matched nothing.
fence = chr(96) * 3
heading = "## Repository Structure"
if heading in readme:
    after = readme.split(heading, 1)[1]
    parts = after.split(fence, 2)
    if len(parts) >= 3:
        # Resolve nested entries against their parent directory, from the tree indentation,
        # so "vessel_mra_config.yaml" is checked as "configs/vessel_mra_config.yaml".
        stack: list[str] = []
        for line in parts[1].splitlines():
            body = line.split("#")[0].rstrip()
            marker = max(body.find("\u251c\u2500\u2500"), body.find("\u2514\u2500\u2500"))
            if marker < 0:
                continue
            depth = marker // 4
            name = body[marker + 3:].strip()
            if not name:
                continue
            stack = stack[:depth]
            rel = "/".join(stack + [name.rstrip("/")])
            if name.endswith("/"):
                stack = stack[:depth] + [name.rstrip("/")]
            named.add(rel)
for rel in sorted(named):
    check(f"README path {rel}", (ROOT / rel).exists())

# ------------------------------------------------- 7. identifiers agree across every file
print("\n7. Identifiers are consistent across README, CITATION.cff, docs and source headers")
ORCID = "0009-0008-1357-9003"
OSF = "osf.io/htswk"
CONCEPT_DOI = "10.5281/zenodo.23599019"
VERSION_DOI = "10.5281/zenodo.23599020"

citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
check("CITATION.cff carries the ORCID", ORCID in citation)
check("CITATION.cff carries the concept DOI", CONCEPT_DOI in citation)
check("CITATION.cff carries the version DOI", VERSION_DOI in citation)
check("CITATION.cff carries repository-code", "repository-code:" in citation)
check("README carries the ORCID", ORCID in readme)
check("README carries the OSF link", OSF in readme)

# No file may cite an ORCID or OSF GUID other than the declared ones.
import itertools
for path in itertools.chain((ROOT / "src").rglob("*.py"), (ROOT / "docs").rglob("*.md"), [ROOT / "README.md"]):
    text = path.read_text(encoding="utf-8", errors="ignore")
    for found in set(re.findall(r"orcid\.org/([0-9X-]{19})", text)):
        check(f"{path.name}: ORCID {found} is the declared one", found == ORCID)
    for found in set(re.findall(r"osf\.io/([a-z0-9]+)", text)):
        check(f"{path.name}: OSF GUID {found} is the declared one", found == "htswk")
    for found in set(re.findall(r"10\.5281/zenodo\.(\d+)", text)):
        check(
            f"{path.name}: zenodo.{found} is a declared DOI",
            found in {CONCEPT_DOI.split(".")[-1], VERSION_DOI.split(".")[-1]},
        )

print(f"\n{'=' * 70}\n{CHECKS - len(FAILURES)}/{CHECKS} continuity checks passed")
if FAILURES:
    print(f"{len(FAILURES)} FAILED:")
    for f in FAILURES:
        print(f"  - {f}")
sys.exit(1 if FAILURES else 0)
