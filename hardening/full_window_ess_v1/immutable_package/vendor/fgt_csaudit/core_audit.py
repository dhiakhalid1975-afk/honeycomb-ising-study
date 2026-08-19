from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Any

import numpy as np

from .provenance import sha256_file

CORE_FILES = (
    "src/fgt_dilution/lattice.py",
    "src/fgt_dilution/disorder.py",
    "src/fgt_dilution/kernels.py",
    "src/fgt_dilution/model.py",
    "src/fgt_dilution/observables.py",
    "src/fgt_dilution/simulation.py",
    "src/fgt_dilution/fss.py",
    "src/fgt_dilution/crossing_bootstrap.py",
    "src/fgt_dilution/quenched_bootstrap.py",
)


def audit_monte_carlo_core(project_root: Path) -> dict[str, Any]:
    """Run a small, deterministic audit of the original v2.4.0 Monte Carlo core.

    This is not a replacement for the original test suite. It verifies the pieces that
    the critical-scaling audit relies on: honeycomb topology, vacancy representation,
    ferromagnetic Wolff dynamics, incremental E/M consistency, and finite observables.
    """
    project_root = Path(project_root).resolve()
    missing = [rel for rel in CORE_FILES if not (project_root / rel).exists()]
    if missing:
        raise FileNotFoundError("missing original Monte Carlo core files: " + ", ".join(missing))

    src = str(project_root / "src")
    inserted = False
    if src not in sys.path:
        sys.path.insert(0, src)
        inserted = True
    try:
        lattice_mod = importlib.import_module("fgt_dilution.lattice")
        disorder_mod = importlib.import_module("fgt_dilution.disorder")
        kernels_mod = importlib.import_module("fgt_dilution.kernels")
        model_mod = importlib.import_module("fgt_dilution.model")
        simulation_mod = importlib.import_module("fgt_dilution.simulation")

        checks: list[dict[str, Any]] = []

        def add(name: str, ok: bool, detail: Any = None) -> None:
            checks.append({"check": name, "pass": bool(ok), "detail": detail})

        L = 8
        lat = lattice_mod.honeycomb_lattice(L, J=1.0)
        degree = np.sum(lat.neighbors >= 0, axis=1)
        add("honeycomb_site_count_N_equals_2L2", lat.n_sites == 2 * L * L, int(lat.n_sites))
        add("honeycomb_coordination_three", bool(np.all(degree == 3)), sorted(set(int(x) for x in degree)))
        add("honeycomb_unique_bond_count_3L2", len(lat.bonds_i) == 3 * L * L, int(len(lat.bonds_i)))

        occ_all = np.ones(lat.n_sites, dtype=np.int8)
        e0 = model_mod.ground_state_energy_per_active_spin(occ_all, lat)
        add("ferromagnetic_ground_state_energy_per_spin_minus_3_over_2", abs(float(e0) + 1.5) < 1e-12, float(e0))

        rng = np.random.default_rng(123456)
        occ = disorder_mod.generate_occupancy(
            lat,
            disorder_mod.DisorderSpec(mode="bernoulli", p=0.80, exact_count=True),
            rng,
        )
        spins = disorder_mod.initialize_spins(occ, np.random.default_rng(98765), ordered=False)
        add("vacancies_are_zero_spins", bool(np.all(spins[occ == 0] == 0)), int(np.sum(occ == 0)))
        active = np.flatnonzero(occ).astype(np.int32)
        kernels_mod.seed_numba(24680)
        marks = np.zeros(lat.n_sites, dtype=np.int32)
        stack = np.empty(lat.n_sites, dtype=np.int32)
        token = 0
        e, m = kernels_mod.energy_magnetization_total(spins, lat.bonds_i, lat.bonds_j, lat.bonds_jval, 0.0)
        max_err = 0.0
        for _ in range(100):
            _, token, dm, de = kernels_mod.wolff_step_incremental(
                spins, active, lat.neighbors, lat.couplings, 1.0 / 1.2, marks, stack, token
            )
            e += float(de)
            m += float(dm)
            er, mr = kernels_mod.energy_magnetization_total(spins, lat.bonds_i, lat.bonds_j, lat.bonds_jval, 0.0)
            max_err = max(max_err, abs(float(er - e)), abs(float(mr - m)))
            e, m = float(er), float(mr)
        add("incremental_wolff_energy_magnetization_consistency", max_err < 1e-10, max_err)

        result = simulation_mod.run_single_realization(
            lat,
            1.5,
            disorder_mod.DisorderSpec(mode="bernoulli", p=1.0, exact_count=True),
            realization=0,
            equilibration_steps=40,
            measurement_steps=80,
            measurement_stride=2,
            algorithm="wolff",
            wolff_step_unit="sweep_equivalent",
            field=0.0,
            ordered_start=False,
            base_seed=13579,
            fast_incremental_observables=True,
            measure_structure_factor=True,
            energy_check_interval=10,
        )
        rec = result.record
        finite_keys = ("energy", "abs_m", "m2", "m4", "chi_abs", "binder")
        add("small_realization_observables_finite", all(np.isfinite(float(rec[k])) for k in finite_keys), {k: rec[k] for k in finite_keys})
        add("small_realization_p_actual_one", abs(float(rec["p_actual"]) - 1.0) < 1e-12, float(rec["p_actual"]))
        mismatch = float(rec.get("energy_max_abs_mismatch", np.nan))
        add("small_realization_incremental_check", np.isfinite(mismatch) and mismatch < 1e-9, mismatch)

        hashes = {rel: sha256_file(project_root / rel) for rel in CORE_FILES}
        passed = all(row["pass"] for row in checks)
        return {
            "status": "PASS" if passed else "FAIL",
            "passed": passed,
            "scope": "targeted_runtime_audit_not_full_revalidation_of_original_project",
            "checks": checks,
            "source_hashes": hashes,
        }
    finally:
        if inserted and sys.path and sys.path[0] == src:
            sys.path.pop(0)
