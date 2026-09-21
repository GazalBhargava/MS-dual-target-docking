"""
Molecular Docking Automation Pipeline
======================================
Unified 2-Phase Screening & Targeted Docking:
  Phase 1 → High-Throughput Virtual Screening (HTVS) against HLA-DRB1.
  Phase 2 → Targeted Interface Docking against CD80 & CD86 CD28-binding pockets.
"""

import argparse
import collections
from datetime import datetime
import os
import pickle
import shutil
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
from Bio.PDB import PDBParser
from vina import Vina


# ==============================================================================
#  CONFIGURATION
# ==============================================================================
class Config:
    # Computational Resources
    CPU_COUNT = min(os.cpu_count() or 4, 8)
    LIGAND_DIR = "ligands_ready"

    # Phase 1: HTVS against primary target
    PRIMARY_RECEPTOR = "HLA-DRB1_prepared.pdbqt"
    PRIMARY_GRID_CENTER = [48.3, -7.509, 159.372]  # Centroid of P4 pocket residues [78, 26, 70, 71, 28, 13]
    PRIMARY_GRID_SIZE = [20.0, 20.0, 20.0]
    PRIMARY_EXHAUSTIVENESS = 8

    # Phase 2: Targeted Interface Docking
    TOP_N = 10
    TARGETED_EXHAUSTIVENESS = 8
    TARGETED_TARGETS = {
        "CD86": {
            "pdb": "CD86Cleaned1NCN.pdb",
            "pdbqt": "CD86Cleaned1NCN.pdbqt",
            "chain": "A",
            "pocket_residues": [56, 58, 67, 69, 122],
            "box_size": [20.0, 20.0, 20.0],
            "poses_dir": "targeted_docked_poses_cd86",
        },
        "CD80": {
            "pdb": "CD80Cleaned1DR9.pdb",
            "pdbqt": "CD80Cleaned1DR9.pdbqt",
            "chain": "A",
            "pocket_residues": [54, 56, 58, 62, 63, 64, 67, 69, 74, 77, 111, 113, 120, 122],
            "box_size": [20.0, 20.0, 20.0],
            "poses_dir": "targeted_docked_poses_cd80",
        },
    }

    # Fallback Candidates (used if Phase 1 is skipped before completion)
    DEFAULT_TOP_CANDIDATES = [
        ("NS-11021", -10.459),
        ("sb-218078", -10.188),
        ("afoxolaner", -10.049),
        ("CCG-63808", -9.928),
        ("AMG900", -9.901),
        ("deltarasin", -9.895),
        ("oxytetracycline", -9.811),
        ("MK-4074", -9.780),
        ("MK-3207", -9.764),
        ("BT-11", -9.753),
    ]

    # Checkpointing & Outputs
    CHECKPOINT_FILE = "docking_checkpoint.pkl"
    CHECKPOINT_INTERVAL_SEC = 30 * 60
    CHECKPOINT_INTERVAL_COUNT = 100
    OUTPUT_DIR = "docking_results"
    HTVS_CSV = os.path.join(OUTPUT_DIR, "htvs_all_affinities.csv")
    TOP10_CSV = os.path.join(OUTPUT_DIR, "top10_candidates.csv")
    TARGETED_CSV = os.path.join(OUTPUT_DIR, "targeted_docking_results.csv")
    COMPREHENSIVE_CSV = os.path.join(OUTPUT_DIR, "comprehensive_docking_summary.csv")


# ==============================================================================
#  PROGRESS TRACKER & CHECKPOINTING
# ==============================================================================
class ProgressTracker:
    """Real-time progress bar with EMA-smoothed ETA and throughput speed."""

    def __init__(self, total: int, completed: int = 0, name: str = "Docking"):
        self.total = max(total, 1)
        self.completed = completed
        self.start_time = time.time()
        self.last_time = self.start_time
        self.recent_speeds: collections.deque = collections.deque(maxlen=30)
        self.bar_len = 30

    def update(self, step: int = 1) -> None:
        now = time.time()
        dur = now - self.last_time
        self.last_time = now
        self.completed += step
        if dur > 0:
            self.recent_speeds.append(dur)

    @staticmethod
    def format_time(seconds: float) -> str:
        if seconds < 0 or not np.isfinite(seconds):
            return "--:--:--"
        s = int(seconds)
        h, s = divmod(s, 3600)
        m, s = divmod(s, 60)
        return f"{h:02d}:{m:02d}:{s:02d}"

    def render(self, label: str = "", score_str: str = "") -> None:
        elapsed = time.time() - self.start_time
        pct = (self.completed / self.total) * 100.0
        filled = min(int(self.bar_len * pct / 100.0), self.bar_len)
        bar = "█" * filled + "░" * (self.bar_len - filled)

        avg_dur = float(np.median(self.recent_speeds)) if self.recent_speeds else 0.0
        remaining = max(0, self.total - self.completed)
        eta_sec = avg_dur * remaining if avg_dur > 0 else 0.0
        speed_str = f"{60.0 / avg_dur:4.1f}/min" if avg_dur > 0 else "--.-/min"

        lbl = (label[:20] + "..") if len(label) > 22 else label
        sys.stdout.write(
            f"\r[{bar}] {pct:5.1f}% | {self.completed}/{self.total} | "
            f"ETA: {self.format_time(eta_sec)} | Speed: {speed_str} | "
            f"Elapsed: {self.format_time(elapsed)} | {lbl:<22} | {score_str}"
        )
        sys.stdout.flush()

    def finish(self) -> None:
        sys.stdout.write("\n")
        sys.stdout.flush()


class CheckpointManager:
    """Atomic save and reload for screening jobs."""

    def __init__(self, path: str = Config.CHECKPOINT_FILE):
        self.path = path
        self.last_time = time.time()
        self.last_count = 0

    def load(self) -> Tuple[Dict[str, float], Set[str]]:
        if os.path.exists(self.path):
            try:
                with open(self.path, "rb") as f:
                    data = pickle.load(f)
                res = data.get("results", {})
                done = set(data.get("completed", set())) | set(res.keys())
                print(f"[Checkpoint] Resumed previous run: {len(res)} ligands completed.")
                return res, done
            except Exception as e:
                print(f"[Checkpoint Warning] Failed to read '{self.path}': {e}. Starting fresh.")
        return {}, set()

    def save(self, results: Dict[str, float], completed: Set[str], force: bool = False) -> bool:
        now = time.time()
        count = len(results)
        if not force:
            if (now - self.last_time < Config.CHECKPOINT_INTERVAL_SEC and
                    count - self.last_count < Config.CHECKPOINT_INTERVAL_COUNT):
                return False

        tmp = f"{self.path}.tmp_{os.getpid()}"
        try:
            with open(tmp, "wb") as f:
                pickle.dump({"results": results, "completed": completed, "count": count, "timestamp": now}, f)
            os.replace(tmp, self.path)
            self.last_time = now
            self.last_count = count
            return True
        except Exception as e:
            print(f"\n[Checkpoint Error] {e}")
            if os.path.exists(tmp):
                os.remove(tmp)
            return False


# ==============================================================================
#  STRUCTURAL & DOCKING HELPERS
# ==============================================================================
def prepare_receptor(pdb_file: str, pdbqt_file: str) -> bool:
    """Ensure receptor PDBQT exists; generate via Meeko if missing."""
    if os.path.exists(pdbqt_file):
        return True

    print(f"  [Meeko] Generating '{pdbqt_file}' from '{pdb_file}'...")
    exe = (shutil.which("mk_prepare_receptor.exe") or
           shutil.which("mk_prepare_receptor") or
           os.path.join(os.path.dirname(sys.executable), "mk_prepare_receptor.exe"))

    if exe and os.path.exists(exe):
        proc = subprocess.run([exe, "--read_pdb", pdb_file, "--write_pdbqt", pdbqt_file],
                              capture_output=True, text=True)
        if proc.returncode == 0 and os.path.exists(pdbqt_file):
            print(f"  [Meeko] Prepared '{pdbqt_file}' successfully.")
            return True
        print(f"  [Meeko Error] {proc.stderr}")
    return False


def compute_pocket_centroid(pdb_file: str, chain: str = "A", residues: Optional[List[int]] = None) -> Optional[List[float]]:
    """Calculate the 3D arithmetic mean of residue coordinates defining a binding site."""
    if not os.path.exists(pdb_file):
        return None

    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure("receptor", pdb_file)
    except Exception:
        return None

    res_set = set(residues) if residues else set()
    coords = []

    # Priority 1: Specified chain
    for model in struct:
        for ch in model:
            if ch.id == chain:
                for res in ch:
                    if not res_set or res.id[1] in res_set:
                        coords.extend(atom.get_coord() for atom in res)

    # Priority 2: Fallback across all chains if residue not matched
    if not coords:
        for model in struct:
            for ch in model:
                for res in ch:
                    if not res_set or res.id[1] in res_set:
                        coords.extend(atom.get_coord() for atom in res)

    if not coords:
        return None

    centroid = np.mean(coords, axis=0)
    return [round(float(c), 3) for c in centroid]


def load_candidate_pool(top_n: int = Config.TOP_N) -> List[Tuple[str, float]]:
    """Retrieve top candidates in priority order: top10 CSV -> checkpoint -> HTVS CSV -> default fallback."""
    for path in [Config.TOP10_CSV, Config.HTVS_CSV]:
        if os.path.exists(path):
            try:
                candidates = []
                with open(path, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("Rank"):
                            continue
                        parts = line.split(",")
                        if len(parts) >= 3:
                            candidates.append((parts[1].strip(' "'), float(parts[2].strip())))
                if candidates:
                    candidates.sort(key=lambda x: x[1])
                    print(f"[Candidates] Loaded {len(candidates[:top_n])} candidates from '{path}'.")
                    return candidates[:top_n]
            except Exception:
                pass

    if os.path.exists(Config.CHECKPOINT_FILE):
        ckpt = CheckpointManager(Config.CHECKPOINT_FILE)
        res, _ = ckpt.load()
        if res:
            sorted_res = sorted(res.items(), key=lambda x: x[1])[:top_n]
            print(f"[Candidates] Loaded Top {len(sorted_res)} candidates from checkpoint.")
            return sorted_res

    print("[Candidates] Using default top candidate pool fallback.")
    return Config.DEFAULT_TOP_CANDIDATES[:top_n]


# ==============================================================================
#  PHASE 1: HIGH-THROUGHPUT VIRTUAL SCREENING
# ==============================================================================
def run_phase1_htvs(
    ligand_dir: str = Config.LIGAND_DIR,
    receptor_file: str = Config.PRIMARY_RECEPTOR,
    center: List[float] = Config.PRIMARY_GRID_CENTER,
    box_size: List[float] = Config.PRIMARY_GRID_SIZE,
    cpu: int = Config.CPU_COUNT,
    exhaustiveness: int = Config.PRIMARY_EXHAUSTIVENESS,
) -> Tuple[Dict[str, float], List[Tuple[str, float]]]:
    """Execute high-throughput screening against the primary target (HLA-DRB1)."""
    print("\n" + "=" * 78)
    print("  PHASE 1 — HIGH-THROUGHPUT VIRTUAL SCREENING (HTVS)")
    print("=" * 78)

    if not os.path.exists(receptor_file):
        raise FileNotFoundError(f"Primary receptor file '{receptor_file}' not found.")
    if not os.path.exists(ligand_dir):
        raise FileNotFoundError(f"Ligand library folder '{ligand_dir}' not found.")

    ligand_files = [os.path.join(ligand_dir, f) for f in os.listdir(ligand_dir) if f.endswith(".pdbqt")]
    total = len(ligand_files)

    print(f"  Target Receptor : {receptor_file}")
    print(f"  Grid Center     : {center}")
    print(f"  Box Dimensions  : {box_size}")
    print(f"  Exhaustiveness  : {exhaustiveness} | CPU Cores: {cpu}")
    print(f"  Ligands to Dock : {total}\n")

    ckpt = CheckpointManager()
    results, completed = ckpt.load()

    v = Vina(sf_name="vina", cpu=cpu)
    v.set_receptor(receptor_file)
    v.compute_vina_maps(center=center, box_size=box_size)

    tracker = ProgressTracker(total=total, completed=len(completed), name="HTVS")
    print("  Screening active ligand library...\n")

    try:
        for lig_path in ligand_files:
            drug = os.path.basename(lig_path)[:-6]  # strip .pdbqt
            if drug in completed:
                continue

            try:
                v.set_ligand_from_file(lig_path)
                v.dock(exhaustiveness=exhaustiveness, n_poses=1)
                energies = v.energies(n_poses=1)
                score = energies[0][0] if energies else 0.0

                results[drug] = score
                completed.add(drug)
                tracker.update()
                tracker.render(label=drug, score_str=f"{score:6.2f} kcal/mol")

                if ckpt.save(results, completed):
                    sys.stdout.write(f"\n[Checkpoint Saved] Total: {len(results)}\n")
            except Exception as e:
                sys.stdout.write(f"\n[Warning] Docking failed for '{drug}': {e}\n")
                tracker.update()
    except KeyboardInterrupt:
        print("\n\n[Warning] Screening interrupted by user. Saving checkpoint...")
        ckpt.save(results, completed, force=True)
        sys.exit(0)

    tracker.finish()
    ckpt.save(results, completed, force=True)
    print("\n  [Done] Phase 1 screening complete.")

    top_candidates = sorted(results.items(), key=lambda x: x[1])[:Config.TOP_N]
    return results, top_candidates


# ==============================================================================
#  PHASE 2: TARGETED INTERFACE DOCKING
# ==============================================================================
def run_phase2_targeted_docking(
    candidates: List[Tuple[str, float]],
    ligand_dir: str = Config.LIGAND_DIR,
    cpu: int = Config.CPU_COUNT,
    exhaustiveness: int = Config.TARGETED_EXHAUSTIVENESS,
    target_filter: Optional[str] = None,
) -> Dict[str, Dict[str, float]]:
    """Perform targeted docking of top candidates into CD80 and CD86 interfaces."""
    print("\n" + "=" * 78)
    print("  PHASE 2 — TARGETED INTERFACE DOCKING (CD80 / CD86)")
    print("=" * 78)

    targets = Config.TARGETED_TARGETS
    if target_filter and target_filter.lower() != "all":
        targets = {k: v for k, v in targets.items() if target_filter.lower() in k.lower()}
        if not targets:
            print(f"  [Warning] Filter '{target_filter}' matched no targets. Defaulting to all.")
            targets = Config.TARGETED_TARGETS

    targeted_results: Dict[str, Dict[str, float]] = {}
    total_runs = len(targets) * len(candidates)
    tracker = ProgressTracker(total=total_runs, name="Targeted")

    for target_name, params in targets.items():
        targeted_results[target_name] = {}
        poses_dir = params["poses_dir"]
        os.makedirs(poses_dir, exist_ok=True)

        print(f"\n  ──── Target: {target_name} ────")
        pdb_file = params["pdb"]
        pdbqt_file = params["pdbqt"]

        if not os.path.exists(pdb_file):
            print(f"  [Skip] PDB structure '{pdb_file}' not found.")
            for _ in candidates:
                tracker.update()
            continue

        if not prepare_receptor(pdb_file, pdbqt_file):
            print(f"  [Skip] Receptor preparation failed for '{pdbqt_file}'.")
            for _ in candidates:
                tracker.update()
            continue

        center = compute_pocket_centroid(pdb_file, chain=params["chain"], residues=params["pocket_residues"])
        if not center:
            print(f"  [Skip] Could not compute centroid for {target_name}.")
            for _ in candidates:
                tracker.update()
            continue

        box = params["box_size"]
        print(f"  Receptor Structure : {pdb_file} -> {pdbqt_file}")
        print(f"  Interface Residues : {params['pocket_residues']}")
        print(f"  Calculated Centroid: {center}")
        print(f"  Grid Box Size      : {box}")
        print(f"  Poses Directory    : {poses_dir}/")
        print(f"  Docking top {len(candidates)} candidates...\n")

        try:
            v = Vina(sf_name="vina", cpu=cpu)
            v.set_receptor(pdbqt_file)
            v.compute_vina_maps(center=center, box_size=box)
        except Exception as e:
            print(f"  [Error] Vina grid initialization failed for {target_name}: {e}")
            for _ in candidates:
                tracker.update()
            continue

        for drug, _ in candidates:
            lig_file = os.path.join(ligand_dir, f"{drug}.pdbqt")
            if not os.path.exists(lig_file):
                sys.stdout.write(f"\n  [Warning] Ligand file '{lig_file}' missing. Skipping.\n")
                tracker.update()
                continue

            try:
                v.set_ligand_from_file(lig_file)
                v.dock(exhaustiveness=exhaustiveness, n_poses=1)
                energies = v.energies(n_poses=1)
                score = energies[0][0] if energies else 0.0
                targeted_results[target_name][drug] = score

                pose_out = os.path.join(poses_dir, f"{drug}_on_{target_name}_targeted.pdbqt")
                v.write_poses(pose_out, n_poses=1, overwrite=True)

                tracker.update()
                tracker.render(label=f"{target_name}: {drug}", score_str=f"{score:6.2f} kcal/mol")
            except Exception as e:
                sys.stdout.write(f"\n  [Warning] Targeted docking error for '{drug}': {e}\n")
                tracker.update()

    tracker.finish()
    print("\n  [Done] Phase 2 targeted docking complete.")
    return targeted_results


# ==============================================================================
#  DATA EXPORT & REPORTING
# ==============================================================================
def export_csv_summaries(
    htvs_results: Optional[Dict[str, float]],
    candidates: List[Tuple[str, float]],
    targeted_results: Dict[str, Dict[str, float]],
) -> None:
    """Save clean, standardized CSV reports for analysis and publication."""
    os.makedirs(Config.OUTPUT_DIR, exist_ok=True)

    # 1. Complete HTVS affinities
    if htvs_results:
        with open(Config.HTVS_CSV, "w", encoding="utf-8") as f:
            f.write("Rank,Ligand_Name,Affinity_kcal_mol\n")
            for r, (d, s) in enumerate(sorted(htvs_results.items(), key=lambda x: x[1]), 1):
                f.write(f'{r},"{d}",{s:.3f}\n')
        print(f"[Export] Saved {len(htvs_results)} HTVS scores to: {Config.HTVS_CSV}")

    # 2. Top-N candidate affinities
    with open(Config.TOP10_CSV, "w", encoding="utf-8") as f:
        f.write("Rank,Ligand_Name,HLA_DRB1_Affinity_kcal_mol\n")
        for r, (d, s) in enumerate(candidates, 1):
            f.write(f'{r},"{d}",{s:.3f}\n')
    print(f"[Export] Saved Top {len(candidates)} candidates to: {Config.TOP10_CSV}")

    # 3. Targeted docking and comprehensive combined summaries
    if targeted_results:
        targets = list(targeted_results.keys())
        for path in [Config.TARGETED_CSV, Config.COMPREHENSIVE_CSV]:
            with open(path, "w", encoding="utf-8") as f:
                header = ["Rank", "Ligand_Name", "HLA_DRB1_Affinity_kcal_mol"] + [f"Targeted_{t}_kcal_mol" for t in targets]
                f.write(",".join(header) + "\n")
                for r, (d, hla) in enumerate(candidates, 1):
                    row = [str(r), f'"{d}"', f"{hla:.3f}"]
                    for t in targets:
                        val = targeted_results[t].get(d)
                        row.append(f"{val:.3f}" if val is not None else "N/A")
                    f.write(",".join(row) + "\n")
        print(f"[Export] Saved targeted summary to: {Config.TARGETED_CSV}")
        print(f"[Export] Saved comprehensive summary to: {Config.COMPREHENSIVE_CSV}")


def print_summary_table(
    candidates: List[Tuple[str, float]],
    targeted_results: Dict[str, Dict[str, float]],
    elapsed_seconds: float,
) -> None:
    """Print an aligned, publication-ready binding affinity summary table."""
    targets = list(targeted_results.keys())
    col_w = 14
    name_w = 22
    total_w = 6 + name_w + (1 + len(targets)) * (col_w + 2)
    sep = "=" * total_w

    # Identify best score per target
    best_hla = min([s for _, s in candidates]) if candidates else None
    best_target: Dict[str, Optional[float]] = {}
    for t in targets:
        vals = [v for v in targeted_results.get(t, {}).values() if v is not None]
        best_target[t] = min(vals) if vals else None

    print("\n" + sep)
    print(f"  BINDING AFFINITY SUMMARY — Top {len(candidates)} Repurposing Candidates")
    print(f"  All values in kcal/mol | * = best score in target column | Elapsed: {ProgressTracker.format_time(elapsed_seconds)}")
    print(sep)

    phase_header = f"  {'':5}{'':>{name_w}}  {'── Phase 1 ──':>{col_w}}"
    if targets:
        phase_header += f"  {'── Phase 2 ──' + ' ' * ((col_w + 2) * len(targets) - 15)}"
    print(phase_header)

    hdr = f"  {'Rank':<5}{'Candidate Ligand':<{name_w}}  {'HLA-DRB1':>{col_w}}"
    for t in targets:
        hdr += f"  {t[:col_w]:>{col_w}}"
    print(hdr)
    print("  " + "-" * (total_w - 2))

    for rank, (drug, hla) in enumerate(candidates, 1):
        hla_mk = "*" if hla == best_hla else " "
        row = f"  {rank:<5}{drug:<{name_w}}{hla:>{col_w - 1}.2f}{hla_mk}"
        for t in targets:
            val = targeted_results.get(t, {}).get(drug)
            if val is not None:
                mk = "*" if val == best_target.get(t) else " "
                row += f"  {val:>{col_w - 1}.2f}{mk}"
            else:
                row += f"  {'N/A':>{col_w}}"
        print(row)

    print("  " + "-" * (total_w - 2))
    best_row = f"  {'':5}{'Best':>{name_w}}{best_hla:>{col_w - 1}.2f}*" if best_hla is not None else ""
    for t in targets:
        bv = best_target.get(t)
        best_row += f"  {bv:>{col_w - 1}.2f}*" if bv is not None else f"  {'N/A':>{col_w}}"
    print(best_row)
    print(sep + "\n")


# ==============================================================================
#  CLI ORCHESTRATOR
# ==============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Unified 2-Phase Screening & Targeted Docking Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ligand-dir", default=Config.LIGAND_DIR, help="Path to prepared ligand library (.pdbqt)")
    parser.add_argument("--receptor", default=Config.PRIMARY_RECEPTOR, help="Path to primary receptor (.pdbqt)")
    parser.add_argument("--cpu", type=int, default=Config.CPU_COUNT, help="Number of CPU cores for Vina")
    parser.add_argument("--exhaustiveness", type=int, default=Config.PRIMARY_EXHAUSTIVENESS, help="Exhaustiveness for Phase 1 HTVS")
    parser.add_argument("--targeted-exhaustiveness", type=int, default=Config.TARGETED_EXHAUSTIVENESS, help="Exhaustiveness for Phase 2")
    parser.add_argument("--target", default="all", help="Target selector for Phase 2 ('cd80', 'cd86', or 'all')")
    parser.add_argument("--top-n", type=int, default=Config.TOP_N, help="Number of candidates to carry into Phase 2")
    parser.add_argument("--skip-phase1", action="store_true", help="Skip Phase 1 and load candidates from CSV/checkpoint")
    parser.add_argument("--skip-phase2", action="store_true", help="Skip Phase 2 targeted docking")
    parser.add_argument("--phase2-only", action="store_true", help="Directly execute Phase 2 (skips Phase 1)")

    args = parser.parse_args()
    if args.phase2_only:
        args.skip_phase1 = True

    start_wall = time.time()
    print("=" * 78)
    print("  MOLECULAR DOCKING AUTOMATION PIPELINE (PHASE 1 & PHASE 2)")
    print(f"  Execution Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 78)

    # 1. Phase 1 HTVS
    htvs_results: Optional[Dict[str, float]] = None
    if not args.skip_phase1:
        htvs_results, candidates = run_phase1_htvs(
            ligand_dir=args.ligand_dir,
            receptor_file=args.receptor,
            center=Config.PRIMARY_GRID_CENTER,
            box_size=Config.PRIMARY_GRID_SIZE,
            cpu=args.cpu,
            exhaustiveness=args.exhaustiveness,
        )
    else:
        print("\n[Phase 1] Skipped. Resolving top candidates from existing data...")
        candidates = load_candidate_pool(top_n=args.top_n)

    if not candidates:
        print("[Error] No candidate ligands found. Exiting.")
        sys.exit(1)

    # 2. Phase 2 Targeted Docking
    targeted_results: Dict[str, Dict[str, float]] = {}
    if not args.skip_phase2:
        targeted_results = run_phase2_targeted_docking(
            candidates=candidates,
            ligand_dir=args.ligand_dir,
            cpu=args.cpu,
            exhaustiveness=args.targeted_exhaustiveness,
            target_filter=args.target,
        )
    else:
        print("\n[Phase 2] Skipped.")

    # 3. Export Data & Print Publication Summary
    print(">>> Exporting summary results...")
    export_csv_summaries(htvs_results, candidates, targeted_results)
    print_summary_table(candidates, targeted_results, elapsed_seconds=time.time() - start_wall)


if __name__ == "__main__":
    main()
