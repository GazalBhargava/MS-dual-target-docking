# Automated Molecular Docking Pipeline (2-Phase Architecture)

Current treatments for Multiple Sclerosis (MS) often rely on single-target immunosuppression. This project introduces a modular, high-throughput computational pipeline designed to identify novel drug alternatives for MS. By screening against the primary HLA-DRB1 presentation cleft and double-verifying those candidates against the CD80 and CD86 co-stimulatory interfaces, this workflow identifies repurposing candidates capable of a dual-action blockade to prevent the autoimmune cascade.

---

## 📌 Architecture Overview

The pipeline operates in two tightly coupled stages:

```
                      [ Prepared Ligand Library ]
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│  PHASE 1: High-Throughput Virtual Screening (HTVS)                 │
│  • Target: HLA-DRB1 (MHC-II antigen presentation cleft)            │
│  • Pocket: P4 pocket centroid [48.300, -7.509, 159.372]             │
│  • Output: All affinities & Top-10 repurposing candidates           │
└──────────────────────────────────┬──────────────────────────────────┘
                                   │  (Top-10 Candidates)
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│  PHASE 2: Targeted Interface Docking (Co-stimulatory Blockade)      │
│  • Target 1: CD86 (PDB: CD86Cleaned1NCN.pdb)                       │
│    Interface residues: [35, 37, 42, 43, 44, 47, 50, 52, 54,        │
│                         80, 82, 85, 87]                             │
│    Grid Centroid: [15.968, 31.116, 37.245]                          │
│                                                                     │
│  • Target 2: CD80 (PDB: CD80Cleaned1DR9.pdb)                       │
│    Interface residues: [54, 56, 58, 62, 63, 64, 67, 69, 74, 77,     │
│                         111, 113, 120, 122]                         │
│    Grid Centroid: [21.363, 30.076, 56.619]                          │
│                                                                     │
│  • Output: Docked pose PDBQT files & Comprehensive 3-target summary │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 🚀 Key Features

- **Real-Time Dynamic Status Bar:** Features a custom CLI progress tracker displaying a smoothed progress bar, live item throughput/min, elapsed time, and ETA calculations. 
  *(Example: `[██████████████████░░░░░░░] 75.0% | ETA: 00:05:30 | Speed: 12.4/min`)*
- **Interruption Recovery (Checkpointing):** Atomic saving to `docking_checkpoint.pkl` allows resuming multi-hour runs without lost progress.
- **Dynamic Centroid Calculation:** Geometric centroid computed directly from 3D alpha/beta atomic coordinates of specific interface residues using `Bio.PDB`.
- **Automatic Receptor Preparation:** On-the-fly PDB-to-PDBQT conversion via Meeko (`mk_prepare_receptor`).

---

## 📂 Project Structure

```text
├── main.py                         # Unified 2-Phase pipeline orchestrator
├── requirements.txt                # Python dependencies
├── .gitignore                      # Excluded data/files config
├── vina.exe                        # (User-provided) AutoDock Vina standalone binary
├── vina.py                         # Custom AutoDock Vina Python wrapper
├── README.md                       # Documentation
│
├── docking_results/                # (Auto-generated) Final CSV reports and rankings
├── targeted_docked_poses_cd86/     # (Auto-generated) Output 3D docked poses for CD86
└── targeted_docked_poses_cd80/     # (Auto-generated) Output 3D docked poses for CD80

*(Note: Target receptors, ligand libraries, and output folders are generated locally during execution and are excluded from version control).*
---

## 📊 Summary of Final Findings (Top 10 Candidates)

| Rank | Candidate Ligand | HLA-DRB1 (kcal/mol) | CD86 (kcal/mol) | CD80 (kcal/mol) | Primary Mechanism / Profile |
| :---: | :--- | :---: | :---: | :---: | :--- |
| **1** | **NS-11021** | **-10.46\*** | -5.77 | -6.76 | Potent HLA-DRB1 primary binder |
| **2** | **sb-218078** | **-10.19** | -6.04 | -7.28 | Chk1 kinase inhibitor |
| **3** | **afoxolaner** | **-10.05** | -6.22 | -6.77 | GABA-gated chloride channel antagonist |
| **4** | **CCG-63808** | -9.93 | -6.31 | -6.18 | RGS4 inhibitor |
| **5** | **AMG900** | -9.90 | -6.41 | -7.09 | Pan-Aurora kinase inhibitor |
| **6** | **deltarasin** | -9.89 | -6.10 | -6.55 | KRAS-PDEδ interaction inhibitor |
| **7** | **oxytetracycline** | -9.81 | -5.16 | -6.04 | Broad-spectrum antibiotic |
| **8** | **MK-4074** | -9.78 | -5.57 | -5.68 | Acetyl-CoA carboxylase inhibitor |
| **9** | **MK-3207** | -9.76 | **-6.85\*** | **-8.08\*** | **Best overall co-stimulation blocker** |
| **10**| **BT-11** | -9.75 | **-6.74** | **-7.93** | **Potent LANCL2 oral immunomodulator** |


*\* Denotes highest affinity in column.*

---

## 🛠️ Installation & Environment Setup

### 1. Requirements
- Python 3.10 or 3.11
- AutoDock Vina binary (`vina.exe` included)
- Microsoft Visual C++ Redistributable (Windows)

### 2. Virtual Environment Setup
```powershell
# Create virtual environment
python -m venv .venv

# Activate virtual environment
.\.venv\Scripts\Activate.ps1

# Install core dependencies
pip install numpy biopython meeko pandas rdkit
```

---

## 💻 Usage & CLI Reference

### 1. Run Complete Pipeline (Phase 1 + Phase 2)
Screens all ligands against HLA-DRB1, ranks the top candidates, and docks them into CD86 and CD80 interfaces:
```powershell
.venv\Scripts\python.exe main.py
```

### 2. Jump Straight to Phase 2 (Targeted Docking)
Carries over existing Top 10 candidates directly to interface docking:
```powershell
.venv\Scripts\python.exe main.py --phase2-only
```

### 3. Run Phase 2 for a Single Target
```powershell
# Dock only against CD86
.venv\Scripts\python.exe main.py --phase2-only --target cd86

# Dock only against CD80
.venv\Scripts\python.exe main.py --phase2-only --target cd80
```

### 4. Custom CPU & Search Parameters
```powershell
.venv\Scripts\python.exe main.py --cpu 8 --targeted-exhaustiveness 16 --top-n 15
```

---

## 📄 License & Attribution
Designed for computational structural biology & drug repurposing workflows. AutoDock Vina is developed by the Center for Computational Structural Biology (CCSB) at Scripps Research.
