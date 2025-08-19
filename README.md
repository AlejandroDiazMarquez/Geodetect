# Geodetect

# CIF Porosity Analyzer

<img width="3405" height="1784" alt="Picture1" src="https://github.com/user-attachments/assets/1385a158-7e32-4d8a-a110-9decad306a98" />

A Python toolkit to quantify porosity in crystalline materials from **CIF** files (Guest MOF CAU-55 with defects).  
It builds a 3D grid inside the (optionally expanded) unit cell, excludes points that fall inside atomic radii + probe size, clusters the remaining **void points** into **pore communities** via graph theory, estimates **pore volumes** with 3D **alpha shapes**, and (optionally) produces a **z-slice GIF** showing pore evolution across the structure.

- **Input:** CIF files in `./CIFs/` and a `radius` database file  
- **Output:** PDB converted from CIF, point grid (`grid.xyz`), OBJ meshes of pores, a summary text report, and optional frames + GIF

> Written by **Alejandro Díaz-Márquez (CNRS Montpellier)** — Version **1.0 (2025-08-01)**

---

## Features

- ✅ **CIF → PDB** conversion using **ASE**
- ✅ **Periodic expansion** (–1, 0, +1 in each direction) to capture border-spanning pores
- ✅ **Grid-based void sampling** with configurable spacing `g` and probe radius `size`
- ✅ **Pore community detection** via **NetworkX** (greedy modularity)
- ✅ **Alpha-shape volumetrics** (3D) to estimate per-community pore volumes
- ✅ **Graph metrics** per pore (clustering coefficient, assortativity, graph diameter, avg. betweenness)
- ✅ **Visualization**: colored OBJ meshes per pore; optional **z-slice GIF**
- ✅ Processes **all** `.cif` files in `./CIFs` automatically

---

## How it works (pipeline)

1. **Load + convert** each `*.cif` to PDB with ASE, read cell lengths/angles.
2. **Optionally expand** atoms to neighbor images (–1…+1 in x,y,z) for PBC awareness.
3. **3D grid sampling** inside the unit cell with spacing `g`; reject points within
   `atom_radius + probe_size` of any atom (KD-tree accelerated).
4. **Graph construction** on remaining “void points” with edges between points within `2g`.
5. **Community detection** (greedy modularity) → pores.
6. **Alpha-shape** (3D) per community → **pore volume** and **OBJ** mesh.
7. **Graph metrics** per community (NetworkX).
8. Optional **z-slice GIF** of pore points colored by community.

---

## Requirements

- **Python** ≥ 3.9 (recommended)
- **Packages**
  - `numpy`, `tqdm`, `matplotlib`, `mdtraj`, `networkx`, `scipy`, `ase`, `alphashape`, `trimesh`, `Pillow`
- **File dependencies**
  - `./CIFs/` directory containing input `.cif` files
  - `radius` text file with element radii (details below)

### Quick install (conda)

```bash
conda create -n cif-porosity python=3.11 -y
conda activate cif-porosity

# scientific stack
pip install numpy tqdm matplotlib mdtraj networkx scipy ase pillow

# geometry libraries
pip install alphashape trimesh
