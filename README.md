# multiloops-classification

Classification of RNA multiloop (junction) structures into structural families using machine learning and geometric approaches.

## Overview

This repository contains the source code and datasets for a computational study on the classification of RNA multiloop junctions. RNA multiloops (also referred to as junctions) are structural motifs where three or more helical stems converge. Their geometry, governed largely by coaxial stacking interactions between helical arms, determines which structural family a junction belongs to. Accurate classification of these motifs is important for understanding RNA tertiary structure, function, and folding.

The project implements and compares three distinct approaches to junction classification:

1. **Thermodynamic feature-based classification** -- reproduction and optimization of the method proposed by Laing et al. (2012)
2. **Rule-based classification** -- reference implementation of the CarTAJ tool by Lamiable et al. (2012)
3. **3D geometry-based classification** -- a novel approach using coarse-grained spatial representations of RNA structure

The study addresses both **3-way junctions** (families A, B, C, defined by which pair of helical arms stacks coaxially) and **4-way junctions** (families H, cH, cL, cK, cW, cX, pi, psi, X).

## Background

RNA multiloops are among the most structurally diverse and functionally significant motifs in RNA molecules. A multiloop is formed at the point where three or more double-stranded helical stems meet, connected by single-stranded loop segments. The way these helical arms arrange themselves in three-dimensional space -- particularly through coaxial stacking, where two helices align end-to-end to form a continuous stack -- defines the junction's structural family.

For 3-way junctions, three families are recognized:
- **Family A** -- coaxial stacking between helices H1 and H2
- **Family B** -- coaxial stacking between helices H2 and H3
- **Family C** -- coaxial stacking between helices H1 and H3

For 4-way junctions, nine families are recognized (H, cH, cL, cK, cW, cX, pi, psi, X), reflecting the more complex combinatorial possibilities of four converging helices.

## Repository Structure

```
multiloops-classification/
├── data/                        # Curated datasets
│   ├── raw_article_data/        # Data manually transcribed from published literature
│   │   ├── laing_2009.csv       # 4-way junction data from Laing et al. (2009)
│   │   ├── laing_2012.csv       # 3-way and 4-way junction data from Laing et al. (2012)
│   │   └── lamiable_2012.csv    # 3-way junction data from Lamiable et al. (2012)
│   ├── 3way_junctions_labeled.csv   # Consolidated labeled dataset for 3-way junctions
│   ├── 4way_junctions_labeled.csv   # Consolidated labeled dataset for 4-way junctions
│   ├── literature_pdbs.txt      # PDB identifiers referenced in published studies
│   └── non_redundant_pdbs.txt   # Non-redundant set of PDB identifiers
│
├── laing_experiment/            # Thermodynamic feature-based classification
│   ├── features.ipynb           # Feature engineering (Turner energies, adenine counts, loop lengths)
│   ├── laing_experiments.ipynb  # Reproduction of the original Random Forest classifier
│   └── optimised_experiment.ipynb  # Hyperparameter optimization and classifier comparison
│
├── lamiable_experiment/         # Rule-based classification (CarTAJ reference)
│   └── cartaj_3w.js             # Original CarTAJ source code for 3-way junction classification
│
├── graph_approach/              # 3D geometry-based classification
│   ├── process_pdb_by_forgi.sh  # PDB-to-coarse-grained conversion pipeline
│   ├── forgi_tofeature_vectors.sh  # Feature vector extraction from forgi output
│   ├── parse_forgi.py           # Parser for forgi coarse-grained junction representations
│   ├── calculate_angles.py      # Geometric feature extraction (angles, triangle heights)
│   ├── graphs_classifier.py     # ML classification pipeline for geometric features
│   └── graph_classification.ipynb  # Interactive notebook for geometric classification
│
├── CITATION.cff                 # Citation metadata
├── LICENSE                      # MIT License
└── requirements.txt             # Python dependencies
```

## Methods

### Thermodynamic Feature-Based Classification (Laing Experiment)

This approach reproduces and extends the method described by Laing et al. (2012). Sixteen features are engineered from junction properties, including Turner (2004) thermodynamic stacking free energies, counts of consecutive adenines in loop segments, and sorted loop lengths. Multiple classifiers are evaluated -- including Random Forest, k-Nearest Neighbors, Support Vector Machines, Logistic Regression, XGBoost, and LightGBM -- with hyperparameter optimization via grid search and class imbalance handling through oversampling strategies (SMOTE, RandomOverSampler). Feature importance is analyzed using SHAP values. Both 3-way and 4-way junction classification are addressed.

### Rule-Based Classification (Lamiable Experiment)

This module contains the reference JavaScript source code of the CarTAJ tool (Lamiable et al., 2012), a deterministic (non-ML) classifier that uses isostericity matrices and base-pair geometry to assign 3-way junctions to families A, B, or C. It serves as a baseline for comparison with the data-driven approaches.

### 3D Geometry-Based Classification (Graph Approach)

This novel approach derives features directly from the three-dimensional coordinates of RNA structures. PDB files are processed through the forgi library to obtain coarse-grained representations of multiloop segments and stems. Junction topologies are then reconstructed by traversing connected structural elements. Six geometric features are extracted per junction: three inter-stem angles (computed by projecting stem direction vectors onto a plane defined by the junction centroid) and three triangle heights derived from the spatial arrangement of stems. These features are used to train classifiers including k-NN, SVM, and Random Forest.

## Requirements

- Python 3.10+
- Dependencies listed in `requirements.txt`

Key dependencies include:

| Package | Purpose |
|---|---|
| `pandas` | Data manipulation |
| `scikit-learn` | Machine learning classifiers and evaluation |
| `xgboost`, `lightgbm` | Gradient boosting classifiers |
| `imbalanced-learn` | Oversampling for class imbalance |
| `forgi` | RNA coarse-grained structure representation |
| `biopython` | PDB file parsing |
| `matplotlib`, `seaborn` | Visualization |

## Installation

```bash
git clone https://github.com/tzok/multiloops-classification.git
cd multiloops-classification
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Known Issues

The `forgi` library (v2.2.3) omits multiloop and interior loop coordinates from its coarse-grained output by default. To enable the graph-based approach, the following lines must be commented out in the installed package:

**File:** `.venv/lib/python3.x/site-packages/forgi/threedee/model/coarse_grain.py`

```python
if key[0] in ["m", "i"]:
    # Bulge coordinates are redundant. They can be deduced from the stem coordinates.
    continue
```

Comment out or remove this block to ensure that multiloop segment coordinates are included in the output.

## Usage

### Thermodynamic Feature-Based Classification

Open and execute the Jupyter notebooks in `laing_experiment/` sequentially:

1. `features.ipynb` -- computes features from the labeled junction dataset
2. `laing_experiments.ipynb` -- runs the original Random Forest classification experiment
3. `optimised_experiment.ipynb` -- runs classifier comparison with hyperparameter optimization

### 3D Geometry-Based Classification

1. Prepare coarse-grained representations from PDB files:
   ```bash
   bash graph_approach/process_pdb_by_forgi.sh
   bash graph_approach/forgi_tofeature_vectors.sh
   ```
2. Parse junctions and extract geometric features:
   ```bash
   python graph_approach/parse_forgi.py
   python graph_approach/calculate_angles.py
   ```
3. Run classification:
   ```bash
   python graph_approach/graphs_classifier.py
   ```
   Or use the interactive notebook `graph_approach/graph_classification.ipynb`.

## Citation

If you use this software in your research, please cite it using the metadata provided in `CITATION.cff`:

```bibtex
@software{mackowiak2026multiloops,
  title   = {multiloops-classification: Classification of RNA multiloops},
  author  = {Mackowiak, Marta and Zok, Tomasz},
  year    = {2026},
  url     = {https://github.com/tzok/multiloops-classification},
  version = {1.0.0},
  license = {MIT}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Authors

- **Marta Mackowiak** ([ORCID: 0009-0000-4777-2786](https://orcid.org/0009-0000-4777-2786))
- **Tomasz Zok** ([ORCID: 0000-0003-4103-9238](https://orcid.org/0000-0003-4103-9238))
