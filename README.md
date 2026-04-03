# High-Throughput, High-Quality: Benchmarking GNINA and AutoDock Vina for Precision Virtual Screening Workflow
R. Buccheri¹ and A. Rescifina¹\
¹*University of Catania*

[![Python](https://img.shields.io/badge/Python-3.11.7-blue.svg)](https://www.python.org/)

## Shared Files
**CSV file 1** “experimental_ligands_data.csv”.\
Complete database used in VS reporting SMILES codes of each compound (SMILES), experimental *K*ᵢ numeric value expressed in nM (*K*), the CNN_VS values (CNN_VS), the p*K*ᵢ derived from Vina output (Vina p*K*ᵢ) and the corresponding target protein (Target protein). Note that Vina failed to analyze 13 ligands reported as #NUM! values in the Vina p*K*ᵢ column.

**CSV file 2** “decoy_data.csv”.\
Complete database used in VS reporting SMILES codes of each decoy molecule (SMILES), the CNN_VS values (CNN_VS), the p*K*ᵢ derived from Vina output (Vina p*K*ᵢ), and the corresponding target protein (Target protein). Note that Vina failed to analyze 172 ligands reported as #NUM! values in the Vina p*K*ᵢ column.

**Python script** “ROC_AUC_EF.py”. Python script used in ROC-AUC and EFs cal-culations.

## Citation
If you find Zinc Package useful in your own research please cite:

**High-Throughput, High-Quality: Benchmarking GNINA and AutoDock Vina for Precision Virtual Screening Workflow**\
[https://doi.org/10.3390/molecules30163361]