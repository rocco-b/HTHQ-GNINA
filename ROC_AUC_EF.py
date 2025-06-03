'''
================================================================================
High-Throughput, High-Quality: unlocking GNINA for preci-sion virtual screening
================================================================================
Rocco Buccheri and Antonio Rescifina
    Department of Drug and Health Sciences, University of Catania, Viale A. Doria 6, 95125 Catania, Italy; roc-co.buccheri@unict.it (R.B.)
    Correspondence: antonio.rescifina@unict.it



Input data file must contain columns separated by tabs (\t).

Must have the following columns:
SMILES(optional)  CNN_VS  pKi_VINA    Activity

Activity column must contain the label 'active' for experimental ligands and the label 'decoy' for deocy molecules.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc



# ===== PARAMETERS =====
file_txt = "data.txt"       # Input data file
top_fraction = 0.01         # Enrichment factor @1%
top_fraction_2 = 0.05       # Enrichment factor @5%
top_fraction_3 = 0.10       # Enrichment factor @10%
# === PARAMETERS END ===



df = pd.read_csv(file_txt, sep='\t')
df["CNN_VS"] = pd.to_numeric(df["CNN_VS"], errors="coerce")
df["pKi_VINA"] = pd.to_numeric(df["pKi_VINA"], errors="coerce")
df = df.dropna(subset=["CNN_VS", "pKi_VINA"])
df["label"] = df["Activity"].apply(lambda x: 1 if x.lower() == "active" else 0)

def calcola_roc_auc(y_true, scores):
    fpr, tpr, _ = roc_curve(y_true, scores)
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc

fpr1, tpr1, auc1 = calcola_roc_auc(df["label"], df["CNN_VS"])
fpr2, tpr2, auc2 = calcola_roc_auc(df["label"], df["pKi_VINA"])

def enrichment_factor(df, score_col, label_col, fraction):
    df_sorted = df.sort_values(by=score_col, ascending=False).reset_index(drop=True)
    top_n = int(len(df) * fraction)
    top_n = max(top_n, 1)
    top_df = df_sorted.head(top_n)
    n_actives_total = df[label_col].sum()
    n_actives_top = top_df[label_col].sum()
    expected_random = n_actives_total * (top_n / len(df))
    ef = n_actives_top / expected_random if expected_random > 0 else 0
    return ef

ef1 = enrichment_factor(df, "CNN_VS", "label", top_fraction)
ef2 = enrichment_factor(df, "pKi_VINA", "label", top_fraction)
ef1_2 = enrichment_factor(df, "CNN_VS", "label", top_fraction_2)
ef2_2 = enrichment_factor(df, "pKi_VINA", "label", top_fraction_2)
ef1_3 = enrichment_factor(df, "CNN_VS", "label", top_fraction_3)
ef2_3 = enrichment_factor(df, "pKi_VINA", "label", top_fraction_3)

plt.figure(figsize=(7, 7))
plt.plot(fpr1, tpr1, label=f'GNINA (AUC = {auc1:.2f})', lw=2, color='blue')
plt.plot(fpr2, tpr2, label=f'Vina (AUC = {auc2:.2f})', lw=2, color='green')
plt.plot([0, 1], [0, 1], 'k--', lw=1)
plt.xlabel('False Positive Rate', fontsize=25)
plt.ylabel('True Positive Rate', fontsize=25)
plt.title('ROC Curve', fontsize=30)           # Graph title
plt.legend(loc="lower right", fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)
plt.tight_layout()
plt.show()

print(f"\nCNN_VS -> Enrichment Factor @1%: {ef1:.2f}")
print(f"pKi_VINA -> Enrichment Factor @1%: {ef2:.2f}")
print(f"\nCNN_VS -> Enrichment Factor @5%: {ef1_2:.2f}")
print(f"pKi_VINA -> Enrichment Factor @5%: {ef2_2:.2f}")
print(f"\nCNN_VS -> Enrichment Factor @10%: {ef1_3:.2f}")
print(f"pKi_VINA -> Enrichment Factor @10%: {ef2_3:.2f}")