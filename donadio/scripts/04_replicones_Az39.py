#!/usr/bin/env python
"""
04 - Fracción de lecturas de Az39 por replicón y por muestra.
Muestra que la cepa Az39-ipdC no tiene lecturas del plásmido p3 (CP007796.1).
Uso: python -I 04_replicones_Az39.py <directorio_de_trabajo>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

wd = Path(sys.argv[1])
cnt = pd.read_csv(wd / "data" / "counts_Az39.tsv", sep="\t", index_col=0)
ann = pd.read_csv(wd / "data" / "anotacion_Az39.tsv", sep="\t", index_col=0, quoting=3)
nombres = {"CP007793.1": "cromosoma", "CP007794.1": "p1", "CP007795.1": "p2",
           "CP007796.1": "p3", "CP007797.1": "p4", "CP007798.1": "p5"}
rep = ann.loc[cnt.index, "replicon"].map(nombres)
grupo = {**{f"{i}-ARG": "WT + E109" for i in (7, 8, 9)},
         **{f"{i}-ARG": "ipdC + E109" for i in (10, 11, 12)},
         **{f"{i}-ARG": "WT sola" for i in (13, 14, 15)},
         **{f"{i}-ARG": "ipdC sola" for i in (16, 17, 18)}}

frac = cnt.groupby(rep).sum().div(cnt.sum(), axis=1) * 100
detect = (cnt[rep == "p3"] >= 5).sum()
tabla = frac.T.round(3)
tabla.insert(0, "grupo", [grupo[s] for s in tabla.index])
tabla["genes_p3_con_>=5_lecturas"] = detect
tabla["genes_p3_total"] = int((rep == "p3").sum())
out = wd / "results" / "Az39"
out.mkdir(parents=True, exist_ok=True)
tabla.to_csv(out / "05_fraccion_lecturas_por_replicon.tsv", sep="\t")
print(tabla.to_string())

orden = [f"{i}-ARG" for i in (13, 14, 15, 7, 8, 9, 16, 17, 18, 10, 11, 12)]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), gridspec_kw={"width_ratios": [3, 2]})
colores = {"cromosoma": "#4D4D4D", "p1": "#0072B2", "p2": "#56B4E9", "p3": "#D55E00",
           "p4": "#009E73", "p5": "#E69F00"}
base = pd.Series(0.0, index=orden)
for r in ["cromosoma", "p1", "p2", "p3", "p4", "p5"]:
    ax[0].bar(range(len(orden)), frac.loc[r, orden], bottom=base, color=colores[r], label=r, width=0.8)
    base += frac.loc[r, orden]
ax[0].set_xticks(range(len(orden)))
ax[0].set_xticklabels([f"{s}\n{grupo[s]}" for s in orden], rotation=60, ha="right", fontsize=7)
ax[0].set_ylabel("% de lecturas de Az39")
ax[0].set_title("Lecturas por replicón")
ax[0].legend(fontsize=7, ncol=6, loc="upper center", bbox_to_anchor=(0.5, -0.38))
ax[1].bar(range(len(orden)), frac.loc["p3", orden],
          color=["#0072B2" if "WT" in grupo[s] else "#CC79A7" for s in orden])
ax[1].set_xticks(range(len(orden)))
ax[1].set_xticklabels([grupo[s] for s in orden], rotation=60, ha="right", fontsize=7)
ax[1].set_ylabel("% de lecturas en p3")
ax[1].set_title("Plásmido p3 (CP007796.1, 541 genes)\nausente en la cepa ipdC")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(wd / "figures" / "Az39" / f"04_plasmido_p3_por_muestra.{ext}", dpi=200)
