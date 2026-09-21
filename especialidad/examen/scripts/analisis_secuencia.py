#!/usr/bin/env python3
"""
Análisis local de la secuencia (datos/secuencia.fasta):
  - Propiedades globales: longitud, MW, pI, GRAVY, % hidrofóbicos
  - Composición aminoacídica
  - Hidrofobicidad Kyte-Doolittle (ventana 19) -> predicción cruda de TM
  - Proxy de baja complejidad / desorden (ventana de complejidad)
Guarda una figura de hidropatía en figuras/ y un resumen en resultados/.

Requiere: biopython, numpy, matplotlib
    pip install biopython numpy matplotlib
Uso:
    python3 scripts/analisis_secuencia.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(HERE, "datos", "secuencia.fasta")
FIGDIR = os.path.join(HERE, "figuras")
RESDIR = os.path.join(HERE, "resultados")

# Escala Kyte-Doolittle
KD = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
    'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
    'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
    'Y': -1.3, 'V': 4.2,
}


def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def hydropathy(seq, window=19):
    half = window // 2
    scores = []
    for i in range(len(seq)):
        lo, hi = max(0, i - half), min(len(seq), i + half + 1)
        vals = [KD.get(a, 0.0) for a in seq[lo:hi]]
        scores.append(sum(vals) / len(vals))
    return np.array(scores)


def predict_tm(scores, thr=1.6, min_len=15):
    """Segmentos TM crudos: tramos con KD medio > thr y longitud >= min_len."""
    segs, start = [], None
    for i, s in enumerate(scores):
        if s > thr and start is None:
            start = i
        elif s <= thr and start is not None:
            if i - start >= min_len:
                segs.append((start + 1, i))
            start = None
    if start is not None and len(scores) - start >= min_len:
        segs.append((start + 1, len(scores)))
    return segs


def low_complexity(seq, window=25):
    """Proxy de baja complejidad: nº de aa distintos en la ventana (bajo = LCR)."""
    half = window // 2
    out = []
    for i in range(len(seq)):
        lo, hi = max(0, i - half), min(len(seq), i + half + 1)
        out.append(len(set(seq[lo:hi])))
    return np.array(out)


def main():
    seq = read_fasta(FASTA)
    n = len(seq)

    # Propiedades globales con Biopython si está disponible
    props = {}
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        pa = ProteinAnalysis(seq)
        props["MW_Da"] = round(pa.molecular_weight(), 1)
        props["pI"] = round(pa.isoelectric_point(), 2)
        props["GRAVY"] = round(pa.gravy(), 3)
        props["aromaticity"] = round(pa.aromaticity(), 3)
        props["instability"] = round(pa.instability_index(), 2)
    except Exception as e:
        props["biopython"] = f"no disponible: {e}"

    hyd = hydropathy(seq)
    tms = predict_tm(hyd)
    lc = low_complexity(seq)

    hydro_res = sum(1 for a in seq if a in "AVILMFWC")
    comp = {a: seq.count(a) for a in sorted(set(seq))}

    # ---- salida de texto ----
    lines = []
    lines.append("=== RESUMEN DE ANALISIS DE SECUENCIA ===")
    lines.append(f"Longitud: {n} aa")
    lines.append(f"Residuos hidrofobicos (AVILMFWC): {hydro_res} ({100*hydro_res/n:.1f}%)")
    for k, v in props.items():
        lines.append(f"{k}: {v}")
    lines.append("")
    lines.append(f"Segmentos TM predichos (KD>1.6, len>=15): {len(tms)}")
    for i, (a, b) in enumerate(tms, 1):
        lines.append(f"  TM{i}: {a}-{b}  ({b-a+1} aa)  {seq[a-1:b]}")
    lines.append("")
    lc_regions = []
    in_lc, st = False, None
    for i, v in enumerate(lc):
        if v <= 11 and not in_lc:
            in_lc, st = True, i
        elif v > 11 and in_lc:
            if i - st >= 15:
                lc_regions.append((st + 1, i))
            in_lc = False
    lines.append(f"Regiones de baja complejidad (proxy desorden): {len(lc_regions)}")
    for a, b in lc_regions:
        lines.append(f"  LCR: {a}-{b}  {seq[a-1:b]}")
    lines.append("")
    lines.append("Composicion aminoacidica:")
    for a, c in comp.items():
        lines.append(f"  {a}: {c} ({100*c/n:.1f}%)")

    os.makedirs(RESDIR, exist_ok=True)
    with open(os.path.join(RESDIR, "resumen_secuencia.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # ---- figura ----
    os.makedirs(FIGDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 4))
    x = np.arange(1, n + 1)
    ax.plot(x, hyd, lw=1, color="steelblue")
    ax.axhline(1.6, color="red", ls="--", lw=0.8, label="umbral TM (1.6)")
    ax.axhline(0, color="grey", lw=0.5)
    for a, b in tms:
        ax.axvspan(a, b, color="orange", alpha=0.25)
    ax.set_xlabel("Posicion (residuo)")
    ax.set_ylabel("Hidropatia KD (ventana 19)")
    ax.set_title("Perfil de hidrofobicidad Kyte-Doolittle")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "hidropatia_kd.png"), dpi=150)
    print(f"\nFigura -> {os.path.join(FIGDIR, 'hidropatia_kd.png')}")


if __name__ == "__main__":
    main()
