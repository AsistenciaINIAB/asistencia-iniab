#!/usr/bin/env python3
"""Calcula conservacion por columna del MSA y la mapea a la query (seq 0)."""
import sys, os
from collections import Counter

ALN = sys.argv[1] if len(sys.argv) > 1 else "msa.aln"

def read_fasta(path):
    names, seqs, name, buf = [], [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(buf))
                name = line[1:]
                names.append(name); buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append("".join(buf))
    return names, seqs

names, seqs = read_fasta(ALN)
n = len(seqs); L = len(seqs[0])
query = seqs[0]  # Moraxella sp. Pampa (con gaps)

qpos = [0]*L
c = 0
for j in range(L):
    if query[j] != "-":
        c += 1
    qpos[j] = c

full_cons = []
high_cons = []
for j in range(L):
    col = [s[j] for s in seqs]
    non_gap = [a for a in col if a != "-"]
    if not non_gap:
        continue
    cnt = Counter(non_gap)
    aa, k = cnt.most_common(1)[0]
    frac = k/len(col)
    if len(set(col)) == 1 and "-" not in col:
        full_cons.append((qpos[j], aa))
    elif frac >= 0.90 and query[j] != "-":
        high_cons.append((qpos[j], aa, round(frac,2)))

print(f"Secuencias: {n} | Columnas: {L} | Query: {names[0].split()[0]}")
print(f"\n== Columnas 100% conservadas: {len(full_cons)} ==")
print(", ".join(f"{p}{a}" for p,a in full_cons))
print(f"\n== >=90% conservados: {len(high_cons)} ==")
print(", ".join(f"{p}{a}({f})" for p,a,f in high_cons))
print("\n== Arg 100% ==", ", ".join(f"{p}R" for p,a in full_cons if a=='R'))
print("== Trp 100% ==", ", ".join(f"{p}W" for p,a in full_cons if a=='W'))
