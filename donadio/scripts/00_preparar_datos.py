#!/usr/bin/env python
"""
00 - Prepara las tablas de entrada del análisis transcriptómico.

1. Lee la matriz de conteos del proveedor (ARG_Completo_Unique reads.xlsx, hojas
   Bj_E109 y Az39) y la escribe como TSV.
2. Arma la tabla de anotación desde las feature tables de NCBI:
     E109 -> GCA_000807315.1_ASM80731v1   (locus RN69_)
     Az39 -> GCA_000632475.2_ASM63247v2   (locus ABAZ39_)
3. CORRIGE LOS NOMBRES DE FILA: la matriz del proveedor usa el símbolo del gen
   (p. ej. 'rplT', 'actP_1') en lugar del locus_tag para 120 genes de Az39 y 139 de
   E109. Como la matriz está en orden genómico, cada fila se asigna al locus_tag que
   ocupa esa posición, verificando que el símbolo (sin el sufijo _N) coincida con
   el de la anotación. Si alguno no coincide, el script se detiene.

Uso: python -I 00_preparar_datos.py <xlsx_conteos> <dir_feature_tables> <directorio_de_trabajo>
"""
import re
import sys
from pathlib import Path

import pandas as pd

xlsx, ref_dir, wd = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
(wd / "data").mkdir(parents=True, exist_ok=True)

ORGS = {
    "E109": ("Bj_E109", "GCA_000807315.1_ASM80731v1", "RN69_"),
    "Az39": ("Az39", "GCA_000632475.2_ASM63247v2", "ABAZ39_"),
}

for org, (sheet, acc, prefix) in ORGS.items():
    # --- anotación ---
    ft = pd.read_csv(ref_dir / f"{acc}_feature_table.txt.gz", sep="\t", dtype=str)
    ft.columns = [c.lstrip("# ") for c in ft.columns]
    genes = ft[ft["feature"] == "gene"][["locus_tag", "genomic_accession", "start", "end",
                                         "strand", "symbol", "class"]]
    prod = (ft[ft["feature"].isin(["CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "misc_RNA"])]
            [["locus_tag", "feature", "product_accession", "name", "attributes"]]
            .rename(columns={"name": "product", "feature": "tipo"})
            .drop_duplicates("locus_tag"))
    ann = (genes.merge(prod, on="locus_tag", how="left")
           .rename(columns={"genomic_accession": "replicon", "symbol": "gene_symbol"}))
    ann.to_csv(wd / "data" / f"anotacion_{org}.tsv", sep="\t", index=False)

    # --- conteos ---
    cnt = pd.read_excel(xlsx, sheet_name=sheet, index_col=0)
    cnt = cnt.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    names = list(cnt.index)
    order = list(ann["locus_tag"])
    if len(names) != len(order):
        sys.exit(f"{org}: la matriz tiene {len(names)} filas y la anotación {len(order)} genes")
    sym = dict(zip(ann["locus_tag"], ann["gene_symbol"].fillna("")))
    fixed = []
    for i, (n, lt) in enumerate(zip(names, order)):
        if str(n).startswith(prefix):
            if n != lt:
                sys.exit(f"{org}: fila {i} '{n}' no está en orden genómico (esperado {lt})")
        else:
            if re.sub(r"_\d+$", "", str(n)) != sym[lt]:
                sys.exit(f"{org}: fila {i} '{n}' no coincide con {lt} ({sym[lt]})")
            fixed.append((n, lt))
            names[i] = lt
    cnt.index = pd.Index(names, name="locus_tag")
    cnt.to_csv(wd / "data" / f"counts_{org}.tsv", sep="\t")
    pd.DataFrame(fixed, columns=["nombre_en_matriz", "locus_tag"]).to_csv(
        wd / "data" / f"correccion_nombres_{org}.tsv", sep="\t", index=False)
    print(f"{org}: {cnt.shape[0]} genes x {cnt.shape[1]} muestras; {len(fixed)} nombres corregidos")
