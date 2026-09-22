# Pregunta 3/4 — Guía paso a paso: modelado estructural (AlphaFold)

Herramientas del taller. Secuencia: YidC de *Moraxella* sp. Pampa (542 aa).

## Paso 1 — Herramienta
- **ColabFold – AlphaFold2_mmseqs2** (la del taller):
  https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- Alternativa: **AlphaFold3 server** https://alphafoldserver.com (single protein).

## Paso 2 — Entrada
- Pegar la secuencia (sin el `>`), `jobname` = `YidC_Moraxella`. Es un **monómero**.

## Paso 3 — Configuración
- msa_mode: `mmseqs2_uniref_env`
- model_type: `auto` (alphafold2_ptm)
- num_recycles: `3` (subir a 6 si pLDDT bajo)
- templates: ON (pdb100); relax/amber: opcional
- Entorno de ejecución → Ejecutar todas.

## Paso 4 — Resultados
- `rank_1.pdb` = mejor modelo.
- **pLDDT** (0–100): >90 muy alta, 70–90 buena, 50–70 baja, <50 desordenado.
- **PAE**: orientación entre dominios.

## Paso 5 — Interpretación (informe)
- Debe mostrar el **β-supersándwich periplásmico + haz de 5 TM**.
- **pLDDT bajo (<50) en el conector N-terminal 1–69** → confirma el desorden de IUPred.

## Paso 6 — Topología sobre el modelo
- **MembraneFold** https://ku.biolib.com/MembraneFold/ (mapea DeepTMHMM en 3D).

## Paso 7 — Mapear conservación
- ChimeraX/PyMOL (o ESPript + PDB): verificar que **R388**, **GNWGW** y los ácidos
  del surco quedan en la cavidad del núcleo TM.

## Paso 8 — Validación
- **ProSA** (Z-score): https://prosa.services.came.sbg.ac.at/prosa.php
- **PROCHECK / SAVES v6** (Ramachandran): https://saves.mbi.ucla.edu/
- **QMEANDisCo**: https://swissmodel.expasy.org/qmean/

## Paso 9 — Comparación estructural (P4)
- **Foldseek** https://search.foldseek.com/search → PDB100 + AFDB50.
- **RCSB TM-align/jFATCAT** https://www.rcsb.org/alignment → contra estructuras YidC
  (E. coli, B. halodurans) → reportar RMSD / TM-score; verificar superposición de R388
  con la Arg del surco (R366 en E. coli).

## Qué guardar como evidencia
- `modelo_rank1.pdb`, captura del modelo coloreado por pLDDT, matriz PAE, salida de
  Foldseek/TM-align, y los puntajes de ProSA/PROCHECK/QMEAN.
