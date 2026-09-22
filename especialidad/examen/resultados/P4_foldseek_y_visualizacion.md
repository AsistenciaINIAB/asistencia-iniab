# Pregunta 4 — Comparación estructural (Foldseek) y visualización del modelo

Modelo: `alphafold/YidC_Moraxella_modelo_rank1.pdb` (AlphaFold rank_001, pLDDT 86,8;
el B-factor de cada átomo = pLDDT).

## A) Foldseek — búsqueda por similitud estructural
Servidor: https://search.foldseek.com/search

1. **Upload query:** subí `YidC_Moraxella_modelo_rank1.pdb` (o pegá el PDB).
2. **Databases:** AFDB50 / AFDB-SwissProt / CATH50 / PDB100 (Foldseek busca en varias).
3. **Mode:** 3Di/AA (por defecto). **Search.**
4. **Resultados a mirar:** los mejores hits deben ser **YidC**, con **probabilidad alta**
   y **E-value muy bajo**; anotá target, organismo, probabilidad, E-value y cobertura.
5. **Interpretación:** que el fold recupere YidC **confirma identidad y plegado** de forma
   estructural (independiente de la secuencia).

> Alternativa del taller: **RCSB pairwise** https://www.rcsb.org/alignment (TM-align /
> jFATCAT) para comparar contra una estructura YidC puntual y obtener **RMSD** y **TM-score**.

## B) Abrir el modelo en un visualizador

### ChimeraX (recomendado)
1. Abrí ChimeraX → `File > Open` → `YidC_Moraxella_modelo_rank1.pdb`.
2. **Colorear por pLDDT** (barra de comandos): `color bfactor palette alphafold`
   → azul = pLDDT alto, naranja/rojo = bajo (conector N-terminal ~30–69).
3. **Dominios:** `cartoon` (se distingue el β-supersándwich del haz de hélices TM).
4. **Residuos funcionales:** `select :388,358,360,451` → `show sel atoms` →
   `color sel red` → `label sel` (R388 = Arg del surco; W358/W360 = GNWGW).

### PyMOL (alternativa)
1. Abrir el .pdb. 2. `spectrum b, blue_white_red` (color por pLDDT).
3. `select sitio, resi 388+358+360+451` → `show sticks, sitio` → `color yellow, sitio`.

## RESULTADOS OBTENIDOS (Foldseek)

Los mejores hits estructurales son todos **"Membrane protein insertase YidC"**
(probabilidad 1.00), con el alineamiento cubriendo casi toda la proteína (query 1–540/542):

| # | Target | Organismo | Prob. | Seq.Id | E-value | Score |
|---|---|---|---|---|---|---|
| 1 | AF-A0A2J9DZ16 | *Psychrobacter* sp. FDAARGOS_221 | 1.00 | 60,2 | **1,22e-72** | 2832 |
| 2 | A0AAD1QMW2 | *Acinetobacter* phage MD-2021a | 1.00 | 52,7 | 1,98e-68 | 2710 |
| 3 | AF-N9PZH5 | *Acinetobacter vivianii* | 1.00 | 52,5 | 1,87e-66 | 2635 |

→ La comparación **estructural** (independiente de la secuencia) recupera **YidC** con
probabilidad 1.00 y E-values ~1e-72, confirmando **identidad y plegado**. El mejor hit es
una YidC de *Psychrobacter* (familia Moraxellaceae), coherente con todo el análisis.

## C) Qué reportar en el informe (sección 3.4)
- **Foldseek:** mejor hit estructural = YidC de *Psychrobacter* sp. FDAARGOS_221
  (AF-A0A2J9DZ16), probabilidad 1.00, E-value 1,22e-72 → confirma fold e identidad.
- **Visualización:** el modelo muestra el β-supersándwich periplásmico + haz de 5 TM; los
  residuos R388 y W358/W360 quedan en el núcleo, delimitando el surco; la zona de bajo
  pLDDT (naranja) coincide con el conector N-terminal desordenado.

## Texto listo para el informe
> La comparación estructural del modelo con Foldseek (bases AlphaFold DB / CATH) recuperó
> como mejores coincidencias, con **probabilidad 1,00 y E-values del orden de 1e-72**,
> estructuras de **insertasa de membrana YidC** (la más significativa, la YidC de
> *Psychrobacter* sp. FDAARGOS_221, AF-A0A2J9DZ16), confirmando por vía estructural
> —independiente de la secuencia— la identidad y el plegado de la proteína. La inspección
> del modelo (coloreado por pLDDT) muestra el dominio periplásmico β-supersándwich y el
> núcleo de cinco hélices transmembrana, con los residuos funcionales conservados (R388 y
> el motivo GNWGW) ubicados en el núcleo, delimitando el surco de translocación.
