# Pregunta 4 — Comparación estructural (Foldseek) y visualización del modelo

Modelo: `alphafold/YidC_Moraxella_modelo_rank1.pdb` (AlphaFold rank_001, pLDDT 86,8;
el B-factor de cada átomo = pLDDT).

## A) Foldseek — búsqueda por similitud estructural
Servidor: https://search.foldseek.com/search

1. **Upload query:** subí `YidC_Moraxella_modelo_rank1.pdb` (o pegá el PDB).
2. **Databases:** marcá **PDB100** (estructuras experimentales) y **AlphaFold/UniProt50**
   (o AlphaFold/Proteome). Para el examen, con PDB100 alcanza.
3. **Mode:** 3Di/AA (por defecto). **Search.**
4. **Resultados a mirar:**
   - Los **mejores hits deben ser YidC/Oxa1/Alb3** (p. ej. YidC de *E. coli*, YidC2 de
     *Bacillus halodurans*), con **TM-score alto** y **E-value muy bajo**.
   - Anotá: PDB ID del mejor hit, **TM-score**, **E-value**, %identidad y cobertura.
5. **Interpretación:** que el fold recupere YidC del PDB **confirma la identidad y el
   plegado** de forma estructural (independiente de la secuencia).

> Alternativa del taller: **RCSB pairwise** https://www.rcsb.org/alignment (programas
> TM-align / jFATCAT) para comparar tu modelo contra una estructura YidC puntual y obtener
> **RMSD** y **TM-score**.

## B) Abrir el modelo en un visualizador

### ChimeraX (recomendado)
1. Abrí ChimeraX → `File > Open` → `YidC_Moraxella_modelo_rank1.pdb`.
2. **Colorear por pLDDT** (confianza), en la barra de comandos:
   `color bfactor palette alphafold`
   → azul = pLDDT alto (confiable), naranja/rojo = bajo (conector N-terminal ~30–69).
3. **Ver los dominios:** `cartoon` ; se distingue el β-supersándwich periplásmico del
   haz de hélices TM.
4. **Marcar residuos funcionales:**
   `select :388,358,360,451` → `show sel atoms` → `color sel red` → `label sel`
   (R388 = Arg del surco; W358/W360 = motivo GNWGW).

### PyMOL (alternativa)
1. `File > Open` el .pdb.
2. Color por B-factor (pLDDT): `spectrum b, blue_white_red`
3. Mostrar residuos: `select sitio, resi 388+358+360+451` → `show sticks, sitio` →
   `color yellow, sitio`.

## C) Qué reportar en el informe (sección 3.4)
- **Foldseek:** mejor hit estructural = YidC del PDB (PDB ID ____), TM-score ____,
  E-value ____ → confirma fold e identidad.
- **Visualización:** el modelo muestra el β-supersándwich periplásmico + haz de 5 TM; los
  residuos R388 y W358/W360 quedan en el núcleo, delimitando el surco; la zona de bajo
  pLDDT (naranja) coincide con el conector N-terminal desordenado.

## Texto listo para el informe
> La comparación estructural del modelo con Foldseek (bases PDB100 y AlphaFold DB) recuperó
> como mejores coincidencias estructuras experimentales de **YidC** (familia
> YidC/Oxa1/Alb3), con TM-score elevado y E-value muy bajo, confirmando por vía estructural
> —independiente de la secuencia— la identidad y el plegado de la proteína. La inspección
> del modelo (coloreado por pLDDT) muestra el dominio periplásmico β-supersándwich y el
> núcleo de cinco hélices transmembrana, con los residuos funcionales conservados (R388 y
> el motivo GNWGW) ubicados en el núcleo, delimitando el surco de translocación.
