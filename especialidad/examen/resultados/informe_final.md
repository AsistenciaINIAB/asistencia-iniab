# Análisis bioinformático de una secuencia proteica

**Alumno:** ____________________  **Materia:** Especialidad — Taller de Bioinformática
**Fecha:** ____ / ____ / ______

---

## 1. Introducción y objetivo

Se analiza una secuencia proteica de **542 aminoácidos** de origen bacteriano con un
enfoque combinado de bioinformática de secuencia y estructura. El objetivo es (i)
identificar la proteína y su organismo de origen, (ii) caracterizarla en cuanto a su
condición de proteína de membrana, sus regiones desordenadas y sus dominios, y (iii)
definir y aplicar una estrategia para localizar posiciones de la secuencia con posible
relevancia funcional y/o estructural.

---

## 2. Materiales y métodos

| Herramienta | Uso | Parámetros principales |
|---|---|---|
| NCBI blastp (web) vs `nr` | Identidad / organismo | BLOSUM62, E=0.05, ~100 hits |
| DeepTMHMM | Topología de membrana | por defecto |
| TMHMM 2.0 / Phobius / SignalP | Topología / señal (vía InterProScan) | por defecto |
| Kyte-Doolittle | Hidrofobicidad | ventana 19 |
| AIUPred (IUPred) | Desorden intrínseco | umbral 0.5 |
| CD-search (CDD) | Dominios | CDD v3.21, E=0.01 |
| InterProScan 5.78 | Dominios / familias | Pfam, TIGRFAM, PANTHER, PRINTS, Gene3D |
| MAFFT L-INS-i | Alineamiento múltiple | `--localpair --maxiterate 1000` (56 homólogos) |
| ESPript + script propio | Conservación por columna | umbral de identidad por posición |

Secuencia de entrada: `secuencia.fasta` (542 aa).

---

## 3. Resultados y discusión

### 3.1 Identificación de la proteína y organismo (Pregunta 1)

La búsqueda por similitud (BLASTp contra `nr`) devolvió como mejor resultado la entrada
**`WP_323843415.1`** con **99,8 % de identidad, 100 % de cobertura y E-value = 0,0**
(bitscore 1116): la secuencia analizada es, por lo tanto, **esencialmente idéntica** a esa
proteína (difieren en un solo residuo). Los ~100 hits siguientes, con identidades entre
54 % y 99 % y E-value = 0,0, corresponden a la misma familia proteica, ampliamente
conservada en bacterias.

La proteína se identifica así como la **insertasa de membrana YidC** (familia
YidC/Oxa1/Alb3), y el organismo de origen es ***Moraxella* sp. Pampa** (NCBI Taxonomy ID
**3111978**), una **γ-proteobacteria Gram-negativa** de la familia Moraxellaceae. YidC es
una proteína esencial que **inserta y pliega proteínas en la membrana plasmática**, de
forma autónoma (insertasa Sec-independiente) o asociada al translocón SecYEG.

### 3.2 Caracterización de la proteína (Pregunta 2)

**3.2.a ¿Es una proteína de membrana? Topología.**
**DeepTMHMM** clasifica la proteína como **α-TM**, es decir, **integral de membrana
politópica** α-helicoidal, y predice **6 segmentos transmembrana (TM)** con topología
**N-in / C-in** (probabilidades a posteriori ≈ 0,99). El resultado fue corroborado por
**TMHMM 2.0** y **Phobius** (InterProScan) y por el perfil de hidrofobicidad de
**Kyte-Doolittle**. La organización es la típica de la YidC de bacterias Gram-negativas:

| Región | Residuos | Localización | Rol |
|---|---|---|---|
| N-terminal | 1–3 | citoplasma | — |
| TM1 | 4–24 | membrana | hélice de anclaje |
| Dominio P1 | 25–362 | periplasma | dominio periplásmico (β-supersándwich) |
| TM2–TM6 | 363–530 | membrana | núcleo de 5 hélices (surco de translocación) |
| C-terminal | 531–542 | citoplasma | cola citoplásmica |

**SignalP** detecta la hélice N-terminal (1–28) como una **señal-ancla**, no como un
péptido señal escindible. Se concluye que es una **proteína integral de membrana de 6 TM**
con un gran dominio periplásmico entre TM1 y el núcleo.

**3.2.b Regiones desordenadas.**
Con **AIUPred (IUPred)** (umbral 0,5), el **desorden global es bajo**, como corresponde a
una proteína de membrana. Se detecta **una única región desordenada en el extremo
N-terminal (residuos 1–69)**, con máximo (score > 0,90) en **~35–60**, un tramo de **baja
complejidad rico en Thr/Ala/Ser/Asp/Asn** que actúa como **conector flexible** entre el
ancla TM1 y el dominio periplásmico plegado. El resto de la proteína es ordenada.

**3.2.c Dominios presentes.**
**CD-search** e **InterProScan** coinciden: la proteína pertenece a la familia **insertasa
de membrana YidC** (gen *yidC*; modelo **PRK01318**, residuos 3–540, E-value = 0). Se
resuelven **dos dominios estructurales**:

| Dominio | Firma (Pfam) | Residuos | Descripción |
|---|---|---|---|
| Dominio periplásmico | **PF14849** (YidC_periplas) | 74–346 | β-supersándwich (CATH 2.70.98.90) |
| Núcleo de membrana | **PF02096** (60KD_IMP) | 351–539 | haz de 5 hélices transmembrana |

Corresponden a InterPro **IPR028053** e **IPR028055**, dentro de **IPR001708**
(YidC/ALB3/OXA1/COX18), con los términos GO:0032977 (actividad insertasa), GO:0051205
(inserción de proteínas en membrana) y GO:0016020 (membrana), y la familia de transporte
TCDB 2.A.9. No hay dominios catalíticos: el sitio funcional es el surco hidrofílico del
núcleo transmembrana.

### 3.3 Estrategia para identificar posiciones con relevancia funcional/estructural (Pregunta 3)

La estrategia combina **conservación evolutiva** y **contexto estructural**:

1. **Recolección de homólogos** por BLASTp/PSI-BLAST y construcción de un conjunto diverso
   (varios géneros de la misma familia).
2. **Alineamiento múltiple** con MAFFT (modalidad L-INS-i, la más precisa para conjuntos
   moderados de secuencias).
3. **Cálculo de conservación por columna** (identidad por posición; visualización con
   ESPript) para detectar los residuos **invariantes**, candidatos a relevancia funcional
   o estructural.
4. **Interpretación en el contexto de la topología y los dominios** (§3.2): las columnas
   conservadas en el núcleo transmembrana señalan el **surco de translocación**; las del
   dominio periplásmico, el **core estructural** del plegado.
5. **Complemento estructural (opcional):** modelado 3D con AlphaFold/ColabFold y comparación
   con estructuras YidC del PDB (Foldseek/TM-align) para mapear la conservación sobre el
   modelo y confirmar que los residuos conservados delimitan la cavidad funcional.

### 3.4 Resultados del análisis y discusión (Pregunta 4)

Se recuperaron **56 homólogos** de YidC (géneros *Moraxella* y *Psychrobacter*,
Moraxellaceae), se alinearon con **MAFFT L-INS-i** y se evaluó la conservación por columna.
Se identificaron **169 posiciones estrictamente conservadas** (idénticas en las 56
secuencias), concentradas en el **núcleo transmembrana** y en el **core del dominio
periplásmico**.

La posición de mayor relevancia funcional es una **arginina invariante, R388**, situada en
la cara citoplásmica del núcleo (motivo `MAKMR`, 384–388). Corresponde a la **arginina
esencial del surco hidrofílico ("greasy slide") de YidC** (equivalente a R366 en
*Escherichia coli*), imprescindible para la actividad de insertasa porque genera el
microambiente hidrofílico por el que transita el segmento del sustrato. Junto a ella se
conservan de forma estricta:

- el **motivo aromático `GNWGW`** (W358 y W360, 100 % conservados), un cinturón aromático
  en la interfaz membrana–periplasma que contacta el sustrato;
- **residuos polares/ácidos que revisten el surco** (p. ej. D404, E419, E446, D461, D466,
  D491; Q/N conservadas);
- un **cinturón de triptófanos** de interfaz (W230, W338, W442, W458, W520);
- el **núcleo estructural del β-supersándwich periplásmico** (motivo `GGAWGTP` 227–232,
  glicinas de giro y aromáticos de empaquetamiento).

**Discusión integradora.** Todas las evidencias convergen: la identificación por
similitud (YidC de *Moraxella* sp. Pampa), la topología (6 TM con dominio periplásmico) y
los dominios (PF14849 + PF02096) definen a la proteína como una **insertasa de membrana
YidC**; y el análisis de conservación localiza sus posiciones críticas en el **surco de
deslizamiento del sustrato** del núcleo transmembrana, cuyo residuo clave es la **Arg
conservada R388**. Este resultado, obtenido a partir de datos reales del alineamiento,
coincide con el mecanismo descrito para la familia y responde de manera directa la
pregunta sobre las posiciones funcional y estructuralmente relevantes.

---

## 4. Conclusiones

- **Identidad:** insertasa de membrana **YidC** (familia YidC/Oxa1/Alb3).
- **Organismo:** *Moraxella* sp. Pampa (NCBI:txid 3111978, Moraxellaceae, γ-proteobacteria
  Gram-negativa).
- **Membrana:** proteína integral politópica de **6 TM**, topología N-in/C-in.
- **Desorden:** una única región desordenada (conector N-terminal de baja complejidad,
  residuos 1–69).
- **Dominios:** periplásmico **PF14849** (74–346) + núcleo de membrana **PF02096** (351–539).
- **Posiciones relevantes:** surco de translocación del núcleo TM; residuo funcional clave
  **R388** (Arg conservada del "greasy slide"), motivo **GNWGW** y residuos polares/ácidos
  del surco; core del β-supersándwich como soporte estructural.

---

## 5. Referencias

- BLAST: Altschul et al., *J. Mol. Biol.* 1990.
- DeepTMHMM: Hallgren et al., 2022 (bioRxiv). TMHMM 2.0; Phobius (Käll et al.). SignalP-6.
- IUPred/AIUPred: Erdős, Dosztányi et al.
- InterProScan/InterPro; Pfam; CDD/CD-search; TIGRFAM; PANTHER; Gene3D; TCDB.
- MAFFT: Katoh & Standley, *Mol. Biol. Evol.* 2013. ESPript: Robert & Gouet, 2014.
- Kumazaki et al., *Sci. Rep.* 2014 (estructura de YidC de *E. coli*); Kuhn & Kiefer,
  *Mol. Microbiol.* 2017 (revisión de YidC).

---

## 6. Anexos

- `datos/secuencia.fasta`, `datos/msa.aln`
- Salidas: `B2YNZVH0016-Alignment-HitTable.csv` (BLAST), DeepTMHMM (`TMRs.gff3`, `plot.png`),
  IUPred, CD-search, InterProScan (`.tsv/.json/.gff`), ESPript (`esp.ps`).
- Scripts: `scripts/analisis_secuencia.py`, `scripts/conservacion_msa.py`,
  `scripts/examen_colab.ipynb`.
- Análisis por pregunta: `resultados/P1_blast.md`, `P2_membrana_deeptmhmm.md`,
  `P2_desorden_iupred.md`, `P2_dominios_cdsearch.md`, `P2_dominios_interpro.md`,
  `P3_conservacion_msa.md`.
