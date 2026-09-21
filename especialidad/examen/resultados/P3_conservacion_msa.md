# Pregunta 3/4 — Posiciones relevantes por conservación (MSA)

**Herramientas:** MAFFT L-INS-i (56 homólogos YidC de Moraxella/Psychrobacter) →
conservación por columna (`scripts/conservacion_msa.py`) + ESPript.
**Entrada:** `datos/msa.aln` (56 secuencias, 595 columnas). **Numeración:** query
*Moraxella* sp. Pampa (542 aa).

## Resumen
- **169 columnas 100% conservadas** (idénticas en las 56 secuencias) y 114 residuos ≥90%.
- La conservación se concentra en el **núcleo transmembrana (surco de translocación)** y
  en el **core del β-supersándwich periplásmico**.

## Posiciones de mayor relevancia FUNCIONAL

1. **R388 — Arg estrictamente conservada del surco hidrofílico ("greasy slide").**
   En el bucle citoplásmico posterior al primer TM del núcleo (motivo `MAKMR`, 384–388).
   Equivale a la **Arg esencial de *E. coli* YidC (R366)**, imprescindible para la función
   de insertasa (crea el microambiente hidrofílico del surco). → **residuo funcional #1.**
2. **Motivo `GNWGW` (N357-W358-G359-W360), Trp 100% conservados.** Cinturón aromático en la
   interfaz membrana–periplasma, antes del núcleo; contacta el sustrato/las hélices.
3. **Residuos polares/ácidos conservados del surco:** 333D, 404D, 419E, 446E, 461D, 466D,
   491D, y Q/N conservados (410Q, 422N, 483Q, 494Q, 523N, 524N).
4. **Cinturón aromático de Trp:** 230W, 338W, 358W, 360W, 442W, 458W, 520W.
5. **R451** (bucle periplásmico `LRHAP`, 450–454) y **R138, R206** (periplásmico).
6. **C427** conservada.

## Posiciones de mayor relevancia ESTRUCTURAL
- **β-supersándwich periplásmico:** motivo `GGAWGTP` (227–232), Gly de giro (128G, 131G,
  166G, 198G, 218G, 257G, 291G, 314G) y aromáticos (120Y, 148Y, 180Y, 237Y, 265Y, 266F,
  312Y, 175F, 242F, 245F).
- **Núcleo TM:** Gly/Pro y residuos hidrofóbicos conservados (402G, 426G, 515G; 392P, 423P,
  429P, 435P, 454P, 502P, 513P).

## Interpretación (discusión, Pregunta 4)
La conservación evolutiva concentra la señal en dos zonas coherentes con la función de
YidC: (i) el **surco de translocación del sustrato** en el núcleo TM, cuyo residuo más
relevante es la **Arg conservada R388** (homóloga de la R366 esencial de *E. coli*),
acompañada del motivo aromático **GNWGW** y de residuos polares/ácidos que forman el
microambiente hidrofílico; y (ii) el **core estructural del dominio periplásmico**
(β-supersándwich). La proteína es una insertasa **YidC**, y sus posiciones críticas
(especialmente **R388**) coinciden con las del mecanismo de la familia.

## Texto listo para el informe (secciones 3.3–3.4)
> A partir de los ~56 homólogos recuperados por BLAST se construyó un alineamiento
> múltiple (MAFFT L-INS-i) y se evaluó la conservación por columna (ESPript). Se
> identificaron **169 posiciones estrictamente conservadas**, concentradas en el núcleo
> transmembrana y en el core del dominio periplásmico. La posición funcional más relevante
> es una **arginina invariante (R388)** en la cara citoplásmica del núcleo, que corresponde
> a la arginina esencial del surco hidrofílico de YidC (R366 en *E. coli*). Se conservan
> además el motivo aromático **GNWGW** (W358/W360), un cinturón de triptófanos de interfaz
> y varios residuos polares/ácidos que revisten la cavidad de translocación, junto con el
> núcleo estructural del β-supersándwich periplásmico. En conjunto, las posiciones
> conservadas delimitan el **surco de deslizamiento del sustrato**, consistente con la
> función de insertasa de membrana de la proteína.
