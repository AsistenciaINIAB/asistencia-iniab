# Pregunta 2 — ¿Proteína de membrana? (DeepTMHMM)

**Herramienta:** DeepTMHMM (https://dtu.biolib.com/DeepTMHMM). **Entrada:** `secuencia.fasta`
(**542 aa**). **Figura:** `plot.png` (Most Likely Topology + Posterior Probabilities).

## Resultado
- **Tipo predicho: `alpha TM`** → proteína **integral de membrana politópica** (α-helicoidal). **Sí es de membrana.**
- **6 segmentos transmembrana**, topología **N-in / C-in**.
- **Dominio periplásmico grande (outside) 25–362 (~338 aa)** entre TM1 y el núcleo.
- Posterior probabilities ≈ 0.99 en casi toda la secuencia → predicción de alta confianza.

## Tabla de topología (de `TMRs.gff3`)

| Región | Posición | Localización | Rol |
|---|---|---|---|
| N-term | 1–3 | inside (citoplasma) | — |
| TM1 | 4–24 | membrana | ancla N-terminal |
| P1 | 25–362 | outside (periplasma) | dominio periplásmico (β-supersándwich) |
| TM2 | 363–376 | membrana | núcleo (H1) |
| bucle | 377–427 | inside | bucle citoplásmico grande |
| TM3 | 428–442 | membrana | núcleo (H2) |
| bucle | 443–469 | outside | — |
| TM4 | 470–485 | membrana | núcleo (H3) |
| bucle | 486–496 | inside | — |
| TM5 | 497–512 | membrana | núcleo (H4) |
| bucle | 513–516 | outside | — |
| TM6 | 517–530 | membrana | núcleo (H5) |
| C-term | 531–542 | inside | cola citoplásmica |

## Interpretación
La arquitectura **TM1 de anclaje + dominio periplásmico P1 + núcleo de 5 TM (TM2–TM6),
con N y C citoplásmicos**, es la topología canónica de la **insertasa YidC de bacterias
Gram-negativas**. El motivo aromático **`GNWGW`** (~356–360) queda en la interfaz
membrana–periplasma, antes de TM2 (cinturón aromático). Esto es coherente con el análisis
de identidad (YidC) y encuadra el resto del examen: el núcleo TM2–TM6 forma el surco de
translocación del sustrato.

## Texto listo para el informe (sección 3.2.a)
> El análisis con DeepTMHMM clasifica a la proteína como α-TM (proteína integral de
> membrana politópica), con **6 hélices transmembrana** y topología N-in/C-in. Presenta
> una hélice de anclaje N-terminal (TM1, 4–24), un **gran dominio periplásmico** (25–362)
> y un **núcleo de 5 TM** (TM2–TM6, 363–530), con extremos N y C citoplásmicos. Las
> probabilidades a posteriori (~0.99) indican alta confianza. Esta organización coincide
> con la topología característica de la insertasa de membrana **YidC** de bacterias
> Gram-negativas.
