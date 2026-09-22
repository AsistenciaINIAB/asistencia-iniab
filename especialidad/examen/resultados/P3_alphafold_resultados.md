# Pregunta 3/4 — Resultados del modelado estructural (AlphaFold / ColabFold)

**Herramienta:** ColabFold — AlphaFold2 (alphafold2_ptm), msa_mode=mmseqs2_uniref_env,
3 recycles. **Entrada:** YidC de *Moraxella* sp. Pampa (542 aa).
**Salida:** `alphafold/YidC_Moraxella_2720e/` (5 modelos + plots pLDDT/PAE/coverage).

## Calidad del modelo
- **Mejor modelo (rank_001, model_3): pLDDT global = 86,8 · pTM = 0,797.**
- Los 5 modelos son consistentes (pLDDT 85,0–86,8; pTM 0,78–0,80) → predicción robusta.

## pLDDT por región (confianza por residuo)
| Región | Residuos | pLDDT medio | Interpretación |
|---|---|---|---|
| Dominio periplásmico | 70–362 | **93,1** | β-supersándwich bien plegado (muy alta confianza) |
| Núcleo transmembrana | 363–542 | **87,3** | haz de 5 hélices bien definido |
| Conector N-terminal | 1–69 | 58,6 | baja confianza |
| — subtramo de baja complejidad | 30–60 | **41,1** | muy baja confianza |

## Resultado clave: convergencia con IUPred
Los residuos con **pLDDT < 50** son **31–51 y 58–70**, es decir, AlphaFold asigna baja
confianza **exactamente a la misma región (~30–69)** que **IUPred** predijo como
desordenada. **Dos métodos independientes (predictor de desorden y pLDDT del modelo)
coinciden** en que la única región flexible/desordenada es el **conector N-terminal de
baja complejidad**, mientras que el resto de la proteína está bien estructurado.

## Residuos funcionales bien modelados (alta confianza)
| Residuo | pLDDT | Rol |
|---|---|---|
| **R388** | 86,2 | Arg conservada del surco ("greasy slide") |
| **W358 / W360** | 94,0 / 96,1 | motivo GNWGW (cinturón aromático) |
| R451 | 88,3 | Arg conservada (bucle periplásmico) |
| C427 | 79,3 | Cys conservada |

El surco de translocación y sus residuos clave se modelan con alta confianza, lo que
refuerza la interpretación funcional de la §3.4.

## Pasos siguientes (opcionales)
- **MembraneFold** para mapear las 6 TM sobre el modelo.
- **Foldseek / RCSB TM-align** contra estructuras YidC del PDB (E. coli, B. halodurans).
- Validación: **ProSA / PROCHECK / QMEANDisCo**.

## Texto listo para el informe (sección 3.4, complemento estructural)
> El modelo tridimensional obtenido con AlphaFold (ColabFold, pLDDT global 86,8; pTM 0,797)
> reproduce la arquitectura esperada de YidC: un **dominio periplásmico β-supersándwich**
> (pLDDT 93,1) y un **núcleo de cinco hélices transmembrana** (pLDDT 87,3). De forma
> notable, la **única región con baja confianza (pLDDT < 50) es el extremo N-terminal
> (residuos ~30–69)**, coincidiendo con la región desordenada predicha por IUPred; este
> acuerdo entre dos métodos independientes confirma que dicho conector de baja complejidad
> es la única zona flexible de la proteína. Los residuos funcionales identificados por
> conservación —la arginina del surco **R388** (pLDDT 86,2) y los triptófanos del motivo
> **GNWGW** (W358/W360, pLDDT 94–96)— se modelan con alta confianza y se ubican en el
> núcleo transmembrana, delimitando el surco de translocación del sustrato.
