# Pregunta 3/4 — Validación del modelo de AlphaFold

**Idea central:** un modelo de AlphaFold NO se valida con ProSA/PROCHECK (esas herramientas
son para modelos por homología / estructuras experimentales). Se valida con las **métricas
de confianza propias de AlphaFold**.

## Métricas de validación (nuestro modelo, rank_001)
| Métrica | Qué mide | Valor | Lectura |
|---|---|---|---|
| **pLDDT** | confianza por residuo (0–100) | global **86,8**; periplásmico 93,1; núcleo TM 87,3 | confiable donde importa; solo el conector 30–69 bajo |
| **pTM** | confianza del plegado global | **0,797** | >0,5 = topología confiable (~0,8 muy bueno) |
| **PAE** | error de posición entre residuos (Å) | intra-dominio 4,6 / 5,4; inter-dominio 8,3; max 31,3 | dominios bien definidos; orientación relativa con algo de flexibilidad |
| **Consistencia 5 modelos** | reproducibilidad | pLDDT 85–86,8; pTM 0,78–0,80 | predicción robusta |

## Por qué ProSA/PROCHECK no son la validación de AlphaFold
- ProSA (Z-score) y PROCHECK (Ramachandran) evalúan estereoquímica/energía frente a
  estructuras reales; se usan para **modelos por homología (Modeller)**.
- AlphaFold ya produce buena geometría y entrega su **propia estimación de confianza**
  (pLDDT/PAE/pTM), que es la métrica correcta a reportar.
- En el taller, ProSA/PROCHECK/QMEAN figuran en el bloque de **modelado por homología**
  (HHpred+Modeller), no para AlphaFold.

## Si de todos modos se pide ProSA/PROCHECK
- Aceptan cualquier `.pdb`: se sube el modelo. Usar el modelo **relaxed (Amber)** para
  chequeos de geometría (el unrelaxed puede tener choques menores).
- Estándar moderno de geometría: **MolProbity** (clashscore, Ramachandran, rotámeros)
  https://molprobity.biochem.duke.edu/
- **QMEANDisCo** sí es apropiado para modelos y correlaciona con la calidad.

## Texto listo para el informe
> La validez del modelo se evaluó mediante las métricas de confianza propias de AlphaFold,
> que son el criterio adecuado para modelos predichos por aprendizaje profundo (a
> diferencia de ProSA/PROCHECK, propios de la validación de modelos por homología). El
> modelo presenta un **pLDDT global de 86,8** (con valores de 93,1 en el dominio
> periplásmico y 87,3 en el núcleo transmembrana), un **pTM de 0,797** —indicativo de una
> topología global confiable— y un **PAE intra-dominio bajo** (4,6–5,4 Å), lo que confirma
> que ambos dominios están bien definidos. Los cinco modelos generados son consistentes
> (pLDDT 85–86,8; pTM 0,78–0,80), lo que respalda la robustez de la predicción.
