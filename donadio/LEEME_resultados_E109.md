# Transcriptómica de *B. japonicum* E109: ¿qué parte de la respuesta a Az39 depende del AIA?

**Datos:** matriz de conteos `ARG_Completo_Unique reads.xlsx`, hoja Bj_E109. Son 12 muestras de E109 cosechadas a las 4 h:

| Condición | Muestras |
|---|---|
| Control | 1–3-ARG |
| + AIA 40 µg/ml | 4–6-ARG |
| + Az39 | 7–9-ARG |
| + Az39-*ipdC* (mutante que no produce AIA) | 10–12-ARG |

**Anotación:** NCBI GCA_000807315.1 (*B. japonicum* E109, CP010313.1). Los 8.621 locus `RN69_` coinciden exactamente con la matriz de conteos.

**Método:**
- DESeq2 1.46 con un único modelo `~ condición` que incluye las 4 condiciones.
- Filtro: ≥10 lecturas en al menos 3 muestras, de modo que quedan 6.637 genes.
- Contrastes con test de Wald. log2FC con y sin shrinkage (apeglm).
- Corte: FDR < 0.05 y |log2FC| ≥ 1.
- Scripts: `scripts/01_deseq2_E109.R` y `scripts/02_isla_AIA_E109.R`. Ambiente: `envs/donadio-multiomica.yml`.

## Resultado principal

**El AIA que produce Az39 induce en E109 una isla genómica de unos 58 kb (RN69_30050–30290, 6,31–6,37 Mb) que contiene el cluster *iac* de degradación de AIA.**

- La comparación **Az39 vs Az39-*ipdC*** da solo **39 genes diferenciales**, 36 de ellos inducidos. **30 de los 39 están en esa isla.** Ahí se encuentran:
  - el cluster *iac*: dioxigenasa de fenilpropionato *iacC/iacD*, ferredoxina, oxidorreductasa;
  - la oxidorreductasa de indolpiruvato (dos subunidades);
  - la gentisato 1,2-dioxigenasa;
  - una nitrilasa;
  - transportadores ABC;
  - reguladores MarR.
- **28 de esos 39 genes también se inducen con AIA exógeno.** El AIA solo alcanza para reproducir la respuesta.
- Con el mutante *ipdC* la isla se sigue induciendo, pero mucho menos: log2FC de alrededor de +3 frente a +8,5 con el WT. Es consistente con que el mutante produce AIA residual, ya descripto para mutantes *ipdC* de *Azospirillum*. Conviene medir AIA por Salkowski en el sobrenadante del mutante para confirmarlo.
- Figuras: `03_operon_iac.png` y `05_isla_AIA_region_genomica.png`.

## Números por contraste (FDR < 0.05)

| Contraste | DE | up | down | con \|log2FC\| ≥ 1 |
|---|---|---|---|---|
| AIA vs Control | 110 | 64 | 46 | 108 |
| Az39 vs Control | 2.676 | 1.270 | 1.406 | 2.586 |
| Az39-*ipdC* vs Control | 2.294 | 1.113 | 1.181 | 2.203 |
| **Az39 vs Az39-*ipdC*** | **39** | **36** | **3** | **39** |
| Az39 vs AIA | 2.623 | 1.253 | 1.370 | 2.527 |
| Az39-*ipdC* vs AIA | 2.263 | 1.117 | 1.146 | 2.182 |

Clasificación de los 2.586 genes que responden a Az39 (`02_clasificacion_respuesta_Az39.tsv`):

| Clase | up | down |
|---|---|---|
| AIA-independiente (responde igual al mutante) | 978 | 1.021 |
| Parcialmente AIA-dependiente | 25 | 3 |
| AIA-dependiente | 11 | 0 |
| Ambiguo (sin potencia para decidir) | 223 | 325 |

## Cómo interpretarlo y qué falta controlar

1. **La gran mayoría de la respuesta de E109 al co-cultivo (unos 2.000 genes) no depende del AIA.** El co-cultivo separa las muestras en el PC1, con el 67% de la varianza. Una parte puede ser biológica, por otras señales o metabolitos de Az39. Otra parte puede ser técnica:
   - el co-cultivo recibe medio agotado de Az39 (MMAB con Trp), mientras que el control recibe medio fresco;
   - E109 aporta pocas lecturas en el co-cultivo (0,5–1,5 M) porque Az39 domina la librería;
   - algunas lecturas de Az39 pueden haberse asignado a genes homólogos de E109, porque el mapeo fue secuencial.

   **El contraste Az39 vs Az39-*ipdC* no tiene ninguno de estos problemas**, porque los dos tratamientos comparten medio y proporción. Por eso es la evidencia más sólida.
2. El AIA exógeno cambia pocos genes (110). La muestra 4-ARG queda cerca de los controles.
3. Los genes "AIA-independientes" son el punto de partida para buscar las otras señales del co-cultivo. Se cruzan con la metabolómica (Trp, malato/fumarato, osmoprotectores) en el paso siguiente.
4. Los log2FC calculados con el modelo de 12 muestras coinciden exactamente con el DESeq2 del proveedor (Exp.3: r = 1,000), lo que valida el procesamiento.

## Archivos

- `resultados/E109/00_QC_librerias.tsv`: lecturas asignadas por muestra.
- `resultados/E109/DE_<contraste>.tsv`: resultados completos con anotación.
- `resultados/E109/02_clasificacion_respuesta_Az39.tsv`: clase de cada gen que responde a Az39.
- `resultados/E109/04_tabla_maestra_log2FC_padj.tsv`: todos los contrastes en una sola tabla.
- `resultados/E109/05_operon_iac.tsv`
- `figuras/E109/`: PCA, correlación entre muestras, cluster *iac*, volcanos y región genómica.
