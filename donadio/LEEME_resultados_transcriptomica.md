# Transcriptómica E109 – Az39: resultados (versión 2)

**Cambios respecto de la versión 1:**
- **Nombres de genes en la matriz de conteos.** La matriz del proveedor (`ARG_Completo_Unique reads.xlsx`) usa el símbolo del gen (p. ej. `rplT`, `actP_1`) en lugar del locus_tag en 139 genes de E109 y 120 de Az39. Se reasignaron por posición genómica, verificando el símbolo de cada uno; no hubo ninguna inconsistencia (`scripts/00_preparar_datos.py`). Las estadísticas no cambian, pero ahora esos genes tienen anotación.
- **Hallazgo nuevo: la cepa Az39-*ipdC* no tiene el plásmido p3.** Ver la sección 2.

**Anotación:**
- E109: NCBI GCA_000807315.1 (CP010313.1).
- Az39: NCBI GCA_000632475.2 (cromosoma CP007793.1 y plásmidos p1–p5, CP007794–CP007798).

**Método:** DESeq2 1.46, FDR < 0.05 y |log2FC| ≥ 1. Ambiente: `envs/donadio-multiomica.yml`.

---

## 1. *B. japonicum* E109: ¿qué parte de la respuesta a Az39 depende del AIA?

Se usaron 12 muestras: control (1–3), + AIA (4–6), + Az39 (7–9) y + Az39-*ipdC* (10–12).

| Contraste | DE (FDR < 0.05) | up | down |
|---|---|---|---|
| AIA vs Control | 110 | 64 | 46 |
| Az39 vs Control | 2.676 | 1.270 | 1.406 |
| Az39-*ipdC* vs Control | 2.294 | 1.113 | 1.181 |
| **Az39 vs Az39-*ipdC*** | **39** | **36** | **3** |
| Az39 vs AIA | 2.623 | 1.253 | 1.370 |

**Resultado:**
- Az39 induce en E109 una **isla de unos 58 kb (RN69_30050–30290)** que contiene el cluster ***iac*** de degradación de AIA.
- 30 de los 39 genes que distinguen el co-cultivo con Az39 WT del co-cultivo con Az39-*ipdC* caen en esa isla.
- **28 de esos 39 genes también se inducen con AIA exógeno solo.**
- Con Az39-*ipdC* la isla se induce menos: log2FC alrededor de +3, frente a +8,5 con el WT.

**Cómo interpretarlo, a la luz de la sección 2:** la cepa "*ipdC*" difiere del WT en *ipdC* y también en todo el plásmido p3. Por lo tanto, la diferencia WT vs *ipdC* no se puede atribuir solo al AIA. El argumento a favor de que el AIA es la causa es que el AIA exógeno, sin Az39, reproduce la inducción de la isla, con magnitudes parecidas (log2FC de +7 a +8). La conclusión que se sostiene es esta: **la isla responde al AIA, y en co-cultivo se induce con una intensidad que depende de la cepa de Az39 usada.**

La mayoría de la respuesta de E109 al co-cultivo, unos 2.000 genes, es igual con las dos cepas de Az39 y no depende del AIA. Hay que tomarla con cautela, porque la comparación contra el control tiene sesgos de diseño:
- el control recibe medio fresco y el co-cultivo recibe medio gastado de Az39;
- E109 tiene pocas lecturas en el co-cultivo (0,5–1,5 M);
- el mapeo fue secuencial, así que algunas lecturas de Az39 pueden haber caído en genes homólogos de E109.

## 2. *A. argentinense* Az39: respuesta a E109 y efecto del mutante

Diseño 2×2 con 12 muestras:

| Genotipo | En co-cultivo con E109 | Sola (MMAB) |
|---|---|---|
| WT | 7–9 | 13–15 |
| *ipdC* | 10–12 | 16–18 |

### 2.1 La cepa *ipdC* perdió el plásmido p3 (resultado crítico)

- **El plásmido p3 (CP007796.1, ~0,69 Mb, 541 genes) no tiene ninguna lectura en las 6 muestras de la cepa *ipdC*.** Tiene 0,00% de las lecturas, con 0 a 2 genes detectados de 541, frente a 1,8–6,1% de las lecturas en las 6 muestras WT (figura `04_plasmido_p3_por_muestra.png`).
- **El transcrito de *ipdC* (ABAZ39_11295) se detecta en el mutante** en niveles parecidos al WT (CPM 95–488 en el mutante, 80–2.020 en el WT). Esto es compatible con una mutación por inserción que trunca la proteína sin eliminar todo el ARN, pero **no confirma la mutación**.
- De los 903 genes DE entre WT e *ipdC* (cultivo solo), **281 son de p3**. Quedan 622 fuera de p3, que pueden reflejar efectos de la falta de AIA o efectos secundarios de haber perdido el plásmido.

**Hay que verificar en el laboratorio:**
1. PCR de genes de p3, o perfil plasmídico (Eckhardt), en la cepa *ipdC* que se usó y en la WT.
2. Confirmar la mutación de *ipdC* por PCR o secuenciación del sitio de inserción.
3. AIA por Salkowski en el sobrenadante de la cepa *ipdC*, para saber cuánto AIA produce realmente.

Desde los FASTQ también se puede ver la cobertura a lo largo de *ipdC* para ubicar la inserción.

### 2.2 Respuesta de Az39 a E109

| Contraste | DE (FDR < 0.05, \|log2FC\| ≥ 1) | up | down |
|---|---|---|---|
| WT: co-cultivo vs sola | 309 | 214 | 95 |
| *ipdC*: co-cultivo vs sola | 69 | 53 | 16 |
| WT vs *ipdC* (sola) | 903 (281 en p3) | 600 | 303 |
| WT vs *ipdC* (co-cultivo) | 815 (300 en p3) | 590 | 225 |
| Interacción genotipo × co-cultivo | 3 | 2 | 1 |

- Con E109, Az39 WT activa:
  - genes de **desnitrificación y respiración de nitrato** (nitrato reductasa, *norD*);
  - un **transportador de C4-dicarboxilatos** (consistente con la caída de malato y fumarato en la metabolómica);
  - genes de síntesis de NAD (quinolinato sintetasa, nicotinato-nucleótido pirofosforilasa).
- Reprime transportadores ABC de aminoácidos, la degradación de histidina (N-formilglutamato amidohidrolasa) y la malonil-CoA sintasa.
- **El plásmido p2 está sobrerrepresentado** entre los genes que responden a E109: 52 genes contra 32 esperados (FDR = 0,002).
- **La expresión de p3 sube en co-cultivo**: 2,8–6,1% de las lecturas contra 1,8–2,1% en Az39 sola.
- La interacción genotipo × co-cultivo tiene poca potencia con n = 3 (solo 3 genes). Con el mutante, la respuesta a E109 es más chica (69 genes contra 309), pero eso no se puede afirmar estadísticamente.
- Las muestras 15-ARG y 18-ARG tienen muy pocas lecturas (0,27 y 0,32 M) y le restan potencia al grupo "sola".

## Archivos

- `resultados/E109/`: igual que en la versión 1 y con todos los genes anotados (`04_tabla_maestra…`, `02_clasificacion…`, `DE_*.tsv`).
- `resultados/Az39/`:
  - `01_expresion_ipdC.tsv`
  - `02_resumen_contrastes.tsv`
  - `03_tabla_maestra_log2FC_padj.tsv`
  - `04_replicones_respuesta_E109.tsv`
  - `05_fraccion_lecturas_por_replicon.tsv`
  - `DE_*.tsv`
- `figuras/E109/` y `figuras/Az39/`: PCA, correlaciones, volcanos, cluster *iac*, isla AIA y plásmido p3.
- `scripts/`: de 00 a 04, en el orden en que se corren.
- `datos_procesados/`: conteos con locus_tag corregido, anotaciones y tablas de corrección de nombres.
