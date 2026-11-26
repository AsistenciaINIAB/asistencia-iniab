# Análisis multi-ómico E109 – Az39: transcriptómica, metabolómica e integración (versión 3)

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

---

## 3. Metabolómica de sobrenadantes (RMN ¹H, 4 h)

**Datos:** `Matriz_asignada_sinQC.csv` del servicio PLABEM, con concentraciones relativas por MCR-ALS:
- 140 señales; 54 asignadas a 26 metabolitos y 86 sin identificar (U_xx);
- 5 grupos × 6 réplicas: E109, E109+AIA, E109+Az39 (EA), Az39 MMAB (AzM), Az39 LB (AzLB).

**Método:**
- log2 de cada señal, modelo `~ grupo + tanda`, donde la tanda es réplicas 1–3 contra 4–6; FDR por BH. Script `05_metabolomica_SN.R`.
- **Los extractos intracelulares no se analizaron.** Su matriz está solo dentro de `ANEXO_Extractos.rar`, que supera el límite de descarga del conector. Si se descomprime en el Drive, se agrega.

**Resultados:**
- **PCA (`01_PCA_sobrenadantes.png`):**
  - los grupos se separan bien;
  - **E1–E3 forman un grupo aparte, mientras que E4–E6 quedan junto a E109+AIA**;
  - el co-cultivo (EA) cae entre E109 y Az39 MMAB, como se espera de mezclar dos cultivos.
- Solo 14 de 139 señales tienen efecto de tanda significativo. Las señales raras de E1–E3 (U_10, U_15, U_16, U_17, U_39, asparagina) aparecen únicamente en esas 3 muestras, que son **las réplicas apareadas con el RNA**. Hay que preguntarle a Evelyn qué pasó con esos frascos.
- **Dos artefactos en E109+AIA:**
  - etanol ×180: casi seguro es el solvente del AIA;
  - la señal "Triptófano3" ×150: probablemente sea el propio AIA (anillo indol), no triptófano. La carpeta `ANEXO_id_AIA_SN` sugiere que el servicio intentó identificarlo; conviene confirmarlo con un estándar.
- **Co-cultivo frente a E109 solo:** los cambios principales se explican por la composición del cultivo de Az39 que se agrega:
  - **fumarato −7,5 log2**: prácticamente 0 en todo lo que tiene Az39 (EA, AzM, AzLB);
  - **malato −1,2**: EA queda igual a AzM, 0,12 contra 0,29 en E109;
  - **triptófano +2 a +3**: EA alcanza el nivel de AzM; el MMAB de Az39 lleva Trp agregado.

  Es decir, **Az39 consume los C4-dicarboxilatos y aporta triptófano**. Las diferencias de los informes 2025/2026 en fumarato, malato y Trp no indican "activación de rutas" en E109: vienen de Az39 y de su medio.
- **Separación del co-cultivo:** si no hubiera interacción, EA = E + AzM − medio fresco. Como no hay espectros de medio fresco, solo se puede afirmar con seguridad producción emergente (EA > E + AzM). **Solo 1 señal sin identificar (U_84_3) cumple ese criterio.** Para el resto hacen falta **blancos de medio** (MMAB+Trp, PSY y AB sin inocular).

## 4. Integración RNA ↔ metabolitos

Script `06_integracion_RNA_metabolitos.R`.

**A) Por vías, comparando contrastes de condición:**

| Vía | E109 (Az39 vs control): DE / genes de la vía | Az39 (co-cultivo vs sola): DE / genes | Metabolitos |
|---|---|---|---|
| Dicarboxilatos C4 | 19 / 49 (sin enriquecimiento) | 4 / 44 | fumarato ↓↓, malato ↓ |
| Trp / AIA (incluye la isla *iac*) | **45 / 56 (FDR = 2·10⁻⁹)** | 2 / 11 | Trp ↑ (aportado por Az39) |
| Osmoprotectores | 8 / 13 | 0 / 10 | trehalosa, betaína y prolina ↑ leve (+0,15 a +1,1) |
| Glu / Asp / Asn | 9 / 22 | 1 / 17 | Glu ↑, Asp ↓, Asn ↓ |

Lectura:
- La única vía que responde de forma coherente en las dos capas es **Trp → AIA → catabolismo en E109**: Az39 aporta Trp y AIA, y E109 enciende la isla *iac*.
- **El consumo de C4-dicarboxilatos lo hace Az39.** En Az39 sube un transportador de C4-dicarboxilatos (ABAZ39_00885, +3,9). En E109 la vía no está enriquecida, así que la caída de malato y fumarato **no** es una respuesta de E109.
- Los cambios de osmoprotectores son chicos y sin respaldo transcripcional en ninguna de las dos bacterias.

**B) Correlaciones gen–metabolito en las 9 muestras apareadas de E109**

Muestras: control 1–3, AIA 4–6 y Az39 7–9, cada una con su frasco de metabolómica. Hay 839 de 7.560 pares con FDR < 0,05. **Con 3 grupos × 3 muestras, la correlación refleja casi solo qué separa al control de los tratados**. Además, el control apareado (E1–E3) es justamente el anómalo. Estas correlaciones no se pueden leer como relaciones mecánicas gen–metabolito.

**C) DIABLO, lado E109 (n = 9):**
- La correlación entre las componentes de RNA y de metabolitos es 0,98.
- La firma de RNA es **la isla *iac***: oxidorreductasa de indolpiruvato, *iacC/iacD*, ferredoxina, gentisato dioxigenasa.
- La firma de metabolitos está dominada por las **señales anómalas de E1–E3**.

Sirve como confirmación visual de que las dos ómicas separan al control de los tratamientos. No aporta una firma metabólica biológica confiable.

## 5. Conclusión integrada y pendientes

1. **El resultado sólido es que el AIA de Az39 activa en E109 la isla de catabolismo de AIA (*iac*).** Se ve en RNA (AIA exógeno ≈ co-cultivo, menor con el mutante) y es coherente con el triptófano y el AIA presentes en el sobrenadante.
2. En el sobrenadante del co-cultivo domina lo que aporta y consume Az39 (Trp ↑, C4-dicarboxilatos ↓). Con estos datos **no hay evidencia de metabolitos producidos por la interacción**, salvo una señal sin identificar.
3. **Pendientes de laboratorio, en orden de prioridad:**
   1. confirmar la pérdida del plásmido p3 en Az39-*ipdC* y la mutación de *ipdC*;
   2. revisar qué pasó con los frascos E1–E3;
   3. medir blancos de medio por RMN (MMAB+Trp, PSY, AB) e identificar la señal "Triptófano3" con un estándar de AIA;
   4. descomprimir `ANEXO_Extractos.rar` en el Drive para analizar los extractos intracelulares.
