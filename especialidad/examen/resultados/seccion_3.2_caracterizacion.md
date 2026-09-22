## 3.2 Caracterización de la proteína (Pregunta 2)

### 3.2.a ¿Es una proteína de membrana? Topología

La secuencia (542 aa) se analizó con predictores de topología de membrana. **DeepTMHMM**
la clasifica como **α-TM**, es decir, una **proteína integral de membrana politópica**
α-helicoidal, y predice **6 segmentos transmembrana (TM)** con una topología **N-in / C-in**
(extremos N y C citoplásmicos). Las probabilidades a posteriori son ≈ 0,99 a lo largo de
casi toda la secuencia, lo que indica una predicción de alta confianza. El resultado fue
corroborado por **TMHMM 2.0** y **Phobius** (ejecutados dentro de InterProScan) y por el
perfil de hidrofobicidad de **Kyte-Doolittle**.

La organización predicha es la característica de la insertasa YidC de bacterias
Gram-negativas:

| Región | Residuos | Localización | Rol |
|---|---|---|---|
| N-terminal | 1–3 | citoplasma (in) | — |
| TM1 | 4–24 | membrana | hélice de anclaje N-terminal |
| Dominio P1 | 25–362 | periplasma (out) | dominio periplásmico (β-supersándwich) |
| TM2–TM6 | 363–530 | membrana | núcleo de 5 hélices (surco de translocación) |
| C-terminal | 531–542 | citoplasma (in) | cola citoplásmica |

Además, **SignalP** detecta la hélice N-terminal (1–28) como una **señal-ancla** y no como
un péptido señal escindible, coherente con que el TM1 es un ancla de membrana permanente.
En conjunto, se concluye que la proteína es una **proteína integral de membrana politópica
de 6 TM**, con un gran dominio periplásmico entre TM1 y el núcleo transmembrana.

### 3.2.b Regiones desordenadas

El desorden intrínseco se evaluó con **AIUPred (IUPred)**, usando el umbral estándar de
0,5. El **contenido de desorden global es bajo**, como es esperable para una proteína
integral de membrana (los segmentos transmembrana y el dominio periplásmico plegado son
estructuralmente ordenados). Se detecta **una única región desordenada, en el extremo
N-terminal (residuos 1–69)**, con un máximo de desorden (score > 0,90) en el segmento
**~35–60**, que corresponde a un tramo de **baja complejidad rico en Thr/Ala/Ser/Asp/Asn**.
Esta región actúa como un **conector flexible** entre la hélice de anclaje (TM1) y el
cuerpo plegado del dominio periplásmico. El resto de la proteína (dominio periplásmico y
núcleo transmembrana) presenta valores por debajo del umbral, es decir, es ordenada.

### 3.2.c Dominios presentes

La composición de dominios se determinó con **CD-search (CDD)** e **InterProScan**, con
resultados concordantes. La proteína se asigna con altísima significancia a la familia de
la **insertasa de membrana YidC** (gen *yidC*): el modelo bacteriano **PRK01318 / NF002352**
cubre casi toda la secuencia (residuos 3–540) con **E-value = 0**. Dentro de esta familia
se resuelven **dos dominios estructurales**:

| Dominio | Firma (Pfam) | Residuos | Descripción |
|---|---|---|---|
| Dominio periplásmico | **PF14849** (YidC_periplas) | 74–346 | β-supersándwich periplásmico (superfamilia CATH 2.70.98.90) |
| Núcleo de membrana | **PF02096** (60KD_IMP) | 351–539 | haz de 5 hélices transmembrana |

Ambos módulos corresponden a las entradas InterPro **IPR028053** (N-terminal/periplásmico)
e **IPR028055** (C-terminal/membrana), dentro de la familia **IPR001708**
(YidC/ALB3/OXA1/COX18). Se asocian los términos **GO:0032977** (actividad insertasa de
membrana), **GO:0051205** (inserción de proteínas en la membrana) y **GO:0016020**
(membrana), y la familia de transporte **TCDB 2.A.9** (YidC/Oxa1/Alb3). No se identifican
dominios catalíticos clásicos: YidC no es una enzima sino una insertasa/chaperona, y su
sitio funcional es el surco hidrofílico del núcleo transmembrana.

**En síntesis (Pregunta 2):** la proteína es una **insertasa de membrana YidC**, integral,
politópica, de **6 hélices transmembrana** (topología N-in/C-in), con **una única región
desordenada** (conector N-terminal de baja complejidad, residuos 1–69) y **dos dominios**:
el **periplásmico** (PF14849) y el **núcleo de membrana** (PF02096).
