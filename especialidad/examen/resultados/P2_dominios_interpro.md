# Pregunta 2 — Dominios (InterProScan)

**Herramienta:** InterProScan 5.78-109.0 (EBI). **Entrada:** 542 aa.

## Familia (cubre la proteína completa)
| Firma | Base | Rango | E-value | Entrada |
|---|---|---|---|---|
| NF002352 / **PRK01318** (membrane protein insertase YidC) | NCBIfam | 3–540 | **0.0** | — |
| MF_01810 *YidC_type1* | HAMAP | 5–538 | — | IPR019998 (Membrane insertase YidC) |
| PTHR12428 (OXA1) | PANTHER | 336–535 | 1.5E-54 | IPR001708 |
| PR00701 (60kDa inner membrane protein) + PR01900 (YidC signature) | PRINTS | varios | 1.1E-77 / 4.1E-22 | IPR001708 |

GO asociados: **GO:0032977** (membrane insertase activity), **GO:0016020** (membrane),
**GO:0051205** (protein insertion into membrane).

## Dominio periplásmico (P1, β-supersándwich)
| Firma | Base | Rango | E-value | InterPro |
|---|---|---|---|---|
| **PF14849** *YidC_periplas* | Pfam | 74–346 | 4.6E-73 | IPR028053 (Membr_insert_YidC_N) |
| TIGR03593 *yidC_nterm* | NCBIfam | 6–356 | 3.0E-97 | IPR028053 |
| G3DSA:2.70.98.90 | Gene3D | 70–338 | 8.7E-72 | IPR038221 (YidC periplasmic domain superfamily) |
| cd19961 *EcYidC-like_peri* | CDD | 72–334 | 1.8E-70 | IPR028053 |

## Núcleo transmembrana (C-terminal, 5 TM)
| Firma | Base | Rango | E-value | InterPro |
|---|---|---|---|---|
| **PF02096** *60KD_IMP* | Pfam | 351–539 | 5.1E-76 | IPR028055 (YidC/Oxa/ALB_C) |
| TIGR03592 *yidC_oxa1_cterm* | NCBIfam | 357–536 | 1.6E-65 | IPR028055 |
| cd20070 *5TM_YidC_Alb3* | CDD | 359–536 | 9.3E-98 | IPR047196 (YidC_ALB_C) |

## Topología y señal (InterPro)
- **SignalP-TM**: 1–28 (score 0.538) → TM1 = **señal-ancla** (no péptido señal escindible).
- TMHMM/Phobius (InterPro): 5 TM (7–24, 354–376, 423–445, 465–487, 500–522) + gran región
  periplásmica 25–357. **DeepTMHMM dio 6 TM** (separa el último segmento); la referencia
  para YidC Gram-negativa es 6 TM del núcleo.

## Conclusión
InterProScan confirma la **arquitectura bidominio de YidC**: **dominio periplásmico**
(PF14849 / IPR028053, fold Gene3D 2.70.98.90) + **núcleo de membrana** (PF02096 / IPR028055),
todo dentro de la familia **YidC/OxaA (IPR001708)** y el modelo bacteriano **PRK01318**
(E-value = 0). Coincide con BLAST, DeepTMHMM y CD-search.

## Texto listo para el informe (sección 3.2.c)
> InterProScan asigna la proteína a la familia **insertasa de membrana YidC** (modelo
> bacteriano PRK01318/NF002352, residuos 3–540, E-value = 0; HAMAP MF_01810 → IPR019998;
> PANTHER/PRINTS → IPR001708). Se resuelven dos dominios estructurales: el **dominio
> periplásmico** (Pfam **PF14849**, 74–346; superfamilia Gene3D **2.70.98.90** →
> IPR038221) y el **núcleo transmembrana** C-terminal (Pfam **PF02096**, 351–539; CDD
> cd20070 *5TM_YidC_Alb3*), correspondientes a IPR028053 e IPR028055. Se asocian las
> funciones GO:0032977 (actividad insertasa), GO:0051205 (inserción de proteínas en
> membrana) y GO:0016020 (membrana). SignalP indica que la hélice N-terminal (1–28) actúa
> como señal-ancla, no como péptido señal escindible.
