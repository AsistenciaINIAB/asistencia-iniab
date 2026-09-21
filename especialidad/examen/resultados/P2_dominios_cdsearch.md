# Pregunta 2 — Dominios (CD-search / CDD)

**Herramienta:** CD-search (NCBI), base **CDD v3.21**, E-value cutoff 0.01,
composition-based adjustment ON, low-complexity filter OFF. **Entrada:** 542 aa.

## Resultado
**Clasificación:** *membrane protein insertase YidC* (ArchID 11479554) —
"funciona como chaperona de membrana y como insertasa independiente de proteínas de membrana".

| Ítem | Valor |
|---|---|
| Hit específico (dominio) | **PRK01318** — *membrane protein insertase* — intervalo **1–540**, **E-value = 0** |
| Gene Symbol | **yidC** |
| Gene Ontology | GO:0032977 (membrane insertase activity) · GO:0016020 (membrane) |
| TCDB | **2.A.9** (familia YidC/Oxa1/Alb3) |
| CATH | 2.70.98.90 (fold del dominio periplásmico) |
| SCOP | 4007001 · 4006662 |
| Feature | "oligomer interface" en el núcleo TM (C-terminal) |

- El modelo **PRK01318** cubre **toda la secuencia (1–540)** con **E-value = 0**: dominio
  único de la familia YidC que engloba el dominio periplásmico + el núcleo transmembrana.
- La anotación *oligomer interface* en el núcleo TM es consistente con la funcionalidad
  dimérica descrita para YidC.

## Equivalencia con Pfam/InterPro
- InterPro: **IPR001708** (Membrane insertion protein OxaA/YidC).
- Pfam: **PF02096** (60KD_IMP, núcleo de membrana) + **PF14849** (YidC_periplas, periplásmico).
- En CDD/CD-search aparece condensado como el cluster **PRK01318**.

## Texto listo para el informe (sección 3.2.c)
> El análisis con CD-search (CDD v3.21) clasifica a la proteína como **insertasa de
> membrana YidC** (gen *yidC*), con un único dominio conservado, **PRK01318** (*membrane
> protein insertase*), que abarca casi toda la secuencia (**1–540, E-value = 0**). La
> anotación asocia las funciones **GO:0032977** (actividad insertasa de membrana) y
> **GO:0016020** (membrana), la familia de transporte **TCDB 2.A.9** (YidC/Oxa1/Alb3) y el
> fold **CATH 2.70.98.90** del dominio periplásmico. En términos de Pfam/InterPro, el
> dominio corresponde a **IPR001708** y a los módulos **PF02096** (núcleo transmembrana)
> y **PF14849** (dominio periplásmico). Se anota además una *oligomer interface* en el
> núcleo, coherente con la funcionalidad dimérica de YidC.
