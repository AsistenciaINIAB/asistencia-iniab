# Pregunta 1 — Identidad y organismo (BLASTp web vs nr)

**Herramienta:** BLASTp (NCBI web, base `nr`). **Entrada:** `secuencia.fasta` (542 aa).
**Archivo:** `B2YNZVH0016-Alignment-HitTable.csv` (~100 hits).

## Mejor hit (identidad resuelta)

| sacc | % id | cobertura | mismatch | gaps | E-value | bitscore | % positives |
|---|---|---|---|---|---|---|---|
| **WP_323843415.1** | **99.815** | 542/542 (100%) | 1 | 0 | **0.0** | **1116** | 100 |

- La proteína query es **prácticamente idéntica** a `WP_323843415.1` (difieren en **1 residuo**).
  → **Identidad de la proteína: resuelta.**
- `WP_` = registro RefSeq **multiespecie** (proteína idéntica en varias cepas/especies).

## Distribución de hits (familia conservada)
- ~100 homólogos con identidad **99.8% → ~54%**, **E-value = 0.0** en todos, cobertura ≈540 aa.
- Indica una **proteína muy conservada y ampliamente distribuida en bacterias** → coherente
  con **YidC** (insertasa de membrana). Además provee el **set de homólogos para el MSA (P3)**.

## Organismo
- **Organismo:** *Moraxella* sp. Pampa — NCBI Taxonomy ID **3111978** (ncbitaxon:3111978).
- Clasificación: **γ-proteobacteria**, familia **Moraxellaceae** → **bacteria Gram-negativa**
  (coherente con la topología YidC: TM1 de anclaje + dominio periplásmico grande + 5 TM).

## Conclusión (integra con DeepTMHMM)
Identidad ≈100% con `WP_323843415.1` (de *Moraxella* sp. Pampa) + topología de 6 TM con
dominio periplásmico grande (DeepTMHMM) ⇒ la proteína es la **insertasa de membrana YidC**
(familia YidC/Oxa1/Alb3) de ***Moraxella* sp. Pampa**, una γ-proteobacteria Gram-negativa.

## Texto listo para el informe (sección 3.1)
> La búsqueda BLASTp contra `nr` identificó como mejor hit a `WP_323843415.1` con
> **99.8 % de identidad, 100 % de cobertura y E-value = 0.0** (bitscore 1116); la proteína
> analizada es, por lo tanto, esencialmente idéntica a esa entrada, anotada como la
> **insertasa de membrana YidC** de ***Moraxella* sp. Pampa** (NCBI:txid 3111978, familia
> Moraxellaceae, γ-proteobacteria Gram-negativa). El resto de los ~100 hits (identidad
> 54–99 %, todos con E-value = 0.0) corresponden a la misma familia proteica, muy conservada
> en bacterias, y constituyen el conjunto de homólogos empleado en el alineamiento múltiple.
