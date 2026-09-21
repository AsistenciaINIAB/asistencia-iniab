# Análisis bioinformático de una secuencia proteica

**Alumno:** ____________________   **Materia:** Especialidad — Taller de Bioinformática
**Fecha:** ____ / ____ / ______

---

## 1. Introducción / objetivo
_(1 párrafo: qué secuencia se analiza y con qué enfoque — identificación + caracterización estructural/funcional.)_

Secuencia analizada: **542 aa** (archivo `secuencia.fasta`).

---

## 2. Materiales y métodos
Herramientas y parámetros usados:

| Herramienta | Uso | Parámetros |
|---|---|---|
| NCBI blastp (web) | Identidad/organismo | db=nr, BLOSUM62, E=0.05 |
| DeepTMHMM | Topología de membrana | defaults |
| AIUPred (IUPred) | Desorden | umbral 0.5 |
| CD-search / InterProScan | Dominios | CDD v3.21 / InterPro 5.78 |
| MAFFT | MSA | L-INS-i (`--localpair --maxiterate 1000`) |
| ESPript | Conservación | similarity global (Risler), thr 0.7 |
| AlphaFold3 / ColabFold | Modelo 3D | recycles=3, templates on |
| Foldseek / TM-align | Comparación estructural | PDB100 + AFDB50 |

---

## 3. Resultados y discusión

### 3.1 Identificación de la proteína y organismo (Pregunta 1)
- **Best hit BLAST:** WP_323843415.1 | 99.8 % identidad | 100 % cobertura | E-value 0.0
- **Proteína:** insertasa de membrana **YidC** (familia YidC/Oxa1/Alb3)
- **Organismo:** *Moraxella* sp. Pampa (NCBI:txid 3111978, Moraxellaceae, γ-proteobacteria Gram-negativa)
- _(Figura 1: captura de BLAST. Tabla 1: hits principales.)_

### 3.2 Caracterización (Pregunta 2)

**3.2.a Membrana / topología (DeepTMHMM)**
- Proteína integral α-TM, **6 TM**, topología N-in/C-in; TM1 ancla (4–24) + dominio
  periplásmico (25–362) + núcleo 5 TM (363–530). _(Figura 2.)_

**3.2.b Regiones desordenadas (AIUPred)**
- Región desordenada única **1–69** (pico 35–60), conector de baja complejidad. _(Figura 3.)_

**3.2.c Dominios (CD-search / InterProScan)**
- Familia YidC (PRK01318, E=0); dominios **PF14849** (periplásmico, 74–346) + **PF02096**
  (núcleo, 351–539); IPR001708; GO:0032977/0051205/0016020. _(Figura 4.)_

### 3.3 Estrategia para posiciones con relevancia funcional/estructural (Pregunta 3)
_(Pipeline: homólogos del BLAST → MSA (MAFFT L-INS-i) → conservación (ESPript) → modelo 3D
(AlphaFold) → mapeo de conservación → comparación estructural (Foldseek/TM-align). Ver
`flujo_de_trabajo.md`.)_

### 3.4 Resultados del análisis y discusión (Pregunta 4)
- **Columnas conservadas del MSA:** ______________________________
- **Posiciones funcionales:** _(Arg conservada del surco; Trp del motivo GNWGW ~356–360;
  residuos hidrofílicos del núcleo TM; núcleo del β-supersándwich.)_
- _(Figura 5: ESPript. Figura 6: modelo AlphaFold. Figura 7: comparación estructural.)_

**Discusión integradora:** _(identidad → topología → dominios → posiciones clave y su rol
en la función de insertasa.)_

---

## 4. Conclusiones
- Identidad: insertasa YidC de *Moraxella* sp. Pampa.
- Membrana: politópica, 6 TM, dominio periplásmico grande.
- Desorden: solo el conector N-terminal (1–69).
- Dominios: PF14849 (periplásmico) + PF02096 (núcleo).
- Posiciones relevantes: ______________________________

---

## 5. Referencias
_(BLAST, DeepTMHMM, IUPred, InterPro, MAFFT, AlphaFold, Foldseek, y bibliografía de YidC.)_

---

## 6. Anexos
- `secuencia.fasta`, `blastp_nr.tsv`, `homologos.fasta`, `msa.aln`, `modelo.pdb/.cif`
- Capturas: `fig01_blast.png`, `fig02_deeptmhmm.png`, ...
