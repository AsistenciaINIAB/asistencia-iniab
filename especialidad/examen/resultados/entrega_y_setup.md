# Qué entregar, cómo presentarlo y setup de cada servidor web

Guía práctica para armar la entrega del examen a partir del análisis de
`datos/secuencia.fasta` (insertasa **YidC**), usando las herramientas del Taller.

---

## PARTE A — Qué tenés que entregar (checklist de entregables)

### A.1 Documento principal (informe, PDF)
Un informe que responda **las 4 preguntas**, con **evidencia** (capturas + tablas) de
cada herramienta. Ver estructura en la Parte B.

### A.2 Evidencia por pregunta (figuras/capturas + archivos)

| Pregunta | Evidencia a incluir | Archivo/figura |
|---|---|---|
| **1. Proteína/organismo** | Captura del BLAST con el best hit; tabla de hits (acc, % id, cobertura, E‑value, organismo) | `blastp_nr.tsv`, captura NCBI/UniProt |
| **2. Membrana** | Plot y tabla de topología | figura DeepTMHMM (+ CCTOP) |
| **2. Desorden** | Perfil de desorden con umbral 0.5 | figura IUPred3 |
| **2. Dominios** | Diagrama de dominios + tabla (Pfam/InterPro/CDD) | captura InterPro + CD‑search |
| **3. Estrategia** | Texto del pipeline (no requiere figura) | — |
| **4. Resultados** | MSA con conservación; modelo 3D; alineamiento estructural | figura ESPript, imagen AlphaFold, salida Foldseek/TM‑align |

### A.3 Archivos anexos (carpeta `anexos/`)
- `secuencia.fasta` — secuencia de entrada.
- `blastp_nr.tsv` — tabla de resultados BLAST.
- `homologos.fasta` — secuencias homólogas usadas para el MSA.
- `msa.aln` — alineamiento múltiple (Clustal/FASTA).
- `modelo.pdb` (o `.cif`) — modelo de AlphaFold.
- Capturas de pantalla de cada servidor (numeradas: `fig01_blast.png`, ...).

> **Regla de oro:** cada afirmación del informe debe poder rastrearse a una figura o
> archivo. "La proteína tiene 6 TM" → figura DeepTMHMM. "El organismo es X" → captura BLAST.

---

## PARTE B — Estructura de presentación del informe

```
1. Portada
   - Título, materia (Especialidad – Taller de Bioinformática), alumno, fecha.

2. Introducción / objetivo (1 párrafo)
   - Qué se analiza y con qué enfoque (identificación + caracterización estructural).

3. Materiales y métodos
   - Secuencia analizada (longitud, archivo).
   - Lista de herramientas y VERSIONES/PARÁMETROS usados (tabla; ver Parte C).

4. Resultados y discusión (una subsección por pregunta)
   4.1 Identificación de la proteína y organismo (Pregunta 1)
       - Best hit BLAST, % identidad, cobertura, E-value, organismo. Figura + tabla.
   4.2 Caracterización (Pregunta 2)
       - 4.2.a ¿Es de membrana? Topología (DeepTMHMM). Tabla de TM + figura.
       - 4.2.b Regiones desordenadas (IUPred). Figura.
       - 4.2.c Dominios (InterPro/Pfam + CD-search). Diagrama + tabla.
   4.3 Estrategia bioinformática para posiciones relevantes (Pregunta 3)
       - Descripción del pipeline (MSA→conservación→estructura→comparación).
   4.4 Resultados del análisis y discusión (Pregunta 4)
       - MSA + conservación (ESPript). Modelo 3D (AlphaFold) + comparación
         estructural (Foldseek/TM-align). Posiciones funcionales identificadas.

5. Conclusiones (bullets)
   - Identidad, topología, dominios, posiciones clave.

6. Referencias (BLAST, DeepTMHMM, IUPred, InterPro, MAFFT, AlphaFold, etc.).

7. Anexos
   - Archivos y capturas.
```

> **Consejo de presentación:** numerá figuras y tablas y referencialas en el texto
> ("Figura 1", "Tabla 2"). Poné en cada figura un pie con la herramienta y los
> parámetros usados.

---

## PARTE C — Setup (parámetros) para cada servidor web

### C.1 BLAST — identidad y organismo (Pregunta 1)

**NCBI blastp** — https://blast.ncbi.nlm.nih.gov (pestaña *Protein BLAST*)
- **Query:** pegar la secuencia (o subir `secuencia.fasta`).
- **Database:** `nr` (non-redundant). Para anotación más limpia probá también
  `refseq_protein` o `swissprot`.
- **Program selection / Algorithm:** `blastp` (quick es suficiente).
- **Algorithm parameters:**
  - Max target sequences: **100**
  - Expect threshold (E-value): **0.05** (dejar el default; para más rigor 1e‑5)
  - Matrix: **BLOSUM62**; Gap costs: **Existence 11, Extension 1**
  - Composition-based statistics: **Conditional (default)**
  - Filtro de baja complejidad: **desactivado** para no perder el TM/LCR (opcional).
- **Para el organismo:** en la vista de resultados abrí **"Taxonomy"** / *Taxonomy reports*
  y mirá la columna *Scientific Name* del best hit.

**UniProt BLAST** — https://www.uniprot.org/blast
- Target database: **UniProtKB** (o Swiss-Prot para curado).
- E-Threshold: **10** (default), Matrix: **Auto/BLOSUM62**, Hits: **50**.
- Ventaja: te da directo el **nombre curado** y el **organismo**.

**PSI-BLAST** (si necesitás más sensibilidad / homólogos remotos)
- http://blast.ncbi.nlm.nih.gov → elegí **PSI-BLAST**.
- PSI-BLAST threshold (inclusión): **0.005**; iterar 2–3 rondas hasta convergencia.
- Sirve además para **reclutar homólogos** para el MSA (Pregunta 3).

**jackhmmer** — https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer
- Database: **UniProtKB** (o Reference Proteomes).
- Significance E-values: seq **0.01**, hit **0.03**; iterar hasta convergencia.

---

### C.2 Topología de membrana (Pregunta 2)

**DeepTMHMM** — https://dtu.biolib.com/DeepTMHMM
- Pegar la secuencia / subir FASTA. **Sin parámetros** (predice TM + dentro/fuera).
- Descargá el **plot de probabilidad** y el archivo de topología (`.gff3`/`.md`).

**CCTOP** (complemento) — http://cctop.enzim.ttk.mta.hu
- Submit sequence; dejar métodos y fuentes por defecto (usa consenso + restricciones).

**MembraneFold** — https://ku.biolib.com/MembraneFold/
- Mapea la topología sobre el modelo AlphaFold (útil para la figura de la Pregunta 4).

---

### C.3 Desorden (Pregunta 2)

**IUPred3** — https://iupred.elte.hu/
- **Analysis type:** IUPred3.
- **Prediction type:** correr **long disorder** y también **short disorder**.
- **ANCHOR2:** activarlo (marca regiones de unión dependientes de contexto).
- **Smoothing:** medium (default). Umbral de desorden: **0.5**.

**DynaMine** — https://bio2byte.be/dynamine/ (flexibilidad del backbone, S²).
**MobiDB** — https://mobidb.bio.unipd.it/ (consultá si ya hay anotación de IDRs).

---

### C.4 Dominios (Pregunta 2)

**InterPro** — https://www.ebi.ac.uk/interpro/ (corre InterProScan completo)
- Pegar secuencia → devuelve **Pfam, PROSITE, PRINTS, CDD, TIGRFAM, SUPERFAMILY,
  Gene3D, PANTHER**, GO y modelo AlphaFold. **Sin parámetros.**
- Reportá los **accession** (p. ej. Pfam PF02096, PF14849; InterPro IPR001708).

**CD-search (CDD)** — https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi
- **Database:** CDD (default) — opcional Pfam/SMART/TIGRFAM/COG.
- **Expect value threshold:** **0.01**.
- **Composition-corrected scoring:** ON. **Results mode:** Full (para ver todos los hits).

**CDART** — https://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi
- Muestra proteínas con la **misma arquitectura de dominios**.

---

### C.5 MSA — para conservación (Pregunta 3)

**MAFFT** (modalidad **L-INS-i**, la más precisa para pocas secuencias)
- Server: https://mafft.cbrc.jp/alignment/server/  → elegir **L-INS-i**.
- Línea de comandos (Colab del taller): `mafft --localpair --maxiterate 1000 in.fasta > msa.aln`
- Matriz: BLOSUM62; gap open 1.53 (defaults de L-INS-i).

**Clustal Omega** — https://www.ebi.ac.uk/Tools/msa/clustalo/
- Input: `homologos.fasta`; Output format: **Clustal** y **FASTA** (ambos).
- Dejar parámetros por defecto (n.º de iteraciones default).

**T-Coffee / MUSCLE** (para comparar): EBI, defaults.

---

### C.6 Visualización de conservación (Pregunta 3/4)

**ESPript 3** — https://espript.ibcp.fr/ESPript/ESPript
- Subir el **MSA alineado** (FASTA/Clustal/MSF).
- (Opcional) subir un **PDB** para superponer estructura secundaria.
- **Similarity:** *global score* con **%Equivalent (Risler)** o BLOSUM62;
  **similarity threshold ≈ 0.7**.
- Render: residuos conservados en **cajas rojas**; exportar a PDF/PNG.

**Jalview** / **MView** (https://www.ebi.ac.uk/Tools/msa/mview): alternativas de coloreo.

---

### C.7 Estructura 3D (Pregunta 3/4)

**AlphaFold3 server** — https://alphafoldserver.com
- Job type: **single protein**; pegar la secuencia. Descargar `.cif` + pLDDT/PAE.

**ColabFold (AlphaFold2_mmseqs2)** —
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- **msa_mode:** `mmseqs2_uniref_env`
- **model_type:** `auto` (alphafold2_ptm para monómero)
- **num_recycles:** `3` (subir a 6 si el pLDDT es bajo)
- **templates:** activar (pdb100) si querés usar molde de YidC conocido
- **amber/relax:** opcional (mejora geometría local)
- Mirar **pLDDT** (confianza; bajo = flexible/desordenado) y **PAE** (orientación P1↔núcleo).

**Validación del modelo** (del taller):
- ProSA — https://prosa.services.came.sbg.ac.at/prosa.php (Z-score)
- PROCHECK / SAVES v6 — https://saves.mbi.ucla.edu/ (Ramachandran)
- QMEANDisCo — https://swissmodel.expasy.org/qmean/

---

### C.8 Comparación estructural (Pregunta 3/4)

**RCSB PDB pairwise** — https://www.rcsb.org/alignment
- Programa: **jFATCAT-rigid** y **jFATCAT-flexible**; también **TM-align**.
- Reportar **RMSD** y **TM-score** contra estructuras YidC (buscar en PDB por "YidC").

**Foldseek** — https://search.foldseek.com/search
- Subir el **modelo AlphaFold** (`.pdb/.cif`).
- **Databases:** PDB100 + AlphaFold/UniProt50 (o AlphaFold/Proteome).
- **Mode:** 3Di/AA (default); sensibilidad default.
- Sirve para **identificar homólogos estructurales** y confirmar la familia/organismo.

**DALI** — http://ekhidna2.biocenter.helsinki.fi/dali/ (alternativa).

---

## Resumen de parámetros clave (para la tabla de Métodos)

| Herramienta | Parámetro principal |
|---|---|
| NCBI blastp | db=nr, BLOSUM62, gap 11/1, E=0.05, max_target=100 |
| PSI-BLAST | inclusión E=0.005, 2–3 iteraciones |
| DeepTMHMM | defaults |
| IUPred3 | long + short, ANCHOR2 on, umbral 0.5 |
| CD-search | CDD, E=0.01, comp-corrected ON |
| MAFFT | L-INS-i (`--localpair --maxiterate 1000`) |
| ESPript | similarity global (Risler), threshold 0.7 |
| ColabFold | msa=mmseqs2_uniref_env, recycles=3, templates on |
| Foldseek | PDB100 + AFDB50, 3Di/AA |
