# Correr el análisis Donadío en local (VS Code + Claude Code)

## 1. Qué bajar del Drive

Carpeta del proyecto: **"Donadio proteomica"**
https://drive.google.com/drive/folders/1rKkyBqEewxLEKrrOtkvIA-5XXznJzsxW

Armá en tu PC esta carpeta, por ejemplo `/mnt/proyectos/Genomica_Donadio/`, y bajá adentro:

| Qué | Dónde está en el Drive | Dónde ponerlo en la PC |
|---|---|---|
| **Todo `Analisis_multiomico/`** (scripts, ambiente, resultados, figuras, datos procesados) | `Analisis_multiomico/` | `/mnt/proyectos/Genomica_Donadio/Analisis_multiomico/` |
| Matriz de conteos del proveedor | `Transcriptomica/RNAseq_ARG_Completo/Mapeamentos/ARG_Completo_Unique reads.xlsx` | `/mnt/proyectos/Genomica_Donadio/crudos/` |
| Resumen de mapeo | `Transcriptomica/RNAseq_ARG_Completo/Resumo Mapeamento ARG.ods` | `/mnt/proyectos/Genomica_Donadio/crudos/` |
| Metadatos | `Transcriptomica/Metadata Transcriptomica.docx` y `Analisis Metabolomico 2024/Metadata de Metabolomica.docx` | `/mnt/proyectos/Genomica_Donadio/crudos/` |
| **FASTQ** (los 36: 18 muestras × R1/R2, ~370 MB cada uno) | `Transcriptomica/ARG-421632644/BCLConvert_.../N-ARG-ds.*/` (**en el Drive solo están 11-ARG y 15-ARG**; pedir el resto a UFPR / BaseSpace) | `/mnt/proyectos/Genomica_Donadio/fastq/` (todos juntos, sin subcarpetas) |

Ya están dentro de `Analisis_multiomico/datos_procesados/`, así que no hace falta bajarlos aparte:
- `counts_E109.tsv` y `counts_Az39.tsv`, con los locus_tag ya corregidos;
- `anotacion_E109.tsv` y `anotacion_Az39.tsv`;
- `metab_SN_matriz_asignada.csv` y `metab_EXT_matriz_MCR.csv`, las matrices de RMN.

Los genomas se descargan solos de NCBI: GCA_000807315.1 (E109) y GCA_000632475.2 (Az39).

Necesitás tener instalados:
- **Miniforge** (conda/mamba): https://github.com/conda-forge/miniforge;
- **VS Code** con la extensión **Claude Code**;
- unos 60 GB libres para los FASTQ, los índices y los BAM;
- Linux, macOS o WSL2 en Windows. bowtie2 y samtools no corren nativos en Windows.

Los scripts también están en GitHub, en el repo `AsistenciaINIAB/asistencia-iniab`, rama `claude/laughing-johnson-z7rk3o`, carpeta `donadio/`.

---

## 2. Prompt para pegar en Claude Code (VS Code)

Abrí `/mnt/proyectos/Genomica_Donadio/` en VS Code (Archivo → Abrir carpeta), abrí Claude Code y pegá esto:

```text
Directorio de trabajo: /mnt/proyectos/Genomica_Donadio (usá rutas relativas a esta carpeta).
Contexto: tesis doctoral (Donadío, UNRC). Interacción Bradyrhizobium japonicum E109 –
Azospirillum argentinense Az39. Ya hay un análisis multi-ómico hecho; los scripts y
resultados están en ./Analisis_multiomico/ (leé primero Analisis_multiomico/LEEME_resultados_multiomica.md
y Analisis_multiomico/hoja_muestras_Donadio.tsv). Trabajá en español. No modifiques los
resultados existentes: escribí todo lo nuevo en ./Analisis_multiomico/v2_local/.

Diseño RNA-seq (Illumina pareado, 4 h, 3 réplicas): 1-3 E109 control | 4-6 E109+AIA |
7-9 E109+Az39 | 10-12 E109+Az39-ipdC | 13-15 Az39 sola (MMAB) | 16-18 Az39-ipdC sola.
FASTQ en ./fastq/ con nombres tipo 11-ARG_S19_L001_R1_001.fastq.gz.

TAREA 1 – Ambiente
- Creá el ambiente con: mamba env create -f Analisis_multiomico/envs/donadio-multiomica.yml
  y agregale: fastp bowtie2 samtools subread (featureCounts) bedtools ncbi-datasets-cli,
  y en R: pathview, KEGGREST. Si DESeq2 no carga por GenomeInfoDbData, instalalo desde Bioconductor.

TAREA 2 – Re-alineamiento competitivo (prioridad)
1. Bajá de NCBI (datasets) FASTA + GFF de GCA_000807315.1 (E109) y GCA_000632475.2 (Az39).
2. Armá UN genoma combinado (E109 + Az39) y un GFF combinado; índice bowtie2.
3. Por muestra: fastp (QC, recorte) → bowtie2 --very-sensitive contra el combinado
   (así cada lectura va a la bacteria donde mejor encaja; no secuencial) → samtools sort/index.
   Reportá también el % de ARNr (el proveedor reportó 58–95%).
4. featureCounts -p --countReadPairs -t gene -g locus_tag (probá -s 0/1/2 en 2 muestras
   y elegí el strandedness con más asignación); MAPQ >= 10 para descartar ambiguas.
5. Salidas: counts_E109_v2.tsv y counts_Az39_v2.tsv (mismo formato que datos_procesados/),
   tabla de QC por muestra, y lecturas por organismo y por replicón.
6. Validaciones:
   a) cobertura por posición en ABAZ39_11295 (ipdC) en muestras WT vs ipdC: ubicar la
      inserción/deleción (gráfico de cobertura);
   b) cobertura del plásmido p3 (CP007796.1) en las 12 muestras de Az39: confirmar que la
      cepa ipdC no lo tiene;
   c) correlación de los nuevos conteos con los del proveedor.
7. Re-corré Analisis_multiomico/scripts/01_deseq2_E109.R y 03_deseq2_Az39.R apuntando a los
   conteos v2 (copiá data/ a v2_local/data/) y compará con los resultados previos: ¿se
   mantienen la isla iac (RN69_30050-30290), la caída de ribosomas en E109 en co-cultivo
   y la clasificación AIA-dependiente/independiente?

TAREA 3 – Vías KEGG con transcriptómica + metabolómica
1. Con KEGGREST: buscá si E109 y Az39 tienen código de organismo KEGG y si sus genes usan
   los locus_tag RN69_/ABAZ39_. Si no, asigná KO con KofamScan (bajar perfiles de genome.jp).
2. Mapeá los metabolitos asignados de la RMN (sobrenadante y extractos) a IDs KEGG
   (C00078 triptófano, C01083 trehalosa, C00149 malato, C00122 fumarato, etc.); guardá la tabla.
3. Enriquecimiento por vía KEGG: genes DE (padj<0.05,|log2FC|>=1) y metabolitos DA (FDR<0.05)
   por contraste (E109: Az39 vs control, AIA vs control, Az39 vs Az39-ipdC; metabolitos:
   EA vs E, Eaia vs E), p hipergeométrico por capa y p combinado (Fisher); FDR BH.
4. Gráfico resumen: dotplot de vías × capa (genes / metabolitos), tamaño = n, color = log2FC medio o -log10 FDR.
5. Mapas pathview (genes + compuestos coloreados por log2FC) para: map00380 (triptófano/AIA),
   map00020 (TCA), map00500 (almidón/sacarosa/trehalosa), map00190 (fosforilación oxidativa),
   map03010 (ribosoma), map00290 (aminoácidos ramificados), y las 5 vías más significativas.

Reglas: scripts numerados en v2_local/scripts/, todo reproducible desde línea de comandos,
semilla fija, sessionInfo. Al final escribí v2_local/LEEME_v2.md con qué cambió respecto
del análisis anterior, números concretos y advertencias. Preguntame antes de borrar algo.
```

---

## 3. Si solo tenés los FASTQ de 11-ARG y 15-ARG

Igual sirve para un chequeo rápido. Pedile a Claude Code que haga solo los pasos 2.1–2.3 y la validación 6b con esas dos muestras. La 11-ARG lleva Az39-*ipdC* y la 15-ARG Az39 salvaje, así que alcanzan para confirmar si falta p3 y para ver la cobertura de *ipdC*.
