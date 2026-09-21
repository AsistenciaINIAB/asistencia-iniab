# Flujo de trabajo para resolver el examen — basado en el Taller de Bioinfo

Este flujo usa **exactamente las herramientas de la práctica** (`Taller_Bioinfo.ipynb`)
y las mapea a las 4 preguntas del examen sobre la secuencia `datos/secuencia.fasta`.

> **Dónde correrlo:** lo más cómodo es el **mismo Colab del taller** (ya instala
> `blast+`, `mafft`, `clustalo`, `hmmer`). Los pasos con `!` son celdas de Colab; los
> pasos "web" son los servidores que usaron en clase. Subí `secuencia.fasta` al Colab.

Correspondencia herramientas del taller → preguntas:

| Pregunta | Bloque del taller | Herramientas |
|---|---|---|
| 1. Qué proteína / organismo | TP5 (BLAST) + TP6 (perfiles) | BLAST web/local, PSI‑BLAST, jackhmmer |
| 2. Membrana / desorden / dominios | Predicción de features + TP6 | DeepTMHMM, CCTOP, IUPred, DynaMine, InterPro/Pfam, CD‑search |
| 3. Estrategia posiciones relevantes | MSA + estructura | MAFFT/ClustalO, ESPript/Jalview, AlphaFold, TM‑align/Foldseek |
| 4. Ejecutar y discutir | todos | integración de resultados |

---

## Paso 0 — Preparación (celda Colab)

```bash
# (en el Colab del taller ya están instalados; si no:)
!sudo apt install -qq ncbi-blast+ mafft clustalo hmmer 2>/dev/null 1>logs.txt
!pip -q install biopython pandas pybiolib
```

```python
from Bio import SeqIO
rec = next(SeqIO.parse("secuencia.fasta","fasta"))
print(len(rec.seq), "aa")   # 534 aa
```

---

## Pregunta 1 — ¿Qué proteína es y de qué organismo? → **BLAST**

### 1a. Vía web (como en el TP5)
- NCBI **blastp** contra `nr`: https://blast.ncbi.nlm.nih.gov  → pegar la secuencia.
- **UniProt BLAST**: https://www.uniprot.org/blast  (da directo la anotación/organismo).
- EBI: https://www.ebi.ac.uk/jdispatcher/sss/ncbiblast
- Mirar: **best hit** (mayor `bitscore`/menor `E-value`), **% identidad**, **cobertura**,
  y el **organismo** (`sscinames`) → responde "qué proteína y de qué organismo".

### 1b. Vía línea de comandos (idéntico al TP5, `blastp -remote`)
```bash
!blastp -query secuencia.fasta -db nr -remote -evalue 1e-5 -max_target_seqs 50 \
  -outfmt "6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames stitle" \
  -out blastp_nr.tsv
```
```python
import pandas as pd
cols = ["qseqid","sacc","pident","length","qcovs","evalue","bitscore",
        "staxids","sscinames","stitle"]
df = pd.read_table("blastp_nr.tsv", names=cols)
df = df.sort_values("bitscore", ascending=False)
df[["sacc","pident","qcovs","evalue","bitscore","sscinames","stitle"]].head(15)
```
> **Criterio (del TP5):** homólogo confiable = `E-value` bajo + **query coverage alto**
> + identidad razonable. El `stitle` del mejor hit da el **nombre** (esperado: *YidC /
> membrane protein insertase / 60 kDa inner-membrane protein*) y `sscinames` el
> **organismo**.

### 1c. Homólogos remotos (TP6) para confirmar familia
Si el hit no es obvio, subir sensibilidad como en el taller:
- **PSI‑BLAST**: http://blast.ncbi.nlm.nih.gov  (usar el *taxonomic report*).
- **jackhmmer**: https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer

---

## Pregunta 2 — ¿Membrana? ¿Desorden? ¿Dominios?

### 2a. Proteína de membrana / topología → **DeepTMHMM** (del taller)
- Web: https://dtu.biolib.com/DeepTMHMM
- Colab:
```bash
!biolib run DTU/DeepTMHMM --fasta secuencia.fasta
```
- Complementar: **CCTOP** (http://cctop.enzim.ttk.mta.hu) y **MembraneFold**
  (https://ku.biolib.com/MembraneFold/, que mapea la topología sobre el modelo AlphaFold).
- **Reportar:** nº y posición de segmentos TM, y la topología dentro/fuera.
  (Esperado: politópica, ~6 TM, con dominio periplásmico grande.)

### 2b. Regiones desordenadas → **IUPred** (del taller)
- Web: https://iupred.elte.hu/  (modos *long* y *short*).
- Complementar: **DynaMine** (https://bio2byte.be/dynamine/) y base **MobiDB**
  (https://mobidb.bio.unipd.it/).
- **Reportar:** tramos con score > 0.5 (candidatos a desorden). En una proteína de
  membrana el desorden global es bajo, con picos en *linkers*/bucles.

### 2c. Dominios → **InterPro / Pfam + CD‑search** (TP6)
- **InterPro** (Pfam): https://www.ebi.ac.uk/interpro/  → pegar la secuencia.
- **CD‑search** (CDD): https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi
- **CDART** (arquitectura de dominios): https://www.ncbi.nlm.nih.gov/Structure/lexington/lexington.cgi
- Opcional local (como el TP6 con hmmsearch/hmmscan):
```bash
# requiere Pfam-A.hmm (hmmpress) descargado
!hmmscan --domtblout dominios.domtbl Pfam-A.hmm secuencia.fasta
```
- **Reportar:** familias/dominios (esperado: núcleo YidC/OxaA *60KD_IMP* PF02096 +
  dominio periplásmico *YidC_periplas* PF14849; InterPro IPR001708).

---

## Pregunta 3 — Estrategia para posiciones con relevancia funcional/estructural

Combina **conservación evolutiva (MSA)** + **estructura 3D**, con las herramientas del
taller:

1. **Reunir homólogos** del BLAST/PSI‑BLAST/jackhmmer (Paso 1) → guardar en
   `homologos.fasta` (secuencias diversas, varios organismos).
2. **Alineamiento múltiple** (bloque de MSA del taller), MAFFT modalidad **L‑INS‑i**:
   ```bash
   !mafft --localpair --maxiterate 1000 homologos.fasta > msa.aln    # L-INS-i
   # alternativa: !clustalo -i homologos.fasta -o msa.aln --outfmt=clustal
   ```
3. **Ver conservación por columna** con los visualizadores del taller:
   - **ESPript**: https://espript.ibcp.fr/ESPript/ESPript (resalta residuos conservados)
   - **Jalview** / **MView** (https://www.ebi.ac.uk/Tools/msa/mview)
   → las **columnas invariantes** = candidatas a relevancia funcional/estructural.
4. **Estructura 3D** (bloque de modelado del taller):
   - **AlphaFold3 server** (https://alphafoldserver.com) o **ColabFold**
     (https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
   - Usar **pLDDT** (confianza por residuo, proxy de desorden) y **PAE** (orientación
     entre dominios).
   - **MembraneFold** para ver los TM sobre el modelo.
5. **Mapear la conservación sobre la estructura** y buscar si los residuos conservados
   forman un **parche funcional** (aquí, el surco de deslizamiento del sustrato).
6. **Comparar con estructuras conocidas** (bloque de similitud estructural del taller):
   - **PDB** (https://www.rcsb.org), **TM‑align/jFATCAT** (https://www.rcsb.org/alignment),
     **DALI**, **Foldseek** (https://search.foldseek.com/search).
   → transferir por homología estructural los **residuos funcionales conocidos** de YidC.
7. **Validar el modelo** (del taller): **ProSA**, **PROCHECK** (SAVES),
   **QMEANDisCo** (SWISS‑MODEL).
8. **Integrar:** una posición es "relevante" cuando **converge** ≥2 evidencias
   (conservación alta + ubicación en surco/interfaz + contacto estructural + equivalencia
   con residuo funcional conocido).

---

## Pregunta 4 — Ejecutar y discutir

Al correr lo anterior, qué mirar y cómo redactar:

- **De BLAST:** nombre de la proteína + organismo del best hit; % identidad y cobertura.
- **De DeepTMHMM/IUPred/InterPro:** tabla con los TM, los tramos flexibles y los dominios
  (Pfam/InterPro) → responde el punto 2 con evidencia.
- **De MSA + ESPript:** listar las **columnas conservadas** y ubicarlas (¿en TM del núcleo?
  ¿en el surco? ¿en el core del dominio periplásmico?).
- **De AlphaFold + comparación estructural:** confirmar el plegado (β‑supersándwich
  periplásmico + haz de 5 TM) y señalar los residuos conservados del **surco hidrofílico**
  (p. ej. la Arg conservada esencial de YidC) como las posiciones funcionales clave.
- **Discusión:** integrar todo → la proteína es una **insertasa YidC** de bacteria
  Gram‑negativa; las posiciones relevantes son los residuos conservados del surco de
  translocación y del cinturón aromático (Trp), más el núcleo hidrofóbico del dominio
  periplásmico (rol estructural).

---

## Checklist / orden de ejecución

1. `secuencia.fasta` al Colab.
2. **BLAST** (web + `blastp -remote`, parsear con pandas) → P1.
3. **DeepTMHMM** + **IUPred** + **InterPro/CD‑search** → P2.
4. Guardar homólogos → **MAFFT L‑INS‑i** → **ESPript** (conservación) → P3.
5. **AlphaFold/ColabFold** + **Foldseek/TM‑align** + validación → P3/P4.
6. Redactar discusión integrando todas las salidas → P4.

> Nota: el `informe_analisis.md` de esta carpeta ya contiene el análisis razonado y las
> respuestas esperadas; este flujo es para **generar la evidencia** ejecutando las
> herramientas del taller.
