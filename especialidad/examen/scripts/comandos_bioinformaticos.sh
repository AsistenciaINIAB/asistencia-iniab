#!/usr/bin/env bash
# =============================================================================
# Pipeline reproducible para el análisis de la secuencia (YidC).
# Requiere red hacia NCBI/EBI/servidores (bloqueada en el sandbox del examen).
# Ejecutar desde: especialidad/examen/
# =============================================================================
set -euo pipefail
Q=datos/secuencia.fasta
mkdir -p resultados

# --- 1. IDENTIDAD Y ORGANISMO: BLASTp contra nr (línea de comandos) ----------
# Local (con base nr descargada):
#   blastp -query "$Q" -db nr -evalue 1e-5 -max_target_seqs 25 \
#     -outfmt "6 sacc pident length evalue stitle staxids sscinames" \
#     -out resultados/blastp_nr.tsv
# Alternativa online: web BLAST (https://blast.ncbi.nlm.nih.gov) o la API URLAPI.
# El best hit da género/especie -> responde el punto 1 (organismo).

# --- 2. DOMINIOS / FAMILIAS: HMMER + Pfam e InterProScan ----------------------
#   hmmscan --domtblout resultados/pfam.domtbl Pfam-A.hmm "$Q"
#   interproscan.sh -i "$Q" -f tsv,gff3 -d resultados/   # Pfam, TIGRFAM, PROSITE...
# Esperado: PF02096 (60KD_IMP, núcleo YidC), PF14849 (YidC_periplas),
#           IPR001708 (OxaA/YidC), TIGR03593 (yidC).

# --- 3. TOPOLOGÍA DE MEMBRANA / PÉPTIDO SEÑAL --------------------------------
# DeepTMHMM (biolib):   biolib run DTU/DeepTMHMM --fasta "$Q"
# Phobius / TOPCONS:     servidores web (TM + dentro/fuera de forma consensuada)
# SignalP-6:             signalp6 --fastafile "$Q" --organism other --output_dir resultados/sp6
#   -> TM1 debe salir como señal-ancla (no escindible), 6 TM totales.

# --- 4. DESORDEN / BAJA COMPLEJIDAD -----------------------------------------
# IUPred3:     python3 iupred3.py "$Q" long > resultados/iupred3.txt
# metapredict: metapredict-predict-disorder "$Q" -o resultados/metapredict.csv
#   -> picos en el conector N-terminal (~29-75) y el bucle citoplásmico (~372-420).

# --- 5. MSA + CONSERVACIÓN (posiciones funcionales) -------------------------
# 5a. Homólogos:
#   jackhmmer -N 3 --tblout resultados/hits.tbl "$Q" uniref90.fasta
#   (o phmmer / DIAMOND) -> extraer secuencias homólogas a homologos.fasta
# 5b. Alineamiento:
#   mafft --auto homologos.fasta > resultados/msa.aln         # o L-INS-i
# 5c. Conservación por residuo:
#   rate4site -s resultados/msa.aln -o resultados/rate4site.txt
#   (o subir el MSA a ConSurf: https://consurf.tau.ac.il)
#   -> residuos invariantes = candidatos funcionales (surco/greasy slide).

# --- 6. COEVOLUCIÓN / CONTACTOS ---------------------------------------------
# EVcouplings (https://evcouplings.org) o GREMLIN sobre el MSA -> pares
#   coevolucionantes (contactos, acople P1-núcleo).

# --- 7. ESTRUCTURA 3D + MAPEO DE CONSERVACIÓN -------------------------------
# AlphaFold2/3 (ColabFold):  colabfold_batch "$Q" resultados/af/
#   -> modelo + pLDDT (proxy de desorden) + PAE (orientación entre dominios).
# Mapear ConSurf/rate4site sobre el modelo (PyMOL/ChimeraX) -> parche funcional.
# Cavidades/electrostática:  fpocket -f modelo.pdb ; APBS para el surco hidrofílico.
# Superposición con YidC de referencia (E. coli, B. halodurans):
#   TM-align modelo.pdb 3wo6.pdb   # transferir residuos funcionales (Arg del surco).

echo "Pipeline documentado. Ejecutar cada bloque según disponibilidad de red/herramientas."
