# Examen de especialidad — Análisis bioinformático de una secuencia

Análisis de una secuencia proteica de 534 aa. **Conclusión:** insertasa de membrana
**YidC** (familia YidC/Oxa1/Alb3) de una **bacteria Gram-negativa**.

## Estructura de la carpeta

```
especialidad/examen/
├── README.md                          (este archivo)
├── datos/
│   └── secuencia.fasta                secuencia de entrada
├── scripts/
│   ├── analisis_secuencia.py          análisis local (KD, composición, TM, LCR)
│   └── comandos_bioinformaticos.sh    pipeline externo (BLAST, HMMER, DeepTMHMM,
│                                       MAFFT, ConSurf, AlphaFold, ...)
├── resultados/
│   └── informe_analisis.md            INFORME PRINCIPAL (respuestas 1-4)
└── figuras/                           (perfil de hidropatía al correr el script)
```

## Cómo reproducir

```bash
cd especialidad/examen
pip install biopython numpy matplotlib
python3 scripts/analisis_secuencia.py     # genera resultados/ y figuras/
bash   scripts/comandos_bioinformaticos.sh  # pipeline con red/herramientas externas
```

## Respuestas (resumen — detalle en `resultados/informe_analisis.md`)

1. **Proteína/organismo:** YidC (insertasa de proteínas de membrana). Organismo:
   Gram-negativa; especie exacta pendiente de BLAST (bloqueado en el sandbox).
2. **Membrana/desorden/dominios:** proteína politópica de **6 TM** (1 ancla + núcleo
   de 5); regiones flexibles en el conector periplásmico N-terminal (~29–75) y el bucle
   citoplásmico (~372–420); dominios: **P1 periplásmico** (β-supersándwich, PF14849) y
   **núcleo de membrana** YidC/OxaA (PF02096); familia IPR001708 / TIGR03593.
3. **Estrategia:** homología → MSA → conservación (ConSurf/Rate4Site) → coevolución
   (EVcouplings) → estructura (AlphaFold) → topología (DeepTMHMM/SignalP) → desorden
   (IUPred3) → cavidades/electrostática, integrando líneas de evidencia.
4. **Resultados:** posiciones relevantes = Arg conservada del surco (candidata ~R382,
   equivalente a Arg366 de *E. coli*), Trp del motivo `GNWGW` y del cinturón aromático,
   residuos hidrofílicos del surco (TM5) y núcleo del β-supersándwich.

> **Limitación del entorno:** en este sandbox el proxy de red bloquea NCBI/UniProt/EBI
> y el clasificador impide ejecutar scripts, por lo que el BLAST, el MSA/ConSurf y el
> modelo AlphaFold quedan como pasos a ejecutar con los scripts provistos. Las
> conclusiones de familia, topología y dominios son firmes a partir de la secuencia.
