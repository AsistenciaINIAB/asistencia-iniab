# Pregunta 2 — Regiones desordenadas (AIUPred)

**Herramienta:** AIUPred v0.1 (IUPred; Erdős & Dosztányi). **Entrada:** secuencia (542 aa).
**Umbral de desorden:** 0.5. **Smoothing:** default.

## Resultado
- **Una sola región desordenada: residuos 1–69** (score > 0.5 de forma continua; cae a
  0.45 en el residuo 70).
- **Pico de desorden (score > 0.90): residuos ~35–60** → segmento de **baja complejidad
  rico en Thr/Ala/Ser/Asp/Asn** (`...ASATTQAATIDLPNANAGDVPTTTNASSDPATTALDGQ...`).
- **Resto de la proteína (70–542): ordenada** (scores < 0.5, mayormente < 0.3): dominio
  periplásmico P1 (β-supersándwich) + núcleo de 5 TM.
- Cola C-terminal (~535–542): leve aumento (0.40–0.43) sin cruzar 0.5 → ligera flexibilidad.

## Interpretación
Desorden **global bajo** (típico de proteína de membrana). Se concentra en un **conector
flexible de baja complejidad** (T/A/S-rico) que une el **ancla TM1** con el **dominio
periplásmico plegado (P1)**. Matiz: la ventana 1–69 solapa el TM1 (hidrofóbico); la señal
de desorden *bona fide* es el **linker ~30–60**, no la hélice de anclaje.

## Texto listo para el informe (sección 3.2.b)
> El análisis con AIUPred (umbral 0.5) muestra un **contenido de desorden global bajo**,
> coherente con una proteína integral de membrana. La única región con score > 0.5 es el
> **extremo N-terminal (residuos 1–69)**, con un máximo de desorden (> 0.90) en el segmento
> **~35–60**, de **baja complejidad y rico en Thr/Ala/Ser/Asp/Asn**. Esta región constituye
> un **conector flexible** entre la hélice de anclaje (TM1) y el dominio periplásmico
> plegado. El resto de la proteína (dominio periplásmico y núcleo transmembrana) es
> estructuralmente ordenada.
