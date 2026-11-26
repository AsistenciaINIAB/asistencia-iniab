#!/usr/bin/env Rscript
## =============================================================================
## 08 - Enriquecimiento funcional (ORA hipergeométrico, FDR BH) con anotación COG
## Anotación: COGclassifier v2 (RPS-BLAST vs CDD/COG, e-value 0.01) sobre los
## proteomas NCBI de E109 (GCA_000807315.1) y Az39 (GCA_000632475.2).
## Dos niveles: categoría funcional COG (letra) y vía COG (cog-24.def.tab).
## Universo: genes testeados en DESeq2 (pasaron el filtro) con anotación COG.
## Conjuntos: genes DE (padj<0.05, |log2FC|>=1) por sentido en cada contraste,
## y clases de la respuesta de E109 a Az39 (AIA-dependiente / independiente).
## Uso: Rscript 08_enriquecimiento_COG.R <wd> <dir_cog> <cog-24.def.tab> <cog-24.fun.tab> <cog-24.pathways.tab>
##   (archivos de https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/)
## =============================================================================
suppressMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
wd <- args[1]; cogdir <- args[2]; deffile <- args[3]; funfile <- args[4]; pathfile <- args[5]
res_dir <- file.path(wd, "results", "enriquecimiento"); fig_dir <- file.path(wd, "figures", "enriquecimiento")
dir.create(res_dir, TRUE, TRUE); dir.create(fig_dir, TRUE, TRUE)
ALPHA <- 0.05; LFC <- 1; MIN_TERM <- 5

def <- read.delim(deffile, header = FALSE, quote = "", col.names = c("COG_ID", "letras", "nombre", "gen", "via", "pubmed", "pdb"))
fl <- strsplit(readLines(funfile), "\t"); fl <- fl[lengths(fl) == 4 & nchar(sapply(fl, `[`, 1)) == 1]
catname <- setNames(paste0(sapply(fl, `[`, 1), ": ", sapply(fl, `[`, 4)), sapply(fl, `[`, 1))
pw <- read.delim(pathfile, header = FALSE, quote = "")[, 1:2]; names(pw) <- c("via", "COG_ID")

leer_cog <- function(org) {
  f <- list.files(file.path(cogdir, paste0("cog_", org)), pattern = "classif.*\\.tsv$", full.names = TRUE)[1]
  x <- read.delim(f, quote = "", check.names = FALSE)
  names(x) <- toupper(names(x))
  x <- x[!duplicated(x$QUERY_ID), ]
  data.frame(locus_tag = x$QUERY_ID, COG_ID = x$COG_ID)
}
anot <- list(E109 = leer_cog("E109"), Az39 = leer_cog("Az39"))
for (o in names(anot)) {
  a <- merge(anot[[o]], def[, c("COG_ID", "letras", "nombre", "via")], by = "COG_ID", all.x = TRUE)
  anot[[o]] <- a
  write.table(a, file.path(res_dir, paste0("00_anotacion_COG_", o, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

## términos por gen: categorías (multi-letra se reparte) y vías COG
terminos <- function(a, nivel) {
  if (nivel == "categoria") {
    l <- strsplit(a$letras, ""); data.frame(locus_tag = rep(a$locus_tag, lengths(l)), termino = catname[unlist(l)])
  } else {
    m <- merge(a[, c("locus_tag", "COG_ID")], pw, by = "COG_ID")
    data.frame(locus_tag = m$locus_tag, termino = m$via)
  }
}

ora <- function(genes, universo, t2g) {
  t2g <- t2g[t2g$locus_tag %in% universo, ]
  genes <- intersect(genes, unique(t2g$locus_tag)); U <- unique(t2g$locus_tag)
  N <- length(U); K <- length(genes)
  if (K < 3 || !nrow(t2g)) return(NULL)
  r <- do.call(rbind, lapply(split(t2g$locus_tag, t2g$termino), function(g) {
    g <- unique(g); n <- length(g); k <- sum(g %in% genes)
    data.frame(n_termino = n, k = k, esperado = n * K / N,
               p = phyper(k - 1, n, N - n, K, lower.tail = FALSE),
               genes = paste(sort(intersect(g, genes)), collapse = ",")) }))
  r <- r[r$n_termino >= MIN_TERM, ]; if (!nrow(r)) return(NULL); r$termino <- rownames(r)
  r$FDR <- p.adjust(r$p, "BH"); r$fold <- r$k / r$esperado; r$K_conjunto <- K; r$N_universo <- N
  r[order(r$p), c("termino", "k", "n_termino", "K_conjunto", "N_universo", "esperado", "fold", "p", "FDR", "genes")]
}

## ------------------------------------------------------- conjuntos -----------
sets <- list()
mE <- read.delim(file.path(wd, "results", "E109", "04_tabla_maestra_log2FC_padj.tsv"), quote = "")
mA <- read.delim(file.path(wd, "results", "Az39", "03_tabla_maestra_log2FC_padj.tsv"), quote = "")
de <- function(m, c1, s) { p <- m[[paste0("padj_", c1)]]; l <- m[[paste0("log2FC_", c1)]]
  m$locus_tag[!is.na(p) & p < ALPHA & abs(l) >= LFC & sign(l) == s] }
for (c1 in c("AIA_vs_Control", "Az39_vs_Control", "ipdC_vs_Control", "Az39_vs_ipdC", "Az39_vs_AIA"))
  for (s in c(1, -1)) sets[[paste0("E109 | ", c1, ifelse(s > 0, " ↑", " ↓"))]] <- list(org = "E109", genes = de(mE, c1, s))
cl <- read.delim(file.path(wd, "results", "E109", "02_clasificacion_respuesta_Az39.tsv"), quote = "")
sets[["E109 | Az39: AIA-dependiente (+parcial) ↑"]] <- list(org = "E109",
  genes = cl$locus_tag[cl$clase %in% c("AIA-dependiente", "Parcialmente AIA-dependiente") & cl$sentido_Az39 == "up"])
sets[["E109 | Az39: AIA-independiente ↑"]] <- list(org = "E109", genes = cl$locus_tag[cl$clase == "AIA-independiente" & cl$sentido_Az39 == "up"])
sets[["E109 | Az39: AIA-independiente ↓"]] <- list(org = "E109", genes = cl$locus_tag[cl$clase == "AIA-independiente" & cl$sentido_Az39 == "down"])
annA <- read.delim(file.path(wd, "data", "anotacion_Az39.tsv"), quote = "")
p3 <- annA$locus_tag[annA$replicon == "CP007796.1"]
for (c1 in c("WT_cocultivo_vs_sola", "ipdC_cocultivo_vs_sola", "WT_vs_ipdC_sola"))
  for (s in c(1, -1)) {
    g <- de(mA, c1, s); if (c1 == "WT_vs_ipdC_sola") g <- setdiff(g, p3)   # excluye plásmido p3 (ausente en ipdC)
    sets[[paste0("Az39 | ", c1, ifelse(c1 == "WT_vs_ipdC_sola", " (sin p3)", ""), ifelse(s > 0, " ↑", " ↓"))]] <- list(org = "Az39", genes = g)
  }
univ <- list(E109 = mE$locus_tag, Az39 = setdiff(mA$locus_tag, character(0)))

## ------------------------------------------------------------ ORA -----------
todo <- list()
for (nivel in c("categoria", "via")) {
  t2g <- lapply(anot, terminos, nivel = nivel)
  for (sn in names(sets)) {
    o <- sets[[sn]]$org
    u <- if (grepl("sin p3", sn)) setdiff(univ[[o]], p3) else univ[[o]]
    r <- ora(sets[[sn]]$genes, u, t2g[[o]])
    if (!is.null(r)) todo[[length(todo) + 1]] <- cbind(nivel = nivel, conjunto = sn, r)
  }
}
R <- do.call(rbind, todo)
write.table(R, file.path(res_dir, "01_ORA_COG_todos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sigR <- R[R$FDR < 0.05, ]
write.table(sigR, file.path(res_dir, "02_ORA_COG_significativos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cob <- sapply(names(anot), function(o) c(genes_con_COG = nrow(anot[[o]]),
  con_via = sum(anot[[o]]$COG_ID %in% pw$COG_ID), universo_testeado = length(univ[[o]]),
  universo_con_COG = sum(univ[[o]] %in% anot[[o]]$locus_tag)))
write.table(cbind(dato = rownames(cob), cob), file.path(res_dir, "00_cobertura_COG.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(cob)
cat(sprintf("términos significativos (FDR<0.05): %d en %d conjuntos\n", nrow(sigR), length(unique(sigR$conjunto))))
print(sigR[, c("nivel", "conjunto", "termino", "k", "n_termino", "fold", "FDR")], row.names = FALSE)

## ------------------------------------------------------------ dotplots ------
for (nivel in c("categoria", "via")) for (o in c("E109", "Az39")) {
  d <- R[R$nivel == nivel & startsWith(R$conjunto, o), ]
  terms <- unique(d$termino[d$FDR < 0.05]); if (!length(terms)) next
  d <- d[d$termino %in% terms, ]
  d$conjunto <- factor(sub("^\\w+ \\| ", "", d$conjunto), levels = unique(sub("^\\w+ \\| ", "", names(sets)[startsWith(names(sets), o)])))
  d$sig <- d$FDR < 0.05
  g <- ggplot(d, aes(conjunto, termino)) +
    geom_point(aes(size = k, colour = -log10(FDR), alpha = sig)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    scale_colour_gradient(low = "#9ECAE1", high = "#B2182B", name = "-log10 FDR") +
    scale_size_continuous(range = c(1, 7), name = "genes") +
    labs(x = NULL, y = NULL, title = sprintf("%s - enriquecimiento COG (%s)", o, ifelse(nivel == "categoria", "categorías funcionales", "vías")),
         subtitle = "puntos opacos = FDR < 0.05; ↑/↓ = genes que suben/bajan") +
    theme_bw() + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8), axis.text.y = element_text(size = 8))
  ggsave(file.path(fig_dir, sprintf("ORA_%s_%s.png", nivel, o)), g,
         width = 4 + 0.55 * nlevels(d$conjunto), height = 2 + 0.28 * length(terms), dpi = 200, limitsize = FALSE)
}
cat("listo\n")
