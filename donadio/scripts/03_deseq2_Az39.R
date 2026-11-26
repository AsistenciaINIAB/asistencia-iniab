#!/usr/bin/env Rscript
## =============================================================================
## 03 - Expresión diferencial de A. argentinense Az39 (DESeq2)
## Diseño factorial 2x2: genotipo (WT / ipdC) x condición (sola / co-cultivo con E109)
##   7-9   Az39 WT   + E109        13-15  Az39 WT   sola (MMAB)
##   10-12 Az39 ipdC + E109        16-18  Az39 ipdC sola (MMAB)
## Preguntas: (1) ¿cómo responde Az39 a E109?  (2) ¿qué cambia en el mutante ipdC?
##            (3) ¿la respuesta a E109 depende del genotipo? (interacción)
## Uso: Rscript 03_deseq2_Az39.R <directorio_de_trabajo>
## =============================================================================
suppressMessages({ library(DESeq2); library(apeglm); library(ggplot2)
                   library(pheatmap); library(ggrepel) })
set.seed(1)
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args)) args[1] else "."
res_dir <- file.path(wd, "results", "Az39"); fig_dir <- file.path(wd, "figures", "Az39")
dir.create(res_dir, TRUE, TRUE); dir.create(fig_dir, TRUE, TRUE)
ALPHA <- 0.05; LFC <- 1

cnt <- as.matrix(read.delim(file.path(wd, "data", "counts_Az39.tsv"), row.names = 1, check.names = FALSE))
ann <- read.delim(file.path(wd, "data", "anotacion_Az39.tsv"), quote = "")
ann <- ann[!duplicated(ann$locus_tag), ]; rownames(ann) <- ann$locus_tag
rep_names <- c(CP007793.1 = "cromosoma", CP007794.1 = "p1", CP007795.1 = "p2",
               CP007796.1 = "p3", CP007797.1 = "p4", CP007798.1 = "p5")
ann$replicon_nombre <- rep_names[ann$replicon]

meta <- data.frame(sample = paste0(7:18, "-ARG"),
  genotipo  = rep(c("WT", "ipdC", "WT", "ipdC"), each = 3),
  condicion = rep(c("cocultivo", "cocultivo", "sola", "sola"), each = 3))
meta$genotipo  <- factor(meta$genotipo,  levels = c("WT", "ipdC"))
meta$condicion <- factor(meta$condicion, levels = c("sola", "cocultivo"))
meta$grupo <- factor(paste(meta$genotipo, meta$condicion, sep = "_"),
                     levels = c("WT_sola", "WT_cocultivo", "ipdC_sola", "ipdC_cocultivo"))
rownames(meta) <- meta$sample
cnt <- cnt[, meta$sample]

qc <- data.frame(sample = meta$sample, grupo = meta$grupo, lecturas_asignadas = colSums(cnt),
                 genes_detectados = colSums(cnt > 0))
write.table(qc, file.path(res_dir, "00_QC_librerias.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ipdC (ABAZ39_11295): validación del mutante a nivel de transcrito
cpm <- t(t(cnt) / colSums(cnt)) * 1e6
ip <- data.frame(meta[, c("sample", "grupo")], ipdC_conteos = cnt["ABAZ39_11295", ],
                 ipdC_CPM = round(cpm["ABAZ39_11295", ], 1))
write.table(ip, file.path(res_dir, "01_expresion_ipdC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

keep <- rowSums(cnt >= 10) >= 3
dds  <- DESeqDataSetFromMatrix(cnt[keep, ], meta, ~ grupo)
dds  <- DESeq(dds, quiet = TRUE)
cat(sprintf("genes de entrada: %d | tras filtro: %d\n", nrow(cnt), sum(keep)))

vsd <- vst(dds, blind = TRUE)
pc  <- prcomp(t(assay(vsd)[order(-apply(assay(vsd), 1, var))[1:1000], ]))
ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pdat <- data.frame(pc$x[, 1:3], meta)
cols <- c(WT_sola = "#0072B2", WT_cocultivo = "#56B4E9", ipdC_sola = "#CC79A7", ipdC_cocultivo = "#F0A3C8")
g <- ggplot(pdat, aes(PC1, PC2, colour = grupo, shape = condicion, label = sample)) +
  geom_point(size = 3.5) + geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = "Grupo") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "Az39 - PCA (VST, 1000 genes más variables)") + theme_bw()
ggsave(file.path(fig_dir, "01_PCA_Az39.png"), g, width = 6.5, height = 5, dpi = 200)
ggsave(file.path(fig_dir, "01_PCA_Az39.pdf"), g, width = 6.5, height = 5)
ann_col <- data.frame(Grupo = meta$grupo, row.names = meta$sample)
pheatmap(cor(assay(vsd)), annotation_col = ann_col, annotation_row = ann_col,
         annotation_colors = list(Grupo = cols), display_numbers = TRUE, number_format = "%.3f",
         fontsize_number = 6, main = "Az39 - correlación entre muestras (VST)",
         filename = file.path(fig_dir, "02_correlacion_muestras_Az39.pdf"), width = 8, height = 7)

## contrastes (numerador vs denominador)
contr <- list(
  WT_cocultivo_vs_sola   = c("WT_cocultivo", "WT_sola"),       # respuesta de Az39 WT a E109
  ipdC_cocultivo_vs_sola = c("ipdC_cocultivo", "ipdC_sola"),   # respuesta del mutante a E109
  WT_vs_ipdC_sola        = c("WT_sola", "ipdC_sola"),          # efecto del mutante, sin E109
  WT_vs_ipdC_cocultivo   = c("WT_cocultivo", "ipdC_cocultivo"))# efecto del mutante, con E109
run_contrast <- function(num, den) {
  d2 <- dds; d2$grupo <- relevel(d2$grupo, ref = den)
  d2 <- nbinomWaldTest(estimateDispersions(estimateSizeFactors(d2), quiet = TRUE), quiet = TRUE)
  coef <- paste0("grupo_", num, "_vs_", den)
  r <- results(d2, name = coef, alpha = ALPHA)
  sh <- lfcShrink(d2, coef = coef, type = "apeglm", quiet = TRUE)
  data.frame(locus_tag = rownames(r), baseMean = r$baseMean, log2FC = r$log2FoldChange,
             lfcSE = r$lfcSE, log2FC_shrunk = sh$log2FoldChange, pvalue = r$pvalue, padj = r$padj)
}
R <- lapply(contr, function(x) run_contrast(x[1], x[2]))

## interacción: ¿la respuesta a E109 difiere entre WT e ipdC?
dds_i <- DESeqDataSetFromMatrix(cnt[keep, ], meta, ~ genotipo + condicion + genotipo:condicion)
dds_i <- DESeq(dds_i, quiet = TRUE)
ri <- results(dds_i, name = "genotipoipdC.condicioncocultivo", alpha = ALPHA)
R$interaccion_genotipo_x_cocultivo <- data.frame(locus_tag = rownames(ri), baseMean = ri$baseMean,
  log2FC = ri$log2FoldChange, lfcSE = ri$lfcSE, log2FC_shrunk = NA, pvalue = ri$pvalue, padj = ri$padj)

addann <- function(df) cbind(df, ann[df$locus_tag, c("gene_symbol", "product", "replicon_nombre", "start", "end", "strand")])
resumen <- data.frame()
for (n in names(R)) {
  df <- addann(R[[n]]); df <- df[order(df$padj), ]
  write.table(df, file.path(res_dir, paste0("DE_", n, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  s <- !is.na(df$padj) & df$padj < ALPHA
  resumen <- rbind(resumen, data.frame(contraste = n, DE_FDR005 = sum(s), up = sum(s & df$log2FC > 0),
    down = sum(s & df$log2FC < 0), DE_absLFC1 = sum(s & abs(df$log2FC) >= LFC),
    up_LFC1 = sum(s & df$log2FC >= LFC), down_LFC1 = sum(s & df$log2FC <= -LFC)))
}
write.table(resumen, file.path(res_dir, "02_resumen_contrastes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(resumen)

## tabla maestra
m <- Reduce(function(a, b) merge(a, b, by = "locus_tag"), lapply(names(R), function(n) {
  x <- R[[n]][, c("locus_tag", "log2FC", "padj")]; names(x)[2:3] <- paste0(c("log2FC_", "padj_"), n); x }))
write.table(addann(m), file.path(res_dir, "03_tabla_maestra_log2FC_padj.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## distribución por replicón de los genes que responden a E109 (WT) - test hipergeométrico
de <- R$WT_cocultivo_vs_sola; de$rep <- ann[de$locus_tag, "replicon_nombre"]
de$sig <- !is.na(de$padj) & de$padj < ALPHA & abs(de$log2FC) >= LFC
N <- nrow(de); K <- sum(de$sig)
rp <- do.call(rbind, lapply(split(de, de$rep), function(x) {
  n <- nrow(x); k <- sum(x$sig)
  data.frame(replicon = x$rep[1], genes_testeados = n, DE = k, esperado = round(n * K / N, 1),
             up = sum(x$sig & x$log2FC > 0), down = sum(x$sig & x$log2FC < 0),
             p_enriq = phyper(k - 1, K, N - K, n, lower.tail = FALSE),
             p_empob = phyper(k, K, N - K, n)) }))
rp$FDR_enriq <- p.adjust(rp$p_enriq, "BH"); rp$FDR_empob <- p.adjust(rp$p_empob, "BH")
write.table(rp, file.path(res_dir, "04_replicones_respuesta_E109.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(rp)

## volcanos
for (n in names(R)) {
  df <- addann(R[[n]]); df <- df[!is.na(df$padj), ]
  df$grupo <- ifelse(df$padj < ALPHA & df$log2FC >= LFC, "up", ifelse(df$padj < ALPHA & df$log2FC <= -LFC, "down", "ns"))
  top <- head(df[df$grupo != "ns", ][order(df$padj[df$grupo != "ns"]), ], 12)
  top$lab <- ifelse(is.na(top$gene_symbol) | top$gene_symbol == "", top$locus_tag, top$gene_symbol)
  g <- ggplot(df, aes(log2FC, -log10(pvalue), colour = grupo)) + geom_point(size = 0.8, alpha = 0.6) +
    geom_text_repel(data = top, aes(label = lab), size = 2.6, max.overlaps = 30, show.legend = FALSE) +
    scale_colour_manual(values = c(up = "#B2182B", down = "#2166AC", ns = "grey70")) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = 2, linewidth = 0.3) +
    labs(title = gsub("_", " ", n), subtitle = sprintf("up=%d  down=%d  (padj<%.2f, |log2FC|>=%g)",
         sum(df$grupo == "up"), sum(df$grupo == "down"), ALPHA, LFC)) + theme_bw()
  ggsave(file.path(fig_dir, paste0("03_volcano_", n, ".png")), g, width = 6, height = 5, dpi = 200)
}
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))
cat("listo:", res_dir, "\n")
