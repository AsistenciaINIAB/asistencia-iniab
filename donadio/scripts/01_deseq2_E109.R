#!/usr/bin/env Rscript
## =============================================================================
## 01 - Expresión diferencial de B. japonicum E109 (DESeq2)
## Pregunta: ¿qué parte de la respuesta de E109 a A. argentinense Az39
##           depende de la producción de AIA por Az39? (mutante Az39-ipdC)
##
## Uso:  Rscript 01_deseq2_E109.R <directorio_de_trabajo>
##   espera  <wd>/data/counts_E109.tsv      (locus_tag x muestras, conteos crudos)
##           <wd>/data/anotacion_E109.tsv   (feature table NCBI GCA_000807315.1)
##   escribe <wd>/results/E109/*.tsv  y  <wd>/figures/E109/*.pdf|png
## Ambiente: donadio-multiomica (envs/donadio-multiomica.yml)
## =============================================================================
suppressMessages({
  library(DESeq2); library(apeglm); library(ggplot2)
  library(pheatmap); library(ggrepel)
})
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
wd   <- if (length(args)) args[1] else "."
res_dir <- file.path(wd, "results", "E109"); fig_dir <- file.path(wd, "figures", "E109")
dir.create(res_dir, TRUE, TRUE); dir.create(fig_dir, TRUE, TRUE)

ALPHA <- 0.05   # FDR
LFC   <- 1      # |log2FC| para considerar un cambio "relevante"

## ---------------------------------------------------------------- datos -----
cnt <- as.matrix(read.delim(file.path(wd, "data", "counts_E109.tsv"),
                            row.names = 1, check.names = FALSE))
ann <- read.delim(file.path(wd, "data", "anotacion_E109.tsv"), quote = "")
ann <- ann[!duplicated(ann$locus_tag), ]
rownames(ann) <- ann$locus_tag

## metadatos (Metadata Transcriptomica.docx; 4 h de interacción, 3 réplicas)
meta <- data.frame(
  sample = paste0(1:12, "-ARG"),
  cond   = rep(c("Control", "AIA", "Az39", "Az39_ipdC"), each = 3),
  rep    = rep(1:3, 4))
meta$cond <- factor(meta$cond, levels = c("Control", "AIA", "Az39", "Az39_ipdC"))
rownames(meta) <- meta$sample
cnt <- cnt[, meta$sample]

## ------------------------------------------------------------------- QC -----
qc <- data.frame(sample = colnames(cnt), cond = meta$cond,
                 lecturas_asignadas = colSums(cnt),
                 genes_detectados = colSums(cnt > 0))
write.table(qc, file.path(res_dir, "00_QC_librerias.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## filtro: >= 10 lecturas en al menos 3 muestras (tamaño mínimo de grupo)
keep <- rowSums(cnt >= 10) >= 3
dds  <- DESeqDataSetFromMatrix(cnt[keep, ], meta, ~ cond)
dds  <- DESeq(dds, quiet = TRUE)
cat(sprintf("genes de entrada: %d | tras filtro: %d\n", nrow(cnt), sum(keep)))

vsd <- vst(dds, blind = TRUE)
## PCA
pc  <- prcomp(t(assay(vsd)[order(-apply(assay(vsd), 1, var))[1:1000], ]))
ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pdat <- data.frame(pc$x[, 1:3], meta)
cols <- c(Control = "#4D4D4D", AIA = "#E69F00", Az39 = "#0072B2", Az39_ipdC = "#CC79A7")
g <- ggplot(pdat, aes(PC1, PC2, colour = cond, label = sample)) +
  geom_point(size = 3.5) + geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = "Tratamiento") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "E109 - PCA (VST, 1000 genes más variables)") + theme_bw()
ggsave(file.path(fig_dir, "01_PCA_E109.pdf"), g, width = 6.5, height = 5)
ggsave(file.path(fig_dir, "01_PCA_E109.png"), g, width = 6.5, height = 5, dpi = 200)
## correlación entre muestras
cm <- cor(assay(vsd), method = "pearson")
ann_col <- data.frame(Tratamiento = meta$cond, row.names = meta$sample)
pheatmap(cm, annotation_col = ann_col, annotation_row = ann_col,
         annotation_colors = list(Tratamiento = cols),
         main = "E109 - correlación entre muestras (VST)", display_numbers = TRUE,
         number_format = "%.3f", fontsize_number = 6,
         filename = file.path(fig_dir, "02_correlacion_muestras_E109.pdf"),
         width = 8, height = 7)

## --------------------------------------------------------- contrastes -------
## numerador vs denominador ; log2FC > 0 = mayor en el numerador
contr <- list(
  AIA_vs_Control       = c("AIA", "Control"),
  Az39_vs_Control      = c("Az39", "Control"),
  ipdC_vs_Control      = c("Az39_ipdC", "Control"),
  Az39_vs_ipdC         = c("Az39", "Az39_ipdC"),   # efecto del AIA de Az39
  Az39_vs_AIA          = c("Az39", "AIA"),         # lo que Az39 aporta más allá del AIA
  ipdC_vs_AIA          = c("Az39_ipdC", "AIA"))

run_contrast <- function(num, den) {
  d2 <- dds; d2$cond <- relevel(d2$cond, ref = den)
  d2 <- nbinomWaldTest(estimateDispersions(estimateSizeFactors(d2), quiet = TRUE), quiet = TRUE)
  coef <- paste0("cond_", num, "_vs_", den)
  r  <- results(d2, name = coef, alpha = ALPHA)
  sh <- lfcShrink(d2, coef = coef, type = "apeglm", quiet = TRUE)
  data.frame(locus_tag = rownames(r), baseMean = r$baseMean,
             log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
             log2FC_shrunk = sh$log2FoldChange,
             pvalue = r$pvalue, padj = r$padj)
}
R <- lapply(contr, function(x) run_contrast(x[1], x[2]))

addann <- function(df) cbind(df, ann[df$locus_tag, c("gene_symbol", "product", "start", "end", "strand")])
resumen <- data.frame()
for (n in names(R)) {
  df <- addann(R[[n]]); df <- df[order(df$padj), ]
  write.table(df, file.path(res_dir, paste0("DE_", n, ".tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
  s <- !is.na(df$padj) & df$padj < ALPHA
  resumen <- rbind(resumen, data.frame(
    contraste = n,
    DE_FDR005 = sum(s),
    up = sum(s & df$log2FC > 0), down = sum(s & df$log2FC < 0),
    DE_FDR005_absLFC1 = sum(s & abs(df$log2FC) >= LFC),
    up_LFC1 = sum(s & df$log2FC >= LFC), down_LFC1 = sum(s & df$log2FC <= -LFC)))
}
write.table(resumen, file.path(res_dir, "01_resumen_contrastes.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(resumen)

## ---------------------- clasificación: ¿la respuesta a Az39 depende del AIA? --
## Gen "responde a Az39": padj<0.05 y |log2FC|>=1 en Az39_vs_Control.
## - AIA-dependiente : NO responde a Az39-ipdC (padj>=0.05 o sentido opuesto)
##                     Y difiere entre Az39 y Az39-ipdC (padj<0.05, mismo sentido que Az39)
## - AIA-independiente: responde igual a Az39-ipdC (padj<0.05, mismo sentido)
##                     Y NO difiere entre Az39 y Az39-ipdC (padj>=0.05)
## - Parcial         : responde a ambos pero con diferencia significativa entre ellos
## - Ambiguo         : el resto (poca potencia para decidir)
## Además: ¿lo reproduce el AIA exógeno? (padj<0.05 en AIA_vs_Control, mismo sentido)
m <- Reduce(function(a, b) merge(a, b, by = "locus_tag"), lapply(names(R), function(n) {
  x <- R[[n]][, c("locus_tag", "log2FC", "padj")]
  names(x)[2:3] <- paste0(c("log2FC_", "padj_"), n); x }))
sig <- function(p) !is.na(p) & p < ALPHA
sgn <- function(x) sign(x)
resp_az <- sig(m$padj_Az39_vs_Control) & abs(m$log2FC_Az39_vs_Control) >= LFC
same_ipdC <- sig(m$padj_ipdC_vs_Control) & sgn(m$log2FC_ipdC_vs_Control) == sgn(m$log2FC_Az39_vs_Control)
diff_wt_mut <- sig(m$padj_Az39_vs_ipdC) & sgn(m$log2FC_Az39_vs_ipdC) == sgn(m$log2FC_Az39_vs_Control)
m$clase <- NA_character_
m$clase[resp_az & !same_ipdC & diff_wt_mut]  <- "AIA-dependiente"
m$clase[resp_az & same_ipdC & !sig(m$padj_Az39_vs_ipdC)] <- "AIA-independiente"
m$clase[resp_az & same_ipdC & diff_wt_mut]   <- "Parcialmente AIA-dependiente"
m$clase[resp_az & is.na(m$clase)]            <- "Ambiguo"
m$reproducido_por_AIA_exogeno <- ifelse(resp_az,
  sig(m$padj_AIA_vs_Control) & sgn(m$log2FC_AIA_vs_Control) == sgn(m$log2FC_Az39_vs_Control), NA)
m$sentido_Az39 <- ifelse(resp_az, ifelse(m$log2FC_Az39_vs_Control > 0, "up", "down"), NA)
cl <- addann(m[resp_az, ])
cl <- cl[order(cl$clase, -abs(cl$log2FC_Az39_vs_Control)), ]
write.table(cl, file.path(res_dir, "02_clasificacion_respuesta_Az39.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
tab <- as.data.frame.matrix(table(cl$clase, cl$sentido_Az39))
tab$reproducido_por_AIA <- tapply(cl$reproducido_por_AIA_exogeno, cl$clase, sum)[rownames(tab)]
write.table(cbind(clase = rownames(tab), tab), file.path(res_dir, "03_resumen_clasificacion.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(tab)
write.table(addann(m), file.path(res_dir, "04_tabla_maestra_log2FC_padj.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## ------------------------------------------------- operón iac (RN69_30270-30305)
iac <- sprintf("RN69_%05d", seq(30270, 30305, 5))
io  <- addann(m[m$locus_tag %in% iac, ])
write.table(io, file.path(res_dir, "05_operon_iac.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
long <- do.call(rbind, lapply(c("AIA_vs_Control", "Az39_vs_Control", "ipdC_vs_Control"), function(n)
  data.frame(locus_tag = io$locus_tag, producto = io$product, contraste = n,
             log2FC = io[[paste0("log2FC_", n)]], padj = io[[paste0("padj_", n)]])))
long$contraste <- factor(long$contraste, levels = c("AIA_vs_Control", "Az39_vs_Control", "ipdC_vs_Control"),
                         labels = c("+ AIA", "+ Az39", "+ Az39-ipdC"))
long$etq <- ifelse(is.na(long$padj), "", ifelse(long$padj < 0.001, "***", ifelse(long$padj < 0.01, "**",
                    ifelse(long$padj < 0.05, "*", ""))))
long$gen <- paste0(long$locus_tag, "\n", substr(long$producto, 1, 28))
g <- ggplot(long, aes(gen, log2FC, fill = contraste)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = etq, y = log2FC + ifelse(log2FC >= 0, 0.3, -0.6)),
            position = position_dodge(0.8), size = 3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("+ AIA" = "#E69F00", "+ Az39" = "#0072B2", "+ Az39-ipdC" = "#CC79A7"),
                    name = "E109 vs control") +
  labs(x = NULL, y = "log2FC vs E109 control",
       title = "Cluster iac (catabolismo de AIA) de E109",
       subtitle = "* padj<0.05  ** <0.01  *** <0.001") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(fig_dir, "03_operon_iac.pdf"), g, width = 9, height = 5)
ggsave(file.path(fig_dir, "03_operon_iac.png"), g, width = 9, height = 5, dpi = 200)

## ------------------------------------------------------------- volcanos -----
for (n in names(R)) {
  df <- addann(R[[n]]); df <- df[!is.na(df$padj), ]
  df$grupo <- ifelse(df$padj < ALPHA & df$log2FC >= LFC, "up",
              ifelse(df$padj < ALPHA & df$log2FC <= -LFC, "down", "ns"))
  top <- head(df[df$grupo != "ns", ][order(df$padj[df$grupo != "ns"]), ], 12)
  top$lab <- ifelse(is.na(top$gene_symbol) | top$gene_symbol == "", top$locus_tag, top$gene_symbol)
  g <- ggplot(df, aes(log2FC, -log10(pvalue), colour = grupo)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_text_repel(data = top, aes(label = lab), size = 2.6, max.overlaps = 30, show.legend = FALSE) +
    scale_colour_manual(values = c(up = "#B2182B", down = "#2166AC", ns = "grey70")) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = 2, linewidth = 0.3) +
    labs(title = gsub("_", " ", n), subtitle = sprintf("up=%d  down=%d  (padj<%.2f, |log2FC|>=%g)",
         sum(df$grupo == "up"), sum(df$grupo == "down"), ALPHA, LFC)) + theme_bw()
  ggsave(file.path(fig_dir, paste0("04_volcano_", n, ".png")), g, width = 6, height = 5, dpi = 200)
}

## comparación con los resultados del proveedor (Exp.3: E109 vs E109+Az39) si está
prov <- file.path(wd, "data", "proveedor_Exp3_DESeq2.csv")
if (file.exists(prov)) {
  p <- read.csv2(prov, row.names = 1)
  a <- R$Az39_vs_Control; rownames(a) <- a$locus_tag
  cg <- intersect(rownames(p), rownames(a))
  ok <- is.finite(p[cg, "log2FoldChange"]) & is.finite(a[cg, "log2FC"])
  cat(sprintf("correlación log2FC con DESeq2 del proveedor (Exp.3): r = %.3f (n=%d)\n",
              cor(p[cg, "log2FoldChange"][ok], a[cg, "log2FC"][ok]), sum(ok)))
}
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))
cat("listo:", res_dir, "\n")
