#!/usr/bin/env Rscript
## =============================================================================
## 06 - Integración transcriptoma (E109, Az39) + metaboloma de sobrenadantes
## Requiere haber corrido 01, 03 y 05.
##
## A) Integración por vías (condición): genes de cada vía (por anotación) con sus
##    log2FC en co-cultivo y los metabolitos de la vía con sus log2FC.
##    Enriquecimiento hipergeométrico de genes DE en la vía.
## B) Correlación gen-metabolito en muestras APAREADAS (mismo frasco, réplicas 1-3):
##      RNA 1-3 <-> E1-3 | 4-6 <-> Eaia1-3 | 7-9 <-> EA1-3   (lado E109, n=9)
##    Spearman, FDR BH. Exploratorio (n=9, 3 grupos: la correlación refleja sobre
##    todo diferencias entre tratamientos).
## C) DIABLO (mixOmics) lado E109: bloques RNA + metabolitos, 3 clases, n=9.
## Uso: Rscript 06_integracion_RNA_metabolitos.R <directorio_de_trabajo>
## =============================================================================
suppressMessages({ library(DESeq2); library(mixOmics); library(ggplot2); library(pheatmap) })
set.seed(1)
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args)) args[1] else "."
res_dir <- file.path(wd, "results", "integracion"); fig_dir <- file.path(wd, "figures", "integracion")
dir.create(res_dir, TRUE, TRUE); dir.create(fig_dir, TRUE, TRUE)
rd <- function(...) read.delim(file.path(wd, ...), quote = "", check.names = FALSE)

## ------------------------------------------------------------ insumos -------
annE <- rd("data", "anotacion_E109.tsv"); annE <- annE[!duplicated(annE$locus_tag), ]
annA <- rd("data", "anotacion_Az39.tsv"); annA <- annA[!duplicated(annA$locus_tag), ]
mE <- rd("results", "E109", "04_tabla_maestra_log2FC_padj.tsv")
mA <- rd("results", "Az39", "03_tabla_maestra_log2FC_padj.tsv")
met <- list(EA_vs_E = rd("results", "metabolomica", "DA_EA_vs_E.tsv"),
            Eaia_vs_E = rd("results", "metabolomica", "DA_Eaia_vs_E.tsv"),
            EA_vs_AzM = rd("results", "metabolomica", "DA_EA_vs_AzM.tsv"))

vias <- list(
  "Dicarboxilatos C4 (malato/fumarato)" = list(
    genes = "C4-dicarboxylate|fumarate hydratase|fumarate reductase|malate dehydrogenase|malic enzyme|succinate dehydrogenase",
    metab = c("Malato", "Ácido Fumarico")),
  "Triptófano / AIA" = list(
    genes = "tryptophan synthase|indole-3-glycerol|anthranilate|indole-3-pyruvate decarboxylase|indolepyruvate|tryptophan 2,3-dioxygenase|indole acetimide|kynurenin",
    metab = c("Triptofano")),
  "Osmoprotectores (trehalosa, betaína, prolina)" = list(
    genes = "trehalose|maltooligosyl|betaine|choline dehydrogenase|proline dehydrogenase|pyrroline-5-carboxylate|glutamate 5-kinase|gamma-glutamyl phosphate",
    metab = c("Trealosa", "Betaína", "Prolina")),
  "Glutamato / aspartato / asparagina" = list(
    genes = "glutamate synthase|glutamate dehydrogenase|glutamine synthetase|aspartate aminotransferase|asparaginase|asparagine synthase",
    metab = c("Ac. Glutamico", "Ac. Aspártico", "Asparagina")))
isla_iac <- sprintf("RN69_%05d", seq(30050, 30290, 5))

sig <- function(p, l) !is.na(p) & p < 0.05 & abs(l) >= 1
hyper <- function(inset, de) { N <- length(de); K <- sum(de); n <- sum(inset); k <- sum(inset & de)
  c(n = n, k = k, esperado = round(n * K / N, 1), p = phyper(k - 1, K, N - K, n, lower.tail = FALSE)) }

## ------------------------------------------------- A) vías por condición ----
filas <- list(); tablas_genes <- list()
for (v in names(vias)) {
  pat <- vias[[v]]$genes
  for (org in c("E109", "Az39")) {
    if (org == "E109") { m <- mE; ann <- annE; cc <- c("Az39_vs_Control", "AIA_vs_Control", "Az39_vs_ipdC") }
    else               { m <- mA; ann <- annA; cc <- c("WT_cocultivo_vs_sola", "WT_vs_ipdC_sola") }
    inset <- m$locus_tag %in% ann$locus_tag[grepl(pat, ann$product, ignore.case = TRUE)]
    if (org == "E109" && v == "Triptófano / AIA") inset <- inset | m$locus_tag %in% isla_iac
    for (c1 in cc) {
      de <- sig(m[[paste0("padj_", c1)]], m[[paste0("log2FC_", c1)]])
      h <- hyper(inset, de)
      up <- sum(inset & de & m[[paste0("log2FC_", c1)]] > 0); dn <- sum(inset & de & m[[paste0("log2FC_", c1)]] < 0)
      filas[[length(filas) + 1]] <- data.frame(via = v, organismo = org, contraste = c1, genes_via = h[["n"]],
        DE = h[["k"]], up = up, down = dn, esperado = h[["esperado"]], p_enriq = h[["p"]])
    }
    g <- m[inset, c("locus_tag", "product", grep("^(log2FC|padj)_", names(m), value = TRUE))]
    if (nrow(g)) tablas_genes[[paste(v, org)]] <- cbind(via = v, organismo = org, g)
  }
}
A <- do.call(rbind, filas); A$FDR_enriq <- p.adjust(A$p_enriq, "BH")
write.table(A, file.path(res_dir, "01_vias_genes_DE_enriquecimiento.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
G <- do.call(rbind, lapply(tablas_genes, function(x) { x[setdiff(c("via","organismo","locus_tag","product",
  grep("^(log2FC|padj)_", names(mE), value = TRUE), grep("^(log2FC|padj)_", names(mA), value = TRUE)), names(x))] <- NA
  x[, c("via","organismo","locus_tag","product", grep("^(log2FC|padj)_", names(mE), value = TRUE),
        grep("^(log2FC|padj)_", names(mA), value = TRUE))] }))
write.table(G, file.path(res_dir, "02_vias_genes_detalle.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

M <- do.call(rbind, lapply(names(vias), function(v) do.call(rbind, lapply(names(met), function(c1) {
  d <- met[[c1]]; d <- d[!is.na(d$metabolito) & d$metabolito %in% vias[[v]]$metab, ]
  if (!nrow(d)) return(NULL)
  data.frame(via = v, contraste = c1, senal = d$senal, metabolito = d$metabolito, log2FC = round(d$log2FC, 2), FDR = signif(d$FDR, 3)) }))))
write.table(M, file.path(res_dir, "03_vias_metabolitos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(A[A$contraste %in% c("Az39_vs_Control", "WT_cocultivo_vs_sola"), ])

## ------------------------------- B) correlación en muestras apareadas -------
cntE <- as.matrix(read.delim(file.path(wd, "data", "counts_E109.tsv"), row.names = 1, check.names = FALSE))
rna_s <- paste0(1:9, "-ARG"); met_s <- c(paste0("E", 1:3), paste0("Eaia", 1:3), paste0("EA", 1:3))
clase <- factor(rep(c("Control", "AIA", "Az39"), each = 3), levels = c("Control", "AIA", "Az39"))
cntE <- cntE[, paste0(1:12, "-ARG")]; cntE <- cntE[rowSums(cntE >= 10) >= 3, ]
dds <- DESeqDataSetFromMatrix(cntE, data.frame(c = factor(rep(1:4, each = 3)), row.names = colnames(cntE)), ~ c)
vs  <- assay(vst(dds, blind = TRUE))[, rna_s]

raw <- read.csv2(file.path(wd, "data", "metab_SN_matriz_asignada.csv"), header = FALSE, check.names = FALSE,
                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
Xm <- as.matrix(apply(raw[-(1:2), -1], 2, function(v) as.numeric(gsub(",", ".", v))))
rownames(Xm) <- trimws(raw[-(1:2), 1]); colnames(Xm) <- as.character(raw[1, -1])
Xm <- Xm[, met_s]; Xm <- Xm[apply(Xm, 1, function(v) sd(v) > 0 & mean(v > 1e-3) >= 0.25), ]
Lm <- log2(Xm + min(Xm[Xm > 0]) / 2)
asig <- !grepl("^U_", rownames(Lm))

genes_sel <- unique(c(isla_iac, unlist(lapply(vias, function(v) annE$locus_tag[grepl(v$genes, annE$product, ignore.case = TRUE)]))))
genes_sel <- intersect(genes_sel, rownames(vs))
cr <- expand.grid(gen = genes_sel, senal = rownames(Lm)[asig], stringsAsFactors = FALSE)
cr$rho <- NA_real_; cr$p <- NA_real_
for (i in seq_len(nrow(cr))) { t <- suppressWarnings(cor.test(vs[cr$gen[i], ], Lm[cr$senal[i], ], method = "spearman", exact = FALSE))
  cr$rho[i] <- t$estimate; cr$p[i] <- t$p.value }
cr$FDR <- p.adjust(cr$p, "BH")
cr$producto <- annE$product[match(cr$gen, annE$locus_tag)]
cr$en_isla_iac <- cr$gen %in% isla_iac
cr <- cr[order(cr$p), ]
write.table(cr, file.path(res_dir, "04_correlacion_gen_metabolito_apareadas_E109.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("pares gen-metabolito testeados: %d | FDR<0.05: %d | FDR<0.10: %d\n", nrow(cr), sum(cr$FDR < 0.05), sum(cr$FDR < 0.10)))

## heatmap: genes de la isla iac x metabolitos asignados
hm <- sapply(rownames(Lm)[asig], function(s) apply(vs[intersect(isla_iac, rownames(vs)), ], 1, function(g)
  cor(g, Lm[s, ], method = "spearman")))
pheatmap(hm, main = "Isla iac de E109 vs metabolitos del sobrenadante\n(Spearman, 9 muestras apareadas)",
         breaks = seq(-1, 1, length.out = 101), fontsize_row = 6, fontsize_col = 7,
         filename = file.path(fig_dir, "01_correlacion_isla_iac_vs_metabolitos.png"), width = 11, height = 8)

## -------------------------------------------- C) DIABLO lado E109 (n=9) -----
## bloques: RNA E109 + metabolitos del sobrenadante + metabolitos intracelulares (extractos)
rawx <- read.csv2(file.path(wd, "data", "metab_EXT_matriz_MCR.csv"), header = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
Xx <- as.matrix(apply(rawx[-(1:2), -1], 2, function(v) as.numeric(gsub(",", ".", v))))
rownames(Xx) <- trimws(rawx[-(1:2), 1]); colnames(Xx) <- as.character(rawx[1, -1])
Xx <- Xx[, met_s]; Xx <- Xx[apply(Xx, 1, function(v) sd(v) > 0 & mean(v > 1e-3) >= 0.25), ]
Lx <- log2(Xx + min(Xx[Xx > 0]) / 2)
top <- order(-apply(vs, 1, var))[1:1000]
X <- list(RNA_E109 = t(vs[top, ]), Metabolitos_SN = t(Lm), Metabolitos_EXT = t(Lx))
for (b in names(X)) rownames(X[[b]]) <- paste0(clase, "_", rep(1:3, 3))
design <- matrix(0.5, 3, 3); diag(design) <- 0
fit <- block.splsda(X, clase, ncomp = 2, design = design,
                    keepX = list(RNA_E109 = c(20, 20), Metabolitos_SN = c(10, 10), Metabolitos_EXT = c(10, 10)))
cc <- cor(sapply(names(X), function(b) fit$variates[[b]][, 1]))
cat("DIABLO: correlación entre componentes 1 de los bloques\n"); print(round(cc, 3))
write.table(round(cc, 3), file.path(res_dir, "06_DIABLO_correlacion_bloques.tsv"), sep = "\t", quote = FALSE)
fs <- do.call(rbind, lapply(names(X), function(b) { v <- selectVar(fit, block = b, comp = 1)[[b]]$value
  data.frame(bloque = b, variable = rownames(v), peso = v$value.var,
             producto = if (b == "RNA_E109") annE$product[match(rownames(v), annE$locus_tag)] else NA) }))
write.table(fs, file.path(res_dir, "05_DIABLO_firma_comp1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
png(file.path(fig_dir, "02_DIABLO_muestras.png"), width = 2400, height = 900, res = 200)
plotIndiv(fit, ind.names = TRUE, legend = TRUE, title = "DIABLO E109: RNA + metabolitos SN + extractos (n=9 apareadas)")
dev.off()
png(file.path(fig_dir, "03_DIABLO_circos.png"), width = 1800, height = 1800, res = 200)
try(circosPlot(fit, cutoff = 0.8, line = TRUE, size.variables = 0.5), silent = TRUE)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))
cat("listo\n")
