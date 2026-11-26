#!/usr/bin/env Rscript
## =============================================================================
## 07 - Gráficos adicionales de DIABLO (lado E109, 9 muestras apareadas)
## Reajusta el mismo modelo que 06 (RNA E109 + metabolitos SN + extractos;
## 3 clases; ncomp=2; keepX RNA 20/20, SN 10/10, EXT 10/10; design 0.5) y genera:
##   01 plotDiablo      - correlación entre bloques en la componente 1
##   02 plotArrow       - consenso de las muestras entre bloques
##   03 plotIndiv       - muestras por bloque (comp 1 vs 2)
##   04 plotLoadings    - pesos de las variables seleccionadas (comp 1 y 2)
##   05 plotVar         - círculo de correlación de las variables seleccionadas
##   06 cimDiablo       - heatmap de las variables seleccionadas
##   07 red             - red de correlaciones entre bloques (|r| >= 0.9)
##   08 perfiles        - valores por muestra de las variables top de cada bloque
## Uso: Rscript 07_diablo_graficos.R <directorio_de_trabajo>
## =============================================================================
suppressMessages({ library(DESeq2); library(mixOmics); library(ggplot2) })
set.seed(1)
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args)) args[1] else "."
fig <- file.path(wd, "figures", "integracion", "diablo"); dir.create(fig, TRUE, TRUE)
res <- file.path(wd, "results", "integracion"); dir.create(res, TRUE, TRUE)

annE <- read.delim(file.path(wd, "data", "anotacion_E109.tsv"), quote = "")
annE <- annE[!duplicated(annE$locus_tag), ]
rna_s <- paste0(1:9, "-ARG")
met_s <- c(paste0("E", 1:3), paste0("Eaia", 1:3), paste0("EA", 1:3))
clase <- factor(rep(c("Control", "AIA", "Az39"), each = 3), levels = c("Control", "AIA", "Az39"))
ids <- paste0(clase, "_", rep(1:3, 3))
cols <- c(Control = "#4D4D4D", AIA = "#E69F00", Az39 = "#0072B2")

## --- bloques (igual que en 06) ---
cnt <- as.matrix(read.delim(file.path(wd, "data", "counts_E109.tsv"), row.names = 1, check.names = FALSE))
cnt <- cnt[, paste0(1:12, "-ARG")]; cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
dds <- DESeqDataSetFromMatrix(cnt, data.frame(c = factor(rep(1:4, each = 3)), row.names = colnames(cnt)), ~ c)
vs  <- assay(vst(dds, blind = TRUE))[, rna_s]
leer <- function(f) {
  r <- read.csv2(file.path(wd, "data", f), header = FALSE, check.names = FALSE,
                 fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  X <- as.matrix(apply(r[-(1:2), -1], 2, function(v) as.numeric(gsub(",", ".", v))))
  rownames(X) <- make.unique(trimws(r[-(1:2), 1])); colnames(X) <- as.character(r[1, -1])
  X <- X[, met_s]; X <- X[apply(X, 1, function(v) sd(v) > 0 & mean(v > 1e-3) >= 0.25), ]
  log2(X + min(X[X > 0]) / 2)
}
Lm <- leer("metab_SN_matriz_asignada.csv"); Lx <- leer("metab_EXT_matriz_MCR.csv")
top <- order(-apply(vs, 1, var))[1:1000]
R <- t(vs[top, ])
## etiquetas legibles para genes: locus + producto abreviado
prod <- annE$product[match(colnames(R), annE$locus_tag)]
colnames(R) <- ifelse(is.na(prod), colnames(R), paste0(colnames(R), " ", substr(prod, 1, 28)))
X <- list(RNA = R, Sobrenadante = t(Lm), Extractos = t(Lx))
for (b in names(X)) rownames(X[[b]]) <- ids

design <- matrix(0.5, 3, 3, dimnames = list(names(X), names(X))); diag(design) <- 0
fit <- block.splsda(X, clase, ncomp = 2, design = design,
                    keepX = list(RNA = c(20, 20), Sobrenadante = c(10, 10), Extractos = c(10, 10)))

png2 <- function(f, w = 2000, h = 1600) png(file.path(fig, f), width = w, height = h, res = 220)

## 01 correlación entre bloques (comp 1): gráfico de pares propio (plotDiablo dibuja r con fuente gigante)
suppressMessages(library(ggrepel))
V <- sapply(names(X), function(b) fit$variates[[b]][, 1])
pp <- do.call(rbind, combn(names(X), 2, simplify = FALSE, FUN = function(p)
  data.frame(par = sprintf("%s vs %s  (r = %.3f)", p[1], p[2], cor(V[, p[1]], V[, p[2]])),
             x = V[, p[1]], y = V[, p[2]], grupo = clase, muestra = ids)))
g <- ggplot(pp, aes(x, y, colour = grupo)) + geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey50", linewidth = 0.4, formula = y ~ x) +
  geom_text_repel(aes(label = muestra), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) + facet_wrap(~ par, scales = "free") +
  labs(x = "componente 1 (bloque 1)", y = "componente 1 (bloque 2)",
       title = "DIABLO - acuerdo entre bloques en la componente 1") + theme_bw()
ggsave(file.path(fig, "01_correlacion_bloques_comp1.png"), g, width = 12, height = 4.3, dpi = 200)

## 02 flechas: cada muestra en cada bloque hacia el centroide de consenso
png2("02_plotArrow_consenso_muestras.png", 2000, 1600)
print(plotArrow(fit, ind.names = TRUE, legend = TRUE, col.per.group = cols,
                title = "DIABLO - consenso entre bloques (punta = centroide)"))
dev.off()

## 03 muestras por bloque
png2("03_plotIndiv_por_bloque.png", 3000, 1100)
plotIndiv(fit, ind.names = FALSE, legend = TRUE, col.per.group = cols, ellipse = FALSE,
          title = "DIABLO - muestras por bloque (comp 1 vs 2)", pch = 16, cex = 3)
dev.off()

## 04 pesos de variables (comp 1 y 2)
for (k in 1:2) {
  png2(sprintf("04_plotLoadings_comp%d.png", k), 3200, 1700)
  plotLoadings(fit, comp = k, contrib = "max", method = "median", legend.color = cols,
               size.name = 0.55, size.title = 1, title = sprintf("Pesos - componente %d (color = grupo con valor máximo)", k))
  dev.off()
}

## 05 círculo de correlación (propio, con etiquetas que no se superponen)
cv <- do.call(rbind, lapply(names(X), function(b) {
  vv <- unique(c(selectVar(fit, block = b, comp = 1)[[b]]$name, selectVar(fit, block = b, comp = 2)[[b]]$name))
  data.frame(bloque = b, variable = vv,
             c1 = as.numeric(cor(X[[b]][, vv, drop = FALSE], fit$variates[[b]][, 1])),
             c2 = as.numeric(cor(X[[b]][, vv, drop = FALSE], fit$variates[[b]][, 2]))) }))
cv$etq <- ifelse(cv$bloque == "RNA", sub("^(RN69_\\d+) ", "\\1\n", cv$variable), cv$variable)
circ <- function(r) data.frame(x = r * cos(seq(0, 2 * pi, length.out = 200)), y = r * sin(seq(0, 2 * pi, length.out = 200)))
bcol <- c(RNA = "#0072B2", Sobrenadante = "#D55E00", Extractos = "#009E73")
g <- ggplot(cv, aes(c1, c2, colour = bloque)) +
  geom_path(data = circ(1), aes(x, y), inherit.aes = FALSE, colour = "grey60") +
  geom_path(data = circ(0.5), aes(x, y), inherit.aes = FALSE, colour = "grey80", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 2) + geom_text_repel(aes(label = etq), size = 2.2, max.overlaps = 60, lineheight = 0.8,
                                         segment.size = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = bcol, name = "Bloque") + coord_equal() +
  labs(x = "correlación con componente 1", y = "correlación con componente 2",
       title = "DIABLO - círculo de correlación de las variables seleccionadas",
       subtitle = "comp 1: separa Az39 de Control/AIA | comp 2: separa AIA de Control") + theme_bw()
ggsave(file.path(fig, "05_circulo_correlacion.png"), g, width = 11, height = 10, dpi = 200)

## 06 heatmap de variables seleccionadas (cimDiablo ordena la leyenda alfabéticamente:
##    se reajusta con niveles alfabéticos para que colores y nombres coincidan)
fit_cim <- fit; lv <- sort(levels(clase))
fit_cim$Y <- factor(as.character(clase), levels = lv)
png2("06_cimDiablo_heatmap.png", 2600, 2600)
cimDiablo(fit_cim, color.blocks = c("#0072B2", "#D55E00", "#009E73"), color.Y = cols[lv],
          comp = 1:2, margins = c(14, 22), size.legend = 0.9, legend.position = "topright",
          trim = 3, transpose = FALSE)
dev.off()

## 07 red de correlaciones entre variables de distintos bloques
ok <- try({
  png2("07_red_correlaciones.png", 2600, 2600)
  network(fit, blocks = 1:3, cutoff = 0.9, color.node = c("#9ECAE1", "#FDAE6B", "#A1D99B"),
          shape.node = c("rectangle", "circle", "circle"), cex.node.name = 0.55,
          size.node = 0.12, graph.scale = 0.6, plot.graph = TRUE)
  dev.off() }, silent = TRUE)
if (inherits(ok, "try-error")) { try(dev.off(), silent = TRUE); cat("network() no disponible:", ok, "\n") }

## 08 perfiles por muestra de las 6 variables con más peso por bloque
perf <- do.call(rbind, lapply(names(X), function(b) {
  v <- selectVar(fit, block = b, comp = 1)[[b]]$value
  sel <- head(rownames(v)[order(-abs(v$value.var))], 6)
  do.call(rbind, lapply(sel, function(s) data.frame(bloque = b, variable = s, muestra = ids,
    grupo = clase, valor = as.numeric(scale(X[[b]][, s]))))) }))
perf$variable <- factor(perf$variable, levels = unique(perf$variable))
g <- ggplot(perf, aes(grupo, valor, colour = grupo)) +
  geom_point(size = 2.5, position = position_jitter(width = 0.08, seed = 1)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "black", linewidth = 0.3) +
  scale_colour_manual(values = cols, guide = "none") +
  facet_wrap(~ bloque + variable, scales = "free_y", ncol = 6, labeller = label_wrap_gen(24)) +
  labs(x = NULL, y = "valor estandarizado (z)", title = "Variables con mayor peso en la componente 1 de DIABLO, por bloque") +
  theme_bw() + theme(strip.text = element_text(size = 6.5), axis.text.x = element_text(size = 7))
ggsave(file.path(fig, "08_perfiles_variables_top.png"), g, width = 13, height = 7.5, dpi = 200)

## tabla de correlaciones entre variables seleccionadas de distintos bloques
sel <- lapply(names(X), function(b) unique(c(selectVar(fit, block = b, comp = 1)[[b]]$name,
                                              selectVar(fit, block = b, comp = 2)[[b]]$name)))
names(sel) <- names(X)
pares <- do.call(rbind, combn(names(X), 2, simplify = FALSE, FUN = function(p)
  do.call(rbind, lapply(sel[[p[1]]], function(a) data.frame(bloque_1 = p[1], variable_1 = a, bloque_2 = p[2],
    variable_2 = sel[[p[2]]], r = as.numeric(cor(X[[p[1]]][, a], X[[p[2]]][, sel[[p[2]]], drop = FALSE])))))))
pares <- pares[order(-abs(pares$r)), ]
write.table(pares, file.path(res, "07_DIABLO_correlaciones_variables_seleccionadas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("figuras en", fig, "\n")
