#!/usr/bin/env Rscript
## =============================================================================
## 05 - Metabolómica RMN 1H de sobrenadantes (4 h), matriz de concentraciones
##      relativas MCR-ALS del servicio PLABEM (Matriz_asignada_sinQC.csv)
## Grupos (6 réplicas): E109, E109+AIA, E109+Az39 (EA), Az39 MMAB (AzM), Az39 LB (ALB)
## Réplicas 1-3 y 4-6 = dos tandas (las 1-3 están apareadas con el RNA).
##
## 1) QC: PCA, efecto de la tanda
## 2) Contrastes por señal: log2(x) ~ grupo + tanda (lm), FDR BH
## 3) Separación del co-cultivo (diseño de volúmenes de los Métodos):
##      E   = 0.75 cultivo E109 + 0.25 MMAB fresco
##      AzM = 0.25 cultivo Az39 + 0.75 PSY fresco
##      EA  = 0.75 cultivo E109 + 0.25 cultivo Az39
##    => sin interacción: EA = E + AzM - B   (B = medio fresco, NO medido)
##    Contraste lineal D = EA - E - AzM (escala lineal). Como B >= 0:
##      D > 0 significativo  -> PRODUCCIÓN EMERGENTE por la interacción (robusto)
##      D < 0                -> consumo neto o aporte del medio: requiere blanco
## Extractos (fracción intracelular, EXT): el pellet de EA es una mezcla de células
##    de ambas bacterias en proporción desconocida => sin interacción, EA queda
##    ENTRE E y AzM. Se testea si EA está fuera de ese rango (por encima de ambos
##    o por debajo de ambos, contrastes unilaterales, se toma el p mayor).
## Uso: Rscript 05_metabolomica_SN.R <directorio_de_trabajo> [SN|EXT]
## =============================================================================
suppressMessages({ library(ggplot2); library(pheatmap); library(ggrepel) })
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args)) args[1] else "."
fr <- if (length(args) > 1) args[2] else "SN"
fr_nombre <- c(SN = "Sobrenadantes", EXT = "Extractos intracelulares")[fr]
infile <- c(SN = "metab_SN_matriz_asignada.csv", EXT = "metab_EXT_matriz_MCR.csv")[fr]
sub <- c(SN = "metabolomica", EXT = "metabolomica_extractos")[fr]
res_dir <- file.path(wd, "results", sub); fig_dir <- file.path(wd, "figures", sub)
dir.create(res_dir, TRUE, TRUE); dir.create(fig_dir, TRUE, TRUE)

raw <- read.csv2(file.path(wd, "data", infile), header = FALSE,
                 check.names = FALSE, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
samples <- as.character(raw[1, -1]); grp_raw <- as.character(raw[2, -1])
X <- as.matrix(apply(raw[-(1:2), -1], 2, function(v) as.numeric(gsub(",", ".", v))))
rownames(X) <- trimws(raw[-(1:2), 1]); colnames(X) <- samples
grupo <- factor(c("Az39 LB" = "AzLB", "Az39 MMAB" = "AzM", "E109" = "E", "E109 + Az39" = "EA",
                  "E109 + AIA" = "Eaia", ALB = "AzLB", Az39 = "AzM", EA = "EA", EAIA = "Eaia")[trimws(grp_raw)],
                levels = c("E", "Eaia", "EA", "AzM", "AzLB"))
rep  <- as.integer(sub("^[A-Za-z]+", "", samples))
tanda <- factor(ifelse(rep <= 3, "T1", "T2"))
meta <- data.frame(sample = samples, grupo, rep, tanda)
write.table(meta, file.path(res_dir, "00_muestras.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## señal asignada vs no identificada; nombre de metabolito limpio
asignado <- !grepl("^U_", rownames(X))
metab <- ifelse(asignado, trimws(sub("\\s*\\(\\?\\)|/ruido|\\s*\\d+$|\\s*\\+ algo", "", rownames(X))), NA)
metab[metab %in% c("Fenilalanina4")] <- "Fenilalanina"   # 'Fenilalanina44' -> idem
metab <- sub("\\d+$", "", metab)
cat(sprintf("señales: %d (asignadas %d a %d metabolitos; no identificadas %d)\n",
            nrow(X), sum(asignado), length(unique(na.omit(metab))), sum(!asignado)))

## ------------------------------------------------------------------ QC ------
eps <- min(X[X > 0]) / 2
L <- log2(X + eps)
keep <- apply(X, 1, function(v) mean(v > 1e-3) >= 0.25)        # señal presente en >=25% de muestras
L <- L[keep, ]; Xk <- X[keep, ]
pc <- prcomp(t(L), center = TRUE, scale. = TRUE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
pd <- data.frame(pc$x[, 1:3], meta)
cols <- c(E = "#4D4D4D", Eaia = "#E69F00", EA = "#0072B2", AzM = "#009E73", AzLB = "#CC79A7")
g <- ggplot(pd, aes(PC1, PC2, colour = grupo, shape = tanda, label = sample)) + geom_point(size = 3) +
  geom_text_repel(size = 2.5, show.legend = FALSE) + scale_colour_manual(values = cols) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = paste(fr_nombre, "- PCA (log2, autoescalado)"), shape = "Tanda (réplicas)") + theme_bw()
ggsave(file.path(fig_dir, paste0("01_PCA_", fr, ".png")), g, width = 7, height = 5.2, dpi = 200)
ggsave(file.path(fig_dir, paste0("01_PCA_", fr, ".pdf")), g, width = 7, height = 5.2)

## efecto de tanda: por señal, F-test de 'tanda' dentro de grupo
fit_tanda <- t(apply(L, 1, function(y) {
  a <- anova(lm(y ~ grupo + tanda)); c(F_tanda = a["tanda", "F value"], p_tanda = a["tanda", "Pr(>F)"]) }))
qt <- data.frame(senal = rownames(L), fit_tanda, FDR_tanda = p.adjust(fit_tanda[, "p_tanda"], "BH"))
write.table(qt[order(qt$p_tanda), ], file.path(res_dir, "01_efecto_tanda_por_senal.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat(sprintf("señales con efecto de tanda (FDR<0.05): %d de %d\n", sum(qt$FDR_tanda < 0.05), nrow(qt)))

## ----------------------------------------------------- contrastes (log2) ----
contr <- list(Eaia_vs_E = c("Eaia", "E"), EA_vs_E = c("EA", "E"), EA_vs_Eaia = c("EA", "Eaia"),
              EA_vs_AzM = c("EA", "AzM"), AzM_vs_AzLB = c("AzM", "AzLB"))
one <- function(num, den) {
  r <- t(apply(L, 1, function(y) {
    g2 <- relevel(grupo, ref = den)
    co <- summary(lm(y ~ g2 + tanda))$coefficients
    k <- paste0("g2", num); c(log2FC = co[k, 1], t = co[k, 3], p = co[k, 4]) }))
  data.frame(senal = rownames(L), metabolito = metab[keep], r, FDR = p.adjust(r[, "p"], "BH"))
}
R <- lapply(contr, function(x) one(x[1], x[2]))
resumen <- do.call(rbind, lapply(names(R), function(n) { d <- R[[n]]; s <- d$FDR < 0.05
  data.frame(contraste = n, senales_FDR005 = sum(s), suben = sum(s & d$log2FC > 0), bajan = sum(s & d$log2FC < 0),
             asignadas_FDR005 = sum(s & !is.na(d$metabolito))) }))
print(resumen)
write.table(resumen, file.path(res_dir, "02_resumen_contrastes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (n in names(R)) write.table(R[[n]][order(R[[n]]$p), ], file.path(res_dir, paste0("DA_", n, ".tsv")),
                                sep = "\t", quote = FALSE, row.names = FALSE)

## ------------------------------------------- separación del co-cultivo -----
## contraste lineal en escala original: D = EA - E - AzM, con tanda como covariable
dec <- t(apply(Xk, 1, function(y) {
  m <- lm(y ~ 0 + grupo + tanda); b <- coef(m); V <- vcov(m)
  w <- setNames(rep(0, length(b)), names(b)); w["grupoEA"] <- 1; w["grupoE"] <- -1; w["grupoAzM"] <- -1
  D <- sum(w * b); se <- sqrt(as.numeric(t(w) %*% V %*% w)); df <- m$df.residual
  c(EA = unname(b["grupoEA"]), E = unname(b["grupoE"]), AzM = unname(b["grupoAzM"]),
    esperado_sin_blanco = unname(b["grupoE"] + b["grupoAzM"]), D = D, t = D / se,
    p_produccion = pt(D / se, df, lower.tail = FALSE), p_dos_colas = 2 * pt(-abs(D / se), df)) }))
dec <- data.frame(senal = rownames(Xk), metabolito = metab[keep], dec)
if (fr == "EXT") {
  ## fuera del rango [E, AzM]: unilateral contra el máximo (arriba) y contra el mínimo (abajo)
  rng <- t(apply(Xk, 1, function(y) {
    m <- lm(y ~ 0 + grupo + tanda); b <- coef(m); V <- vcov(m); df <- m$df.residual
    ct <- function(w) { w2 <- setNames(rep(0, length(b)), names(b)); w2[names(w)] <- w
      d <- sum(w2 * b); se <- sqrt(as.numeric(t(w2) %*% V %*% w2)); c(d, pt(d / se, df, lower.tail = FALSE)) }
    a1 <- ct(c(grupoEA = 1, grupoE = -1)); a2 <- ct(c(grupoEA = 1, grupoAzM = -1))
    b1 <- ct(c(grupoE = 1, grupoEA = -1)); b2 <- ct(c(grupoAzM = 1, grupoEA = -1))
    c(p_arriba = max(a1[2], a2[2]), p_abajo = max(b1[2], b2[2])) }))
  dec$p_arriba <- rng[, "p_arriba"]; dec$p_abajo <- rng[, "p_abajo"]
  dec$FDR_arriba <- p.adjust(dec$p_arriba, "BH"); dec$FDR_abajo <- p.adjust(dec$p_abajo, "BH")
}
dec$FDR_produccion <- p.adjust(dec$p_produccion, "BH")
dec$FDR_dos_colas  <- p.adjust(dec$p_dos_colas, "BH")
dec$clase <- ifelse(dec$FDR_produccion < 0.05, "Producción emergente (EA > E + AzM)",
             ifelse(dec$FDR_dos_colas < 0.05 & dec$D < 0, "EA < E + AzM (consumo o aporte del medio; requiere blanco)",
                    "Compatible con suma de cultivos"))
if (fr == "EXT") dec$clase <- ifelse(dec$FDR_arriba < 0.05, "EA por encima de E109 y de Az39 (emergente)",
                          ifelse(dec$FDR_abajo < 0.05, "EA por debajo de E109 y de Az39 (emergente)",
                                 "Dentro del rango de las dos bacterias"))
dec <- dec[order(dec$clase, dec$p_dos_colas), ]
write.table(dec, file.path(res_dir, "03_separacion_cocultivo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table(dec$clase))
print(head(dec[grepl("emergente|Producción", dec$clase), c("senal", "clase", "EA", "E", "AzM")], 40))

## ------------------------------------------------ resumen por metabolito -----
## metabolitos asignados: media de sus señales (escala lineal), por grupo
ma <- Xk[!is.na(metab[keep]), ]; mm <- metab[keep][!is.na(metab[keep])]
agg <- rowsum(ma, mm) / as.vector(table(mm)[rownames(rowsum(ma, mm))])
gm <- sapply(levels(grupo), function(gn) rowMeans(agg[, grupo == gn, drop = FALSE]))
write.table(data.frame(metabolito = rownames(gm), round(gm, 4)), file.path(res_dir, "04_medias_por_metabolito.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
z <- t(scale(t(log2(agg + eps))))
ann_c <- data.frame(Grupo = grupo, Tanda = tanda, row.names = samples)
pheatmap(z, annotation_col = ann_c, annotation_colors = list(Grupo = cols, Tanda = c(T1 = "grey30", T2 = "grey80")),
         cluster_cols = TRUE, main = paste(fr_nombre, "- metabolitos asignados (z-score log2)"),
         filename = file.path(fig_dir, "02_heatmap_metabolitos_asignados.pdf"), width = 10, height = 6.5)
pheatmap(z, annotation_col = ann_c, annotation_colors = list(Grupo = cols, Tanda = c(T1 = "grey30", T2 = "grey80")),
         cluster_cols = TRUE, main = paste(fr_nombre, "- metabolitos asignados (z-score log2)"),
         filename = file.path(fig_dir, "02_heatmap_metabolitos_asignados.png"), width = 10, height = 6.5)

## figura de separación para metabolitos asignados
dm <- dec[!is.na(dec$metabolito), ]
long <- do.call(rbind, lapply(seq_len(nrow(dm)), function(i) data.frame(
  senal = dm$senal[i], clase = dm$clase[i],
  tipo = c("E109 (E)", "Az39 MMAB (AzM)", "E + AzM (suma)", "Co-cultivo (EA)"),
  valor = c(dm$E[i], dm$AzM[i], dm$esperado_sin_blanco[i], dm$EA[i]))))
long$tipo <- factor(long$tipo, levels = c("E109 (E)", "Az39 MMAB (AzM)", "E + AzM (suma)", "Co-cultivo (EA)"))
long$senal <- factor(long$senal, levels = dm$senal[order(dm$D)])
g <- ggplot(long, aes(valor, senal, colour = tipo, shape = tipo)) + geom_point(size = 2) +
  scale_colour_manual(values = c("#4D4D4D", "#009E73", "grey60", "#0072B2"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 4, 15), name = NULL) +
  facet_grid(clase ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(22)) +
  labs(x = "concentración relativa (MCR-ALS, media de grupo ajustada por tanda)", y = NULL,
       title = paste(fr_nombre, "- co-cultivo vs cultivos individuales (señales asignadas)")) +
  theme_bw() + theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5),
                     strip.text.y = element_text(size = 7))
ggsave(file.path(fig_dir, "03_separacion_cocultivo.png"), g, width = 9, height = 10, dpi = 200)
ggsave(file.path(fig_dir, "03_separacion_cocultivo.pdf"), g, width = 9, height = 10)
cat("listo\n")
