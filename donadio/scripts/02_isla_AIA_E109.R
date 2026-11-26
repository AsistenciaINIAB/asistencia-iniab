#!/usr/bin/env Rscript
## =============================================================================
## 02 - Región genómica de E109 que responde al AIA de Az39
## Grafica log2FC (vs E109 control) a lo largo del cromosoma alrededor del
## cluster iac (RN69_30050-30290, ~6.31-6.37 Mb) para AIA exógeno, Az39 y Az39-ipdC.
## Uso: Rscript 02_isla_AIA_E109.R <directorio_de_trabajo>   (después de 01_deseq2_E109.R)
## =============================================================================
suppressMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE); wd <- if (length(args)) args[1] else "."
m <- read.delim(file.path(wd, "results", "E109", "04_tabla_maestra_log2FC_padj.tsv"), quote = "")
fig_dir <- file.path(wd, "figures", "E109")

win <- m[!is.na(m$start) & m$start > 6.25e6 & m$end < 6.43e6, ]
long <- do.call(rbind, lapply(c("AIA_vs_Control", "Az39_vs_Control", "ipdC_vs_Control", "Az39_vs_ipdC"), function(n)
  data.frame(locus_tag = win$locus_tag, pos = (win$start + win$end) / 2e6, contraste = n,
             log2FC = win[[paste0("log2FC_", n)]], padj = win[[paste0("padj_", n)]])))
long$contraste <- factor(long$contraste,
  levels = c("AIA_vs_Control", "Az39_vs_Control", "ipdC_vs_Control", "Az39_vs_ipdC"),
  labels = c("E109 + AIA vs control", "E109 + Az39 vs control",
             "E109 + Az39-ipdC vs control", "Az39 vs Az39-ipdC (efecto del AIA de Az39)"))
long$sig <- !is.na(long$padj) & long$padj < 0.05
isla <- range((m$start[m$locus_tag %in% c("RN69_30050", "RN69_30290")]) / 1e6) + c(-0.002, 0.003)
iac  <- range(m$start[m$locus_tag %in% c("RN69_30270", "RN69_30290")] / 1e6) + c(-0.001, 0.002)

g <- ggplot(long, aes(pos, log2FC)) +
  annotate("rect", xmin = isla[1], xmax = isla[2], ymin = -Inf, ymax = Inf, fill = "#FFF3C4") +
  annotate("rect", xmin = iac[1], xmax = iac[2], ymin = -Inf, ymax = Inf, fill = "#F9D58A") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_segment(aes(xend = pos, y = 0, yend = log2FC, colour = sig), linewidth = 0.4) +
  geom_point(aes(colour = sig), size = 1) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey60"), labels = c("ns", "padj<0.05"), name = NULL) +
  facet_wrap(~ contraste, ncol = 1) +
  labs(x = "posición en el cromosoma de E109 (Mb, CP010313.1)", y = "log2FC",
       title = "Isla de genes de E109 inducida por el AIA de Az39",
       subtitle = "amarillo: RN69_30050-30290 (~58 kb) | naranja: cluster iac (RN69_30270-30290)") +
  theme_bw() + theme(strip.text = element_text(size = 8))
ggsave(file.path(fig_dir, "05_isla_AIA_region_genomica.pdf"), g, width = 9, height = 8)
ggsave(file.path(fig_dir, "05_isla_AIA_region_genomica.png"), g, width = 9, height = 8, dpi = 200)
cat("listo\n")
