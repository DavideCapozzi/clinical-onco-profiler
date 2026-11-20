#!/usr/bin/env Rscript
# tests/test_panel_figures.R
# ==============================================================================
# Unit test of the whole-panel figures on SYNTHETIC values (R/modules_pub_figures.R;
# R/modules_pub_style.R marker labels):
#   S13 LMM forest: Satterthwaite df recovered from (t, P), CI vs P agreement, the report
#       reader with and without the bootstrap sheet, the figure (pub_lmm_t_ci(),
#       pub_read_lmm_report(), pub_fig_lmm_forest_all());
#   S14 cell frequencies: the Step-01 reader (patients, back-transform, the Step-01 rds on
#       disk), box statistics equal to stat_boxplot(), the 3 × IQR display window, marker
#       families / pages, group order / colours shared with Figure 1B, the figure (hidden values
#       and lines, boxes from all values, page size) (pub_panel_values(), pub_read_panel_values(),
#       pub_log_box_stats(), pub_log_window(), pub_marker_family(), pub_group_scale(),
#       pub_fig_panel_frequencies()).
#
# Run:  Rscript tests/test_panel_figures.R      (from the repository root, ~5 s)
# Exit: 0 on pass, 1 on failure.
# ==============================================================================
suppressPackageStartupMessages({
  library(here); library(jsonlite)
  source(here("R/utils_io.R")); source(here("R/utils_figures.R"))
  source(here("R/modules_pub_style.R")); source(here("R/modules_pub_figures.R"))
})

fails <- character(0)
check <- function(ok, what) if (!isTRUE(ok)) fails <<- c(fails, what)
errs  <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")

# Synthetic interaction table: known df (incl. a singular model's n - 4 = 130), P from the t test.
mk <- function(marker, beta, se, df) {
  t <- beta / se
  data.frame(Marker = marker, Estimate_Interaction = beta, Std_Error = se, T_Value_Interaction = t,
             P_Value_Interaction = 2 * pt(-abs(t), df), N_Observations = 134L)
}
df_true <- c(46.37, 130, 70.2, 52.8, 58)
d <- mk(c("KI67NAIVE", "CD28KI67", "KI67CD8", "PD1EFF", "CM"),
        c(-1.022, -0.866, -0.591, 0.445, 0.027), c(0.315, 0.264, 0.2877, 0.3712, 0.1495), df_true)
d$FDR_Interaction <- p.adjust(d$P_Value_Interaction, "BH")

# --- 1. pub_lmm_t_ci(): df recovery, CI, agreement with P --------------------------
ci <- pub_lmm_t_ci(d)
check(all(abs(ci$DF - df_true) / df_true < 1e-6), sprintf("df recovered: %s", paste(round(ci$DF, 4), collapse = ",")))
q <- qt(0.975, df_true)
check(all(abs(ci$CI_Lower - (d$Estimate_Interaction - q * d$Std_Error)) < 1e-8) &&
      all(abs(ci$CI_Upper - (d$Estimate_Interaction + q * d$Std_Error)) < 1e-8), "CI = beta -/+ qt * SE")
check(identical(ci$CI_Lower > 0 | ci$CI_Upper < 0, d$P_Value_Interaction < 0.05), "CI excludes 0 iff P < 0.05")
d_df <- d; d_df$DF_Interaction <- df_true
check(identical(pub_lmm_t_ci(d_df)$DF, df_true), "DF_Interaction column used as is")
d_df$DF_Interaction <- 2                                        # df contradicting the stored P
check(errs(pub_lmm_t_ci(d_df)), "DF_Interaction inconsistent with P: error expected")
z <- mk("NAIVE", 0, 0.2, 50)                                    # t = 0: df unidentified
check(is.infinite(pub_lmm_t_ci(z)$DF) && abs(pub_lmm_t_ci(z)$CI_Upper - qnorm(0.975) * 0.2) < 1e-12,
      "t = 0 takes the normal quantile")
bad <- d; bad$P_Value_Interaction[1] <- 0.5                     # P inconsistent with t
check(errs(pub_lmm_t_ci(bad)), "P inconsistent with t: error expected")
check(errs(pub_lmm_t_ci(d[, c("Marker", "Estimate_Interaction")])), "missing columns: error expected")

# --- 2. pub_read_lmm_report(): sheets present / absent ---------------------------------
work <- tempfile("lmm_forest_"); dir.create(work)
wr <- function(f, sheets) {
  wb <- openxlsx::createWorkbook()
  for (s in names(sheets)) { openxlsx::addWorksheet(wb, s); openxlsx::writeData(wb, s, sheets[[s]]) }
  openxlsx::saveWorkbook(wb, f, overwrite = TRUE); f
}
boot <- data.frame(Marker = c("KI67NAIVE", "CD28KI67"), Median_Beta_Boot = c(-1, -0.8), CI_Lower_2.5 = c(-1.6, -1.3),
                   CI_Upper_97.5 = c(-0.4, -0.4), Bootstrap_P = 0, Pct_FDR_Significant = c(55, 59), N_Valid_Iterations = 500L)
r_full <- pub_read_lmm_report(wr(file.path(work, "full.xlsx"), list(LMM_Interaction_Results = d, LMM_Bootstrap_CI = boot)))
check(!is.null(r_full$boot) && nrow(r_full$all) == nrow(d) && identical(names(r_full$obs),
      c("Marker", "Estimate_Interaction", "FDR_Interaction")), "reader: both sheets")
check(all(c("Std_Error", "T_Value_Interaction", "P_Value_Interaction") %in% names(r_full$all)), "reader: all has the test columns")
r_noboot <- pub_read_lmm_report(wr(file.path(work, "noboot.xlsx"), list(LMM_Interaction_Results = d)))
check(is.null(r_noboot$boot) && nrow(r_noboot$all) == nrow(d), "reader: no bootstrap sheet keeps the panel")
check(is.null(pub_read_lmm_report(wr(file.path(work, "none.xlsx"), list(Other = data.frame(a = 1))))), "reader: neither sheet -> NULL")
check(is.null(pub_read_lmm_report(file.path(work, "absent.xlsx"))), "reader: missing file -> NULL")

# --- 3. pub_fig_lmm_forest_all() -------------------------------------------------------
gate <- c("KI67NAIVE", "CD28KI67")
p <- pub_fig_lmm_forest_all(r_full$all, gate = gate, ref_label = "SD/PD")
check(inherits(p, "ggplot"), "figure is a ggplot")
b <- ggplot2::ggplot_build(p)
lev <- levels(p$data$Marker)
check(identical(lev, pub_relabel_marker(d$Marker[order(d$Estimate_Interaction)])), "rows sorted by beta, display names")
pts <- b$data[[3]]                                               # geom_point layer
sel <- pts$colour[match(c("Ki67⁺NAÏVE", "CD28⁺Ki67⁺"), lev[pts$y])]
check(all(sel == pub_palette[["combined"]]) && all(pts$colour[!lev[pts$y] %in% c("Ki67⁺NAÏVE", "CD28⁺Ki67⁺")] == "grey55"),
      "gate markers blue, others grey")
check(identical(p$labels$x, "Time×Group β\n(logit scale; ref. SD/PD)"), "x title as Figure 1A")
check(is.null(pub_fig_lmm_forest_all(NULL)), "NULL input -> NULL")
na_fdr <- r_full$all; na_fdr$FDR_Interaction[2] <- NA
check(errs(pub_fig_lmm_forest_all(na_fdr)), "missing FDR: error expected")
f <- save_pub_figure(p, file.path(work, "FigureS13_test"), PUB_W1, 34 + 4.2 * nrow(d))
check(!is.null(f) && file.exists(f) && file.size(f) > 0, "figure saved")

# --- 4. marker labels ------------------------------------------------------------------
check(!anyDuplicated(PUB_MARKER_LABELS), "display names are unique")
check(identical(pub_relabel_marker(c("KI657EMRA", "CD137tot", "UNKNOWN")), c("Ki67⁺EMRA", "CD137⁺", "UNKNOWN")),
      "relabel: mapped codes + passthrough")
check(identical(pub_fmt_p(c(0.0419, 0.9, 0.0004)), c("0.042", "0.90", "< 0.001")), "pub_fmt_p")

# --- 5. S14 panel values: Step-01 reader -----------------------------------------------
set.seed(1)
ids  <- sprintf("P%02d", 1:8)
meta <- data.frame(Patient_ID = c(ids, ids[1:6]), Timepoint = c(rep("T0", 8), rep("T1", 6)),
                   Group = c(rep(c("RP", "SD_PD"), 4), rep(c("RP", "SD_PD"), 3)), stringsAsFactors = FALSE)
pct  <- matrix(runif(nrow(meta) * 3, 0.5, 60), ncol = 3, dimnames = list(NULL, c("KI67NAIVE", "CD28KI67", "CM")))
raw  <- cbind(meta, as.data.frame(qlogis(pct / 100)))          # Step-01 layout: metadata columns + logit markers
DL   <- list(metadata = meta, hybrid_markers = colnames(pct), hybrid_data_raw = raw)
pv   <- pub_panel_values(DL, ids[1:6])
check(nrow(pv) == 6 * 2 * 3 && setequal(pv$Patient_ID, ids[1:6]), "panel values: only the requested patients, every marker")
check(isTRUE(all.equal(pv$Value_pct[pv$Marker == "CM"], pct[meta$Patient_ID %in% ids[1:6], "CM"])),
      "panel values: logit back to %")
check(identical(levels(pv$Timepoint), c("T0", "T1")), "panel values: timepoint levels")
check(is.null(pub_panel_values(DL, "nobody")), "panel values: no patient -> NULL")
cfg1 <- list(output_root = file.path(work, "run", "Synthetic"), project_name = "Synthetic")
rds1 <- step_output_path(cfg1, 1, "data_processed_Synthetic_longitudinal", "rds")
dir.create(dirname(rds1), recursive = TRUE, showWarnings = FALSE); saveRDS(DL, rds1)
check(identical(pub_read_panel_values(cfg1, ids[1:6]), pv), "panel values: read from the Step-01 rds")
check(is.null(pub_read_panel_values(list(output_root = file.path(work, "none"), project_name = "X"), ids)),
      "panel values: no rds -> NULL")

# --- 6. S14 helpers: box statistics, display window, marker families, group scale ------
v <- c(10^rnorm(40, 0, 0.4), 1e-4, 300)                         # bulk + a floor zero + a high value
bs <- pub_log_box_stats(v)
gb <- ggplot2::layer_data(ggplot2::ggplot(data.frame(x = 1, y = v), ggplot2::aes(x, y)) +
                          ggplot2::geom_boxplot() + ggplot2::scale_y_log10())
check(isTRUE(all.equal(unname(log10(bs)), unname(unlist(gb[1, c("ymin", "lower", "middle", "upper", "ymax")])))),
      "box stats = stat_boxplot under scale_y_log10")
q  <- quantile(log10(v), c(0.25, 0.75), names = FALSE)
w3 <- pub_log_window(v)
check(isTRUE(all.equal(log10(w3), c(q[1] - 3 * diff(q), q[2] + 3 * diff(q)))), "window = quartiles -/+ 3 IQR (log)")
check(w3[1] > 1e-4 && w3[2] < 300, "window excludes the far-out values")
check(isTRUE(all.equal(pub_log_window(v, box_ends = c(1e-3, 1e3)), c(1e-3, 1e3))), "window widened to the box ends")
fam <- pub_marker_family(c("KI67NAIVE", "KI657EMRA", "CD28KI67", "CD137KI67", "PD1KI67", "CD3tot", "CM", "EMRA"))
check(identical(as.character(fam), c("Ki67", "Ki67", "CD28", "CD137", "PD1", rep("Lineage and memory", 3))) &&
      identical(levels(fam), c("Ki67", "CD28", "CD137", "PD1", "Lineage and memory")), "marker families and page order")
gsc <- pub_group_scale(c("SD/PD", "PR", "PR"), "PR")
check(identical(gsc$levels, c("PR", "SD/PD")) && identical(unname(gsc$colours), unname(pub_palette[c("combined", "unpen")])),
      "group scale: responder first, blue")
check(is.null(pub_group_scale(c("A", "B"), "PR")$levels), "group scale: responder absent keeps the data order")

# --- 7. S14 figure: pages, window, boxes from all values ---------------------------------
# 40 paired patients (20 per group): a realistic cell size, so the far-out value is beyond
# the whiskers of its own box as well as beyond the marker's window.
ids2  <- sprintf("Q%02d", 1:40)
meta2 <- data.frame(Patient_ID = rep(ids2, 2), Timepoint = rep(c("T0", "T1"), each = 40),
                    Group = rep(rep(c("RP", "SD_PD"), 20), 2), stringsAsFactors = FALSE)
pct2  <- sapply(c(KI67NAIVE = 1, CD28KI67 = 5, CM = 25), function(b) b * 10^rnorm(80, 0, 0.2))
pv2   <- pub_panel_values(list(metadata = meta2, hybrid_markers = colnames(pct2),
                               hybrid_data_raw = cbind(meta2, as.data.frame(qlogis(pct2 / 100)))), ids2)
i_out <- which(pv2$Marker == "CM")[1]; pv2$Value_pct[i_out] <- 1e-4                # one far-out value
lmm3 <- r_full$all[r_full$all$Marker %in% colnames(pct), ]
pg <- suppressMessages(pub_fig_panel_frequencies(pv2, lmm3, resp_label = "RP", ncol = 1, nrow = 1))
check(is.list(pg) && length(pg) == 3 && all(vapply(pg, inherits, TRUE, "ggplot")), "S14: one page per family")
check(identical(vapply(pg, function(p) p$labels$tag, ""), c("A", "B", "C")), "S14: pages tagged A, B, C")
lab_of <- function(m) sprintf("%s\n%s", pub_relabel_marker(m), pub_lab_p(d$FDR_Interaction[d$Marker == m], "FDR"))
check(identical(vapply(pg, function(p) levels(p$data$Panel), ""), vapply(c("KI67NAIVE", "CD28KI67", "CM"), lab_of, "", USE.NAMES = FALSE)),
      "S14: Ki67, CD28, lineage/memory pages with beta-ordered FDR strips")
check(identical(levels(pg[[1]]$data$x), c("PR\nT0", "PR\nT1", "SD/PD\nT0", "SD/PD\nT1")), "S14 x: responder pair first")
check(identical(unname(attr(pg, "hidden")), c(1L, nrow(pv2))), sprintf("S14 hidden count: %s", paste(attr(pg, "hidden"), collapse = "/")))
cm <- pg[[3]]
check(!any(cm$data$Value_pct <= 1e-4), "S14: the far-out value is not drawn")
check(!(pv2$Patient_ID[i_out] %in% cm$layers[[1]]$data$Patient_ID[cm$layers[[1]]$data$Timepoint == pv2$Timepoint[i_out]]),
      "S14: its patient line is not drawn")
cell <- pv2[pv2$Marker == "CM" & pv2$Group == pv2$Group[i_out] & pv2$Timepoint == pv2$Timepoint[i_out], "Value_pct"]
bx <- cm$layers[[2]]$data
bx <- bx[as.character(bx$x) == sprintf("%s\n%s", pub_relabel_group(pv2$Group[i_out]), pv2$Timepoint[i_out]), ]
check(isTRUE(all.equal(unname(unlist(bx[, c("ymin", "lower", "middle", "upper", "ymax")])), unname(pub_log_box_stats(cell)))),
      "S14: the box of the hidden value's cell uses all its values")
pn <- suppressMessages(pub_fig_panel_frequencies(pv, ncol = 1, nrow = 1))   # no LMM table: names, no FDR
check(identical(vapply(pn, function(p) levels(p$data$Panel), ""), pub_relabel_marker(c("KI67NAIVE", "CD28KI67", "CM"))),
      "S14 without LMM table: names only")
check(is.null(pub_fig_panel_frequencies(NULL)), "S14: NULL input -> NULL")
pd3 <- suppressMessages(pub_fig_panel_frequencies(pv, lmm3, resp_label = "RP"))  # default 3 x 3: short pages padded
check(all(vapply(pd3, inherits, TRUE, "patchwork")), "S14: short pages keep the 3-row panel height (spacer)")
f14 <- save_pub_figure(pd3, file.path(work, "FigureS14_test"), PUB_W2, 190)
np  <- suppressWarnings(system2("pdfinfo", f14, stdout = TRUE))
check(!is.null(f14) && any(grepl("^Pages:\\s+3$", np)), sprintf("S14 saved as a 3-page pdf: %s", paste(grep("Pages", np, value = TRUE), collapse = "")))

unlink(work, recursive = TRUE)
if (length(fails)) {
  cat(sprintf("[FAIL] %d check(s):\n", length(fails))); cat(sprintf("  - %s\n", fails), sep = "")
  quit(status = 1)
}
cat("[OK] whole-panel figures (S13, S14): every check passed.\n")
