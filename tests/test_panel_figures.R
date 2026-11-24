#!/usr/bin/env Rscript
# tests/test_panel_figures.R
# ==============================================================================
# Unit test of the whole-panel forest and table on SYNTHETIC values (R/modules_pub_figures.R;
# R/modules_pub_style.R marker labels):
#   S13 LMM forest: Satterthwaite df recovered from (t, P), CI vs P agreement, the report
#       reader with and without the bootstrap sheet, the figure (pub_lmm_t_ci(),
#       pub_read_lmm_report(), pub_fig_lmm_forest_all());
#   whole-panel LMM table: columns, order, rounding, xlsx round trip (pub_table_lmm_all());
#   group order / colours shared with Figure 1B (pub_group_scale()).
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

# --- 5. group scale shared with Figure 1B ------------------------------------------------
gsc <- pub_group_scale(c("SD/PD", "PR", "PR"), "PR")
check(identical(gsc$levels, c("PR", "SD/PD")) && identical(unname(gsc$colours), unname(pub_palette[c("combined", "unpen")])),
      "group scale: responder first, blue")
check(is.null(pub_group_scale(c("A", "B"), "PR")$levels), "group scale: responder absent keeps the data order")

# --- 6. whole-panel LMM table: Marker, Beta, P, FDR ----------------------------------------
tb <- pub_table_lmm_all(r_full$all)
check(identical(names(tb), c("Marker", "Beta", "P", "FDR")), "table: exactly Marker, Beta, P, FDR")
check(nrow(tb) == nrow(d), "table: one row per marker")
o <- order(d$FDR_Interaction, d$P_Value_Interaction)
check(identical(tb$Marker, pub_relabel_marker(d$Marker[o])), "table: sorted by FDR then P, display names")
check(identical(tb$Beta, round(d$Estimate_Interaction[o], 2)) && identical(tb$P, signif(d$P_Value_Interaction[o], 2)) &&
      identical(tb$FDR, signif(d$FDR_Interaction[o], 2)), "table: beta 2 decimals, P / FDR 2 significant digits")
check(is.null(pub_table_lmm_all(NULL)), "table: NULL input -> NULL")
check(errs(pub_table_lmm_all(na_fdr)), "table: missing FDR: error expected")
ft <- file.path(work, "Table_LMM_AllMarkers_test.xlsx")
openxlsx::write.xlsx(tb, ft, overwrite = TRUE)
rb <- as.data.frame(readxl::read_excel(ft))
check(identical(readxl::excel_sheets(ft), "Sheet 1") && isTRUE(all.equal(rb, tb, check.attributes = FALSE)),
      "table: one sheet, read back unchanged")

unlink(work, recursive = TRUE)
if (length(fails)) {
  cat(sprintf("[FAIL] %d check(s):\n", length(fails))); cat(sprintf("  - %s\n", fails), sep = "")
  quit(status = 1)
}
cat("[OK] whole-panel forest (S13) and LMM table: every check passed.\n")
