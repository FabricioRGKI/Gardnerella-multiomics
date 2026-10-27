
# ============================================================
# UpSet — combined CSV, split into UP and DOWN plots
#   * Labels = top-3 MODULES per intersection
#   * DOWN plot: intersection bars grow DOWNWARD
# ============================================================

run_upset_analysis_4way <- function(csv_files,
                                    colors_up   = c(Lc = "#E3377E", Li = "#B02CBA",
                                                    Fv = "#4497AF", Pb = "#56BC37"),
                                    colors_down = c(Lc = "#f19ec1", Li = "#e5a6ea",
                                                    Fv = "#b0d5e0", Pb = "#b8e6aa"),
                                    top_n       = 3,
                                    width       = 10,
                                    height      = 6,
                                    base_size   = 10,
                                    label_size  = 3,
                                    count_size  = 3.2,
                                    title_size  = 14,
                                    out_prefix  = "UpSet") {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2)
    library(openxlsx); library(cowplot)
  })
  
  if (length(csv_files) != 1) stop("Provide exactly one combined CSV file.")
  csv_file <- csv_files[[1]]
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  required_cols <- c("locus_tag", "Log2FC", "condition", "regulation", "modules")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0)
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  
  cond_map <- c(Gv_LcGv = "Lc", Gv_LiGv = "Li", Gv_FvGv = "Fv", Gv_PbGv = "Pb")
  if (!all(df$condition %in% names(cond_map)))
    stop("Unexpected condition values.")
  df$cond_short <- cond_map[df$condition]
  df$regulation <- tolower(df$regulation)
  
  cond_levels <- c("Lc", "Li", "Fv", "Pb")
  
  # ---- Helper: build one UpSet panel ----
  # flipped = TRUE  -> intersection bars grow DOWN, placed under the matrix
  # flipped = FALSE -> classic layout (bars up, above the matrix)
  make_upset <- function(regulation_dir, palette, plot_title, flipped = FALSE) {
    sub_df <- df[df$regulation == regulation_dir, , drop = FALSE]
    if (nrow(sub_df) == 0) return(NULL)
    
    set_names <- cond_levels
    df_list <- list()
    for (s in set_names) {
      d <- sub_df[sub_df$cond_short == s, , drop = FALSE]
      d <- d[!duplicated(d$locus_tag), , drop = FALSE]
      df_list[[s]] <- d
    }
    names(df_list) <- set_names
    
    all_tags <- unique(unlist(lapply(df_list, `[[`, "locus_tag")))
    if (length(all_tags) == 0) return(NULL)
    
    memb <- data.frame(locus_tag = all_tags, stringsAsFactors = FALSE)
    for (nm in set_names) memb[[nm]] <- all_tags %in% df_list[[nm]]$locus_tag
    memb$pattern <- apply(memb[, set_names, drop = FALSE], 1,
                          function(r) paste(set_names[r], collapse = "&"))
    
    inter_sizes <- memb %>%
      count(pattern, name = "size") %>%
      arrange(desc(size)) %>%
      mutate(pattern_id = factor(pattern, levels = pattern))
    
    # ---- Top-N MODULE labels ----
    top_modules <- function(tags, df_list, n = 3, show_all = FALSE) {
      rows <- do.call(rbind, lapply(names(df_list), function(nm) {
        d <- df_list[[nm]]
        if (nrow(d) == 0 || !"modules" %in% colnames(d)) return(NULL)
        data.frame(locus_tag = d$locus_tag,
                   module    = d$modules,
                   Log2FC    = d$Log2FC,
                   stringsAsFactors = FALSE) %>%
          filter(locus_tag %in% tags, !is.na(module), module != "")
      }))
      if (is.null(rows) || nrow(rows) == 0) return("")
      # rank modules by #features in intersection, tie-break by max |Log2FC|
      ranked <- rows %>%
        group_by(module) %>%
        summarise(n_feat   = n(),
                  max_abs  = max(abs(Log2FC), na.rm = TRUE),
                  .groups  = "drop") %>%
        arrange(desc(n_feat), desc(max_abs))
      mods <- if (show_all) ranked$module else head(ranked$module, n)
      paste(mods, collapse = "\n")
    }
    inter_sizes$labels <- vapply(inter_sizes$pattern, function(p) {
      tags <- memb$locus_tag[memb$pattern == p]
      show_all <- length(strsplit(p, "&")[[1]]) == length(set_names)
      top_modules(tags, df_list, n = top_n, show_all = show_all)
    }, character(1))
    
    # ---- Color per intersection (average of set colors) ----
    inter_sizes$bar_color <- vapply(inter_sizes$pattern, function(p) {
      ids <- strsplit(p, "&")[[1]]
      rgb_mat <- col2rgb(palette[ids])
      avg <- rowMeans(rgb_mat)
      rgb(avg[1], avg[2], avg[3], maxColorValue = 255)
    }, character(1))
    
    # ---- Intersection-size panel ----
    if (flipped) {
      # Bars grow DOWN: put y=0 at top; add label above each bar
      p_top <- ggplot(inter_sizes, aes(pattern_id, size, fill = pattern_id)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = labels), vjust = 1.15,   # labels inside/below bar tops
                  size = label_size, lineheight = 0.9) +
        scale_fill_manual(values = setNames(inter_sizes$bar_color,
                                            as.character(inter_sizes$pattern_id))) +
        scale_y_reverse(expand = expansion(mult = c(0.45, 0))) +  # grow downward
        labs(y = "Intersection size", x = NULL) +
        theme_minimal(base_size = base_size) +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank())
    } else {
      p_top <- ggplot(inter_sizes, aes(pattern_id, size, fill = pattern_id)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = labels), vjust = -0.2,
                  size = label_size, lineheight = 0.9) +
        scale_fill_manual(values = setNames(inter_sizes$bar_color,
                                            as.character(inter_sizes$pattern_id))) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.45))) +
        labs(y = "Intersection size", x = NULL) +
        theme_minimal(base_size = base_size) +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank())
    }
    
    # ---- Matrix ----
    mat_df <- expand.grid(pattern_id = inter_sizes$pattern_id,
                          set = set_names, stringsAsFactors = FALSE)
    mat_df$present <- mapply(function(p, s)
      grepl(paste0("(^|&)", s, "($|&)"), p),
      mat_df$pattern_id, mat_df$set)
    mat_df$set <- factor(mat_df$set, levels = rev(set_names))
    
    segments_df <- mat_df %>% filter(present) %>%
      mutate(set_idx = match(as.character(set), rev(set_names))) %>%
      group_by(pattern_id) %>%
      summarise(y_min = min(set_idx), y_max = max(set_idx), .groups = "drop")
    
    p_matrix <- ggplot(mat_df, aes(pattern_id, set)) +
      geom_segment(data = segments_df,
                   aes(x = pattern_id, xend = pattern_id,
                       y = rev(set_names)[y_min], yend = rev(set_names)[y_max]),
                   inherit.aes = FALSE, color = "grey40", linewidth = 0.7) +
      geom_point(aes(color = set, alpha = present), size = 4.5) +
      scale_color_manual(values = palette) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.15)) +
      scale_y_discrete(limits = rev(set_names)) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = base_size) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    # ---- Set-size bars ----
    set_sizes <- data.frame(
      set  = factor(set_names, levels = rev(set_names)),
      size = vapply(df_list, function(d) length(unique(d$locus_tag)), integer(1)),
      stringsAsFactors = FALSE
    )
    p_left <- ggplot(set_sizes, aes(size, set, fill = set)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = size), hjust = 1.1, size = count_size) +
      scale_fill_manual(values = palette) +
      scale_x_reverse(expand = expansion(mult = c(0.25, 0))) +
      labs(x = "Set size", y = NULL) +
      theme_minimal(base_size = base_size) +
      theme(legend.position = "none",
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank())
    
    # ---- Assemble ----
    p_blank <- ggplot() + theme_void()
    
    if (flipped) {
      # Order (top -> bottom): set-size+matrix  |  intersection bars (grow down)
      mid_row <- plot_grid(p_left, p_matrix, ncol = 2, rel_widths = c(1, 4),
                           align = "h", axis = "tb")
      bot_row <- plot_grid(p_blank, p_top, ncol = 2, rel_widths = c(1, 4))
      combined <- plot_grid(mid_row, bot_row, ncol = 1,
                            rel_heights = c(1, 2), align = "v", axis = "lr")
    } else {
      # Classic: intersection bars on top, matrix+set-sizes under them
      top_row <- plot_grid(p_blank, p_top, ncol = 2, rel_widths = c(1, 4))
      bot_row <- plot_grid(p_left, p_matrix, ncol = 2, rel_widths = c(1, 4),
                           align = "h", axis = "tb")
      combined <- plot_grid(top_row, bot_row, ncol = 1,
                            rel_heights = c(2, 1), align = "v", axis = "lr")
    }
    
    combined <- ggdraw(combined) +
      draw_label(plot_title, x = 0.5, y = 0.98,
                 hjust = 0.5, vjust = 1, size = title_size)
    
    list(plot = combined, patterns = inter_sizes,
         membership = memb, datasets = df_list, set_sizes = set_sizes)
  }
  
  # ---- Run both directions ----
  res_up   <- make_upset("up",   colors_up,   "UpSet - UPregulated proteins (top modules)",
                         flipped = FALSE)
  res_down <- make_upset("down", colors_down, "UpSet - DOWNregulated proteins (top modules)",
                         flipped = TRUE)
  
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if (!is.null(res_up)) {
    print(res_up$plot)
    f_up <- paste0(out_prefix, "_UP_", ts, ".pdf")
    ggsave(f_up, res_up$plot, width = width, height = height, device = "pdf")
    cat("UP PDF saved:  ", f_up, "\n")
  }
  if (!is.null(res_down)) {
    print(res_down$plot)
    f_dn <- paste0(out_prefix, "_DOWN_", ts, ".pdf")
    ggsave(f_dn, res_down$plot, width = width, height = height, device = "pdf")
    cat("DOWN PDF saved:", f_dn, "\n")
  }
  
  # ---- Excel (unchanged structure) ----
  wb <- createWorkbook()
  as_chr <- function(d) { d[] <- lapply(d, as.character); d }
  
  add_sheets <- function(res, tag) {
    if (is.null(res)) return(invisible())
    addWorksheet(wb, paste0("Patterns_", tag))
    writeData(wb, paste0("Patterns_", tag),
              res$patterns %>% select(pattern, size, labels))
    addWorksheet(wb, paste0("Membership_", tag))
    writeData(wb, paste0("Membership_", tag), res$membership)
    
    sn <- names(res$datasets)
    shared_all_tags <- Reduce(intersect,
                              lapply(res$datasets, function(d) unique(d$locus_tag)))
    shared_all_df <- dplyr::bind_rows(lapply(sn, function(nm) {
      res$datasets[[nm]] %>% filter(locus_tag %in% shared_all_tags) %>%
        mutate(Dataset = nm) %>% as_chr()
    }))
    sh_all <- paste0("Shared_All_", tag)
    addWorksheet(wb, sh_all)
    if (nrow(shared_all_df) > 0) writeData(wb, sh_all, shared_all_df) else
      writeData(wb, sh_all, data.frame(Message = "None shared by all"))
    
    shared_2p <- res$membership$locus_tag[
      rowSums(res$membership[, sn, drop = FALSE]) >= 2]
    shared_2p_df <- dplyr::bind_rows(lapply(sn, function(nm) {
      res$datasets[[nm]] %>% filter(locus_tag %in% shared_2p) %>%
        mutate(Dataset = nm) %>% as_chr()
    }))
    sh_2p <- paste0("Shared_2plus_", tag)
    addWorksheet(wb, sh_2p)
    if (nrow(shared_2p_df) > 0) writeData(wb, sh_2p, shared_2p_df) else
      writeData(wb, sh_2p, data.frame(Message = "None shared by >=2"))
    
    for (nm in sn) {
      others <- setdiff(sn, nm)
      other_tags <- unique(unlist(lapply(res$datasets[others], `[[`, "locus_tag")))
      uniq <- setdiff(unique(res$datasets[[nm]]$locus_tag), other_tags)
      sh <- substr(paste0("Unique_", nm, "_", tag), 1, 31)
      addWorksheet(wb, sh)
      writeData(wb, sh, res$datasets[[nm]] %>% filter(locus_tag %in% uniq))
    }
  }
  add_sheets(res_up, "UP")
  add_sheets(res_down, "DOWN")
  
  xlsx_file <- paste0(out_prefix, "_", ts, ".xlsx")
  saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  cat("Excel saved:", xlsx_file, "\n")
  
  invisible(list(up = res_up, down = res_down))
}

___
#Running example
                                     
result <- run_upset_analysis_4way(
  csv_files   = "Gv_regulated_modules_all.csv",
  colors_up   = c(Lc = "#E3377E", Li = "#B02CBA",
                  Fv = "#4497AF", Pb = "#56BC37"),
  colors_down = c(Lc = "#f19ec1", Li = "#e5a6ea",
                  Fv = "#b0d5e0", Pb = "#b8e6aa"),
  top_n       = 3,
  width       = 5,
  height      = 3
)
