#script for using "locus_tag" column for matching conditions

run_venn_analysis_4way <- function(csv_files, colors = c("#b0d5e0", "#91a0ff", "#b8e6aa", "#ffb6c1")) {
  # Load required libraries
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggvenn)
  library(openxlsx)
  library(svglite)
  
  # Validate inputs
  if (length(csv_files) != 4) {
    stop("Please provide exactly 4 CSV files")
  }
  
  if (length(colors) != 4) {
    stop("Please provide exactly 4 colors")
  }
  
  # Read the four files and store their original names
  file1 <- csv_files[1]
  file2 <- csv_files[2]
  file3 <- csv_files[3]
  file4 <- csv_files[4]
  
  df1 <- read.csv(file1, stringsAsFactors = FALSE)
  df2 <- read.csv(file2, stringsAsFactors = FALSE)
  df3 <- read.csv(file3, stringsAsFactors = FALSE)
  df4 <- read.csv(file4, stringsAsFactors = FALSE)
  
  # Extract base names for dataset labeling (remove .csv extension)
  dataset1_name <- tools::file_path_sans_ext(basename(file1))
  dataset2_name <- tools::file_path_sans_ext(basename(file2))
  dataset3_name <- tools::file_path_sans_ext(basename(file3))
  dataset4_name <- tools::file_path_sans_ext(basename(file4))
  
  # Create shortened names for Excel sheets (max 31 characters)
  short_name1 <- ifelse(nchar(dataset1_name) > 28, 
                        paste0(substr(dataset1_name, 1, 28), "..."), 
                        dataset1_name)
  short_name2 <- ifelse(nchar(dataset2_name) > 28, 
                        paste0(substr(dataset2_name, 1, 28), "..."), 
                        dataset2_name)
  short_name3 <- ifelse(nchar(dataset3_name) > 28, 
                        paste0(substr(dataset3_name, 1, 28), "..."), 
                        dataset3_name)
  short_name4 <- ifelse(nchar(dataset4_name) > 28, 
                        paste0(substr(dataset4_name, 1, 28), "..."), 
                        dataset4_name)
  
  cat("Loaded datasets:\n")
  cat("1.", dataset1_name, "\n")
  cat("2.", dataset2_name, "\n")
  cat("3.", dataset3_name, "\n")
  cat("4.", dataset4_name, "\n\n")
  
  # Verify that locus_tag column exists in all dataframes
  for(i in 1:4) {
    df_name <- paste0("df", i)
    if(!"locus_tag" %in% colnames(get(df_name))) {
      stop(paste("Column 'locus_tag' not found in", csv_files[i]))
    }
  }
  
  # Helper function to get top 3 labels based on highest absolute Log2FC
  get_top_labels <- function(df) {
    if("Log2FC" %in% colnames(df) && "label_edit" %in% colnames(df)) {
      # Calculate absolute Log2FC and sort
      df_temp <- df %>%
        mutate(abs_Log2FC = abs(Log2FC)) %>%
        arrange(desc(abs_Log2FC))
      
      # Get top 3 unique labels (in case of duplicates)
      top_labels <- unique(df_temp$label_edit[!is.na(df_temp$label_edit) & df_temp$label_edit != ""])
      top_labels <- head(top_labels, 3)
      
      if(length(top_labels) > 0) {
        return(paste(top_labels, collapse = "\n"))
      } else {
        return("")
      }
    } else {
      return("")
    }
  }
  
  # Get top labels for each dataset
  top_labels_1 <- get_top_labels(df1)
  top_labels_2 <- get_top_labels(df2)
  top_labels_3 <- get_top_labels(df3)
  top_labels_4 <- get_top_labels(df4)
  
  # Create names for Venn diagram with top labels
  venn_name1 <- if(top_labels_1 != "") paste0(dataset1_name, "\n", top_labels_1) else dataset1_name
  venn_name2 <- if(top_labels_2 != "") paste0(dataset2_name, "\n", top_labels_2) else dataset2_name
  venn_name3 <- if(top_labels_3 != "") paste0(dataset3_name, "\n", top_labels_3) else dataset3_name
  venn_name4 <- if(top_labels_4 != "") paste0(dataset4_name, "\n", top_labels_4) else dataset4_name
  
  # Create a list for Venn diagram with all four datasets using locus_tag
  venn_data <- list(
    unique(df1$locus_tag),
    unique(df2$locus_tag),
    unique(df3$locus_tag),
    unique(df4$locus_tag)
  )
  
  # Set names for the Venn diagram with top labels
  names(venn_data) <- c(venn_name1, venn_name2, venn_name3, venn_name4)
  
  # Create Venn diagram and save as SVG
  venn_plot <- ggvenn(venn_data, 
                      fill_color = colors,
                      stroke_size = 0.4,
                      set_name_size = 4)
  
  # Print the Venn diagram to the screen
  print(venn_plot)
  
  # Save the Venn diagram as SVG
  svg_filename <- paste0("Venn_Diagram_4way_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg")
  svglite(svg_filename, width = 6, height = 4)
  print(venn_plot)
  dev.off()
  
  cat("✅ Venn diagram saved as:", svg_filename, "\n\n")
  
  # Print top labels info
  cat("=== TOP LABELS BY ABSOLUTE Log2FC ===\n")
  cat(dataset1_name, ":", top_labels_1, "\n")
  cat(dataset2_name, ":", top_labels_2, "\n")
  cat(dataset3_name, ":", top_labels_3, "\n")
  cat(dataset4_name, ":", top_labels_4, "\n\n")
  
  # Find all possible overlaps based on locus_tag
  common_all_four <- Reduce(intersect, list(
    df1$locus_tag,
    df2$locus_tag,
    df3$locus_tag,
    df4$locus_tag
  ))
  
  # 4-way intersections
  common_df1_df2_df3 <- intersect(intersect(df1$locus_tag, df2$locus_tag), df3$locus_tag)
  common_df1_df2_df4 <- intersect(intersect(df1$locus_tag, df2$locus_tag), df4$locus_tag)
  common_df1_df3_df4 <- intersect(intersect(df1$locus_tag, df3$locus_tag), df4$locus_tag)
  common_df2_df3_df4 <- intersect(intersect(df2$locus_tag, df3$locus_tag), df4$locus_tag)
  
  # Remove elements that are in all four from the 3-way intersections
  common_df1_df2_df3 <- setdiff(common_df1_df2_df3, common_all_four)
  common_df1_df2_df4 <- setdiff(common_df1_df2_df4, common_all_four)
  common_df1_df3_df4 <- setdiff(common_df1_df3_df4, common_all_four)
  common_df2_df3_df4 <- setdiff(common_df2_df3_df4, common_all_four)
  
  # 2-way intersections (excluding elements in higher-order overlaps)
  common_df1_df2 <- setdiff(intersect(df1$locus_tag, df2$locus_tag), 
                            union(common_all_four, union(common_df1_df2_df3, common_df1_df2_df4)))
  common_df1_df3 <- setdiff(intersect(df1$locus_tag, df3$locus_tag), 
                            union(common_all_four, union(common_df1_df2_df3, common_df1_df3_df4)))
  common_df1_df4 <- setdiff(intersect(df1$locus_tag, df4$locus_tag), 
                            union(common_all_four, union(common_df1_df2_df4, common_df1_df3_df4)))
  common_df2_df3 <- setdiff(intersect(df2$locus_tag, df3$locus_tag), 
                            union(common_all_four, union(common_df1_df2_df3, common_df2_df3_df4)))
  common_df2_df4 <- setdiff(intersect(df2$locus_tag, df4$locus_tag), 
                            union(common_all_four, union(common_df1_df2_df4, common_df2_df3_df4)))
  common_df3_df4 <- setdiff(intersect(df3$locus_tag, df4$locus_tag), 
                            union(common_all_four, union(common_df1_df3_df4, common_df2_df3_df4)))
  
  # Find unique elements for each dataset
  unique_df1 <- setdiff(df1$locus_tag, union(union(df2$locus_tag, df3$locus_tag), df4$locus_tag))
  unique_df2 <- setdiff(df2$locus_tag, union(union(df1$locus_tag, df3$locus_tag), df4$locus_tag))
  unique_df3 <- setdiff(df3$locus_tag, union(union(df1$locus_tag, df2$locus_tag), df4$locus_tag))
  unique_df4 <- setdiff(df4$locus_tag, union(union(df1$locus_tag, df2$locus_tag), df3$locus_tag))
  
  # Print summary information
  cat("=== DATASET SUMMARY ===\n")
  cat("Number of unique locus tags in", dataset1_name, ":", length(unique(df1$locus_tag)), "\n")
  cat("Number of unique locus tags in", dataset2_name, ":", length(unique(df2$locus_tag)), "\n")
  cat("Number of unique locus tags in", dataset3_name, ":", length(unique(df3$locus_tag)), "\n")
  cat("Number of unique locus tags in", dataset4_name, ":", length(unique(df4$locus_tag)), "\n")
  cat("Number of locus tags common to all four datasets:", length(common_all_four), "\n")
  cat("Number of locus tags common to", dataset1_name, dataset2_name, "&", dataset3_name, "only:", length(common_df1_df2_df3), "\n")
  cat("Number of locus tags common to", dataset1_name, dataset2_name, "&", dataset4_name, "only:", length(common_df1_df2_df4), "\n")
  cat("Number of locus tags common to", dataset1_name, dataset3_name, "&", dataset4_name, "only:", length(common_df1_df3_df4), "\n")
  cat("Number of locus tags common to", dataset2_name, dataset3_name, "&", dataset4_name, "only:", length(common_df2_df3_df4), "\n")
  cat("Number of locus tags common to", dataset1_name, "&", dataset2_name, "only:", length(common_df1_df2), "\n")
  cat("Number of locus tags common to", dataset1_name, "&", dataset3_name, "only:", length(common_df1_df3), "\n")
  cat("Number of locus tags common to", dataset1_name, "&", dataset4_name, "only:", length(common_df1_df4), "\n")
  cat("Number of locus tags common to", dataset2_name, "&", dataset3_name, "only:", length(common_df2_df3), "\n")
  cat("Number of locus tags common to", dataset2_name, "&", dataset4_name, "only:", length(common_df2_df4), "\n")
  cat("Number of locus tags common to", dataset3_name, "&", dataset4_name, "only:", length(common_df3_df4), "\n")
  cat("Number of unique locus tags only in", dataset1_name, ":", length(unique_df1), "\n")
  cat("Number of unique locus tags only in", dataset2_name, ":", length(unique_df2), "\n")
  cat("Number of unique locus tags only in", dataset3_name, ":", length(unique_df3), "\n")
  cat("Number of unique locus tags only in", dataset4_name, ":", length(unique_df4), "\n\n")
  
  # Create comprehensive overlap table
  overlap_table <- data.frame()
  
  # Helper function to add rows to overlap table
  add_overlap_row <- function(element, overlap_type, datasets_present, overlap_id_counter) {
    # Determine which datasets have this locus_tag
    has_df1 <- element %in% df1$locus_tag
    has_df2 <- element %in% df2$locus_tag
    has_df3 <- element %in% df3$locus_tag
    has_df4 <- element %in% df4$locus_tag
    
    # Get matches from each dataset
    df1_matches <- if(has_df1) df1 %>% filter(locus_tag == element) else data.frame()
    df2_matches <- if(has_df2) df2 %>% filter(locus_tag == element) else data.frame()
    df3_matches <- if(has_df3) df3 %>% filter(locus_tag == element) else data.frame()
    df4_matches <- if(has_df4) df4 %>% filter(locus_tag == element) else data.frame()
    
    # Helper function to safely extract column values
    safe_extract <- function(df_matches, column_name) {
      if(nrow(df_matches) > 0 && column_name %in% colnames(df_matches)) {
        # Get unique values and collapse them
        values <- unique(df_matches[[column_name]])
        values <- values[!is.na(values) & values != ""]
        if(length(values) > 0) {
          return(paste(values, collapse = "; "))
        } else {
          return("")
        }
      } else {
        return("")
      }
    }
    
    # Extract Log2FC values if column exists
    get_log2fc <- function(df_matches) {
      if(nrow(df_matches) > 0 && "Log2FC" %in% colnames(df_matches)) {
        return(paste(df_matches$Log2FC, collapse = ", "))
      } else {
        return("")
      }
    }
    
    # Get protein and label information from the first available dataset
    protein <- ""
    label <- ""
    
    # Try to get protein and label from df1 first, then df2, etc.
    if(has_df1 && nrow(df1_matches) > 0) {
      protein <- safe_extract(df1_matches, "protein")
      label <- safe_extract(df1_matches, "label")
    }
    
    # If protein/label not found in df1, try other dataframes
    if(protein == "" && has_df2 && nrow(df2_matches) > 0) {
      protein <- safe_extract(df2_matches, "protein")
      label <- safe_extract(df2_matches, "label")
    }
    
    if(protein == "" && has_df3 && nrow(df3_matches) > 0) {
      protein <- safe_extract(df3_matches, "protein")
      label <- safe_extract(df3_matches, "label")
    }
    
    if(protein == "" && has_df4 && nrow(df4_matches) > 0) {
      protein <- safe_extract(df4_matches, "protein")
      label <- safe_extract(df4_matches, "label")
    }
    
    # Create row with Locus_Tag included
    row <- data.frame(
      Overlap_ID = overlap_id_counter,
      Locus_Tag = element,
      Protein = protein,
      Label = label,
      DF1_Log2FC_Values = if(has_df1) get_log2fc(df1_matches) else "",
      DF2_Log2FC_Values = if(has_df2) get_log2fc(df2_matches) else "",
      DF3_Log2FC_Values = if(has_df3) get_log2fc(df3_matches) else "",
      DF4_Log2FC_Values = if(has_df4) get_log2fc(df4_matches) else "",
      Count_DF1 = if(has_df1) nrow(df1_matches) else 0,
      Count_DF2 = if(has_df2) nrow(df2_matches) else 0,
      Count_DF3 = if(has_df3) nrow(df3_matches) else 0,
      Count_DF4 = if(has_df4) nrow(df4_matches) else 0,
      Overlap_Type = overlap_type,
      stringsAsFactors = FALSE
    )
    
    return(row)
  }
  
  # Add elements common to all four
  overlap_id <- 1
  for(element in common_all_four) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, "All_Four", c(1,2,3,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  # Add 3-way overlaps
  for(element in common_df1_df2_df3) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name2, "_", short_name3, "_only"), c(1,2,3), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df1_df2_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name2, "_", short_name4, "_only"), c(1,2,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df1_df3_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name3, "_", short_name4, "_only"), c(1,3,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df2_df3_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name2, "_", short_name3, "_", short_name4, "_only"), c(2,3,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  # Add 2-way overlaps
  for(element in common_df1_df2) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name2, "_only"), c(1,2), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df1_df3) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name3, "_only"), c(1,3), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df1_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name1, "_", short_name4, "_only"), c(1,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df2_df3) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name2, "_", short_name3, "_only"), c(2,3), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df2_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name2, "_", short_name4, "_only"), c(2,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  for(element in common_df3_df4) {
    overlap_table <- rbind(overlap_table, 
                           add_overlap_row(element, paste0(short_name3, "_", short_name4, "_only"), c(3,4), overlap_id))
    overlap_id <- overlap_id + 1
  }
  
  # Update column names to use actual dataset names
  if (nrow(overlap_table) > 0) {
    colnames(overlap_table) <- c("Overlap_ID", "Locus_Tag", "Protein", "Label",
                                 paste0(short_name1, "_Log2FC_Values"),
                                 paste0(short_name2, "_Log2FC_Values"),
                                 paste0(short_name3, "_Log2FC_Values"),
                                 paste0(short_name4, "_Log2FC_Values"),
                                 paste0("Count_", short_name1),
                                 paste0("Count_", short_name2),
                                 paste0("Count_", short_name3),
                                 paste0("Count_", short_name4),
                                 "Overlap_Type")
  }
  
  # Create a workbook
  wb <- createWorkbook()
  
  # Add the main overlap table
  addWorksheet(wb, "Overlap_Analysis")
  if (nrow(overlap_table) > 0) {
    writeData(wb, "Overlap_Analysis", overlap_table)
  } else {
    writeData(wb, "Overlap_Analysis", data.frame(Message = "No overlapping elements found between datasets"))
  }
  
  # Apply formatting
  headerStyle <- createStyle(
    fontColour = "white", 
    fgFill = "#1F497D",
    halign = "center",
    textDecoration = "bold",
    border = "TopBottomLeftRight"
  )
  
  if (nrow(overlap_table) > 0) {
    addStyle(wb, "Overlap_Analysis", headerStyle, rows = 1, cols = 1:ncol(overlap_table))
    setColWidths(wb, "Overlap_Analysis", cols = 1:ncol(overlap_table), widths = "auto")
  } else {
    addStyle(wb, "Overlap_Analysis", headerStyle, rows = 1, cols = 1)
    setColWidths(wb, "Overlap_Analysis", cols = 1, widths = "auto")
  }
  
  # Add a summary sheet
  addWorksheet(wb, "Summary")
  
  summary_data <- data.frame(
    Metric = c(
      paste("Total unique locus tags in", dataset1_name),
      paste("Total unique locus tags in", dataset2_name),
      paste("Total unique locus tags in", dataset3_name),
      paste("Total unique locus tags in", dataset4_name),
      "Locus tags common to all four datasets",
      paste("Locus tags common to", dataset1_name, dataset2_name, "&", dataset3_name, "only"),
      paste("Locus tags common to", dataset1_name, dataset2_name, "&", dataset4_name, "only"),
      paste("Locus tags common to", dataset1_name, dataset3_name, "&", dataset4_name, "only"),
      paste("Locus tags common to", dataset2_name, dataset3_name, "&", dataset4_name, "only"),
      paste("Locus tags common to", dataset1_name, "&", dataset2_name, "only"),
      paste("Locus tags common to", dataset1_name, "&", dataset3_name, "only"),
      paste("Locus tags common to", dataset1_name, "&", dataset4_name, "only"),
      paste("Locus tags common to", dataset2_name, "&", dataset3_name, "only"),
      paste("Locus tags common to", dataset2_name, "&", dataset4_name, "only"),
      paste("Locus tags common to", dataset3_name, "&", dataset4_name, "only"),
      paste("Unique locus tags only in", dataset1_name),
      paste("Unique locus tags only in", dataset2_name),
      paste("Unique locus tags only in", dataset3_name),
      paste("Unique locus tags only in", dataset4_name)
    ),
    Value = c(
      length(unique(df1$locus_tag)),
      length(unique(df2$locus_tag)),
      length(unique(df3$locus_tag)),
      length(unique(df4$locus_tag)),
      length(common_all_four),
      length(common_df1_df2_df3),
      length(common_df1_df2_df4),
      length(common_df1_df3_df4),
      length(common_df2_df3_df4),
      length(common_df1_df2),
      length(common_df1_df3),
      length(common_df1_df4),
      length(common_df2_df3),
      length(common_df2_df4),
      length(common_df3_df4),
      length(unique_df1),
      length(unique_df2),
      length(unique_df3),
      length(unique_df4)
    )
  )
  
  writeData(wb, "Summary", summary_data)
  addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:2)
  setColWidths(wb, "Summary", cols = 1:2, widths = "auto")
  
  # Add sheets for unique locus tags in each dataset (full row data)
  unique_sheets <- list(
    list(data = unique_df1, df = df1, name = short_name1, full_name = dataset1_name),
    list(data = unique_df2, df = df2, name = short_name2, full_name = dataset2_name),
    list(data = unique_df3, df = df3, name = short_name3, full_name = dataset3_name),
    list(data = unique_df4, df = df4, name = short_name4, full_name = dataset4_name)
  )
  
  for(sheet_info in unique_sheets) {
    if(length(sheet_info$data) > 0) {
      sheet_name <- paste0("Unique_", sheet_info$name)
      sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name limit
      addWorksheet(wb, sheet_name)
      unique_data <- sheet_info$df %>% filter(locus_tag %in% sheet_info$data)
      writeData(wb, sheet_name, unique_data)
      addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(unique_data))
      setColWidths(wb, sheet_name, cols = 1:ncol(unique_data), widths = "auto")
    }
  }
  
  # Save the workbook
  output_file <- paste0("Four_Dataset_Locus_Tag_Overlap_Analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Print confirmation and summary
  cat("✅ Excel file created successfully:", output_file, "\n\n")
  cat("=== FILE CONTENTS ===\n")
  cat("1. 'Overlap_Analysis' sheet: Compact table view with Locus_Tag, Protein, and Label\n")
  cat("2. 'Summary' sheet: Overall statistics\n")
  for(sheet_info in unique_sheets) {
    if(length(sheet_info$data) > 0) {
      cat("3. 'Unique_", sheet_info$name, "' sheet: Locus tags found only in", sheet_info$full_name, "\n")
    }
  }
  cat("\n")
  
  # Show a preview of the main table
  if (nrow(overlap_table) > 0) {
    cat("=== PREVIEW OF OVERLAP TABLE ===\n")
    print(head(overlap_table, 10))
  } else {
    cat("=== OVERLAP TABLE ===\n")
    cat("No overlapping locus tags found between datasets\n")
  }
  
  # Return the results invisibly
  invisible(list(
    venn_plot = venn_plot,
    overlap_table = overlap_table,
    summary = summary_data,
    datasets = list(
      df1 = df1,
      df2 = df2,
      df3 = df3,
      df4 = df4
    ),
    overlaps = list(
      common_all_four = common_all_four,
      common_df1_df2_df3 = common_df1_df2_df3,
      common_df1_df2_df4 = common_df1_df2_df4,
      common_df1_df3_df4 = common_df1_df3_df4,
      common_df2_df3_df4 = common_df2_df3_df4,
      common_df1_df2 = common_df1_df2,
      common_df1_df3 = common_df1_df3,
      common_df1_df4 = common_df1_df4,
      common_df2_df3 = common_df2_df3,
      common_df2_df4 = common_df2_df4,
      common_df3_df4 = common_df3_df4,
      unique_df1 = unique_df1,
      unique_df2 = unique_df2,
      unique_df3 = unique_df3,
      unique_df4 = unique_df4
    )
  ))
}

_______
#Running examples

rm(list = ls())
setwd("path/where/your/files/are")
getwd()

#LcBv vs Bv upregulation comparisons
result <- run_venn_analysis_4way(
  csv_files = c(
    "LcLi_vs_Li_sig_up.csv",
    "LcGv_vs_Gv_sig_up.csv",
    "LcFv_vs_Fv_sig_up.csv",
    "LcPb_vs_Pb_sig_up.csv"
  ),
  colors = c("#B02CBA", "#0121F5", "#4497AF", "#56BC37")  # colours for diagram
)

#LcBv vs Bv downregulation comparisons
result <- run_venn_analysis_4way(
  csv_files = c(
    "LcLi_vs_Li_sig_down.csv",
    "LcGv_vs_Gv_sig_down.csv",
    "LcFv_vs_Fv_sig_down.csv",
    "LcPb_vs_Pb_sig_down.csv"  # Added fourth file
  ),
  colors = c("#e5a6ea", "#91a0ff", "#b0d5e0", "#b8e6aa")  # colours for diagram
)
