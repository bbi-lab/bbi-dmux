#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(argparse)
  library(jsonlite)
  library(shiny)
  library(stringr)
})

parser = argparse::ArgumentParser(description='Script to generate demux dashboard.')
parser$add_argument('input_folder', help='input folder.')
parser$add_argument('--p7_rows', required=TRUE, help='p7 rows')
parser$add_argument('--p5_cols', required=TRUE, help='p5 cols')
parser$add_argument('--p7_wells', required=TRUE, help='p7 wells')
parser$add_argument('--p5_wells', required=TRUE, help='p5 wells')
parser$add_argument('--level', required=TRUE, help='2 or 3 level')
parser$add_argument('--project_name', required=TRUE, help='Name of the project')
parser$add_argument('--sample_sheet', required=TRUE, help='sample_sheet')
args = parser$parse_args()

lane_list <- gsub(".stats.json", "", list.files(args$input_folder, pattern = ".json"))
lane_names <- gsub("L00", "Lane ", lane_list)
lane_nums <- gsub("L00", "", lane_list)

well_fix <- list("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

if (args$p7_rows != 0) {
  if (args$p7_rows != "none") {
    p7_rows <- unlist(stringr::str_split(args$p7_rows, " "))
  } else {
    p7_rows <- "none"
  } 
  if (args$p5_cols != "none") {
    p5_cols <- unlist(stringr::str_split(args$p5_cols, " "))
  } else {
    p5_cols <- "none"
  }
  
  well_wise <- FALSE
} else {
  if (args$p7_wells != "none") {
    p7_rows <- unlist(stringr::str_split(args$p7_wells, " "))
    p7_rows <- paste0(substr(p7_rows, start=1, stop=1), well_fix[as.numeric(substr(p7_rows, start = 2, stop=3))])
  } else {
    p7_rows <- "none"
  }
  if (args$p5_wells != "none") {
    p5_cols <- unlist(stringr::str_split(args$p5_wells, " "))
    p5_cols <- paste0(substr(p5_cols, start=1, stop=1), well_fix[as.numeric(substr(p5_cols, start = 2, stop=3))])
  } else {
    p5_cols <- "none"
  }
  
  well_df <- data.frame(p5 = p5_cols, p7=p7_rows)
  well_wise <- TRUE
  
}

samp <- read.csv(args$sample_sheet, header = T)

if(length(samp$Sample.ID == "Sentinel") > 0) {
  sent_barcs <- as.character(samp[samp$Sample.ID == "Sentinel",]$RT.Barcode)
} else {
  sent_barcs <- NULL
}
if(length(samp$Sample.ID == "Barnyard") > 0) {
  barn_barcs <- as.character(samp[samp$Sample.ID == "Barnyard",]$RT.Barcode)
} else {
  barn_barcs <- NULL
}

barn_norm_pics <- c()
sent_norm_pics <- c()
bad_well_list <- list()
barcode_list <- c()
well_df_list <- list()
include_norm = "false"

for (lane in lane_list) {
  rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
  cols <- 1:12
  plate <- data.frame(expand.grid(cols, rows))
  levels(plate$Var2) <- rev(levels(plate$Var2))
  
  rt_counts <- read.csv(paste0(lane, ".rt_counts.csv"), header = F, stringsAsFactors = FALSE)
  rt_counts$plate <- stringr::str_split_fixed(rt_counts$V1, "-", 2)[,1]
  rt_counts$well <- stringr::str_split_fixed(rt_counts$V1, "-", 2)[,2]
  
  rt_counts$rows <- substring(rt_counts$well, first = 1, last = 1)
  rt_counts$cols <- suppressWarnings(as.numeric(substring(rt_counts$well, first = 2,
                                                          last = nchar(rt_counts$well))))
  rt_counts$ReadCount <- rt_counts$V2
  
  bad_well <- rt_counts[is.na(rt_counts$rows) | is.na(rt_counts$cols) | !rt_counts$rows %in% rows | !rt_counts$cols %in% cols,]
  rt_counts <- rt_counts[!(is.na(rt_counts$rows) | is.na(rt_counts$cols) | !rt_counts$rows %in% rows | !rt_counts$cols %in% cols),]
  
  for(p in unique(rt_counts$plate)) {
    sub_df <- rt_counts[rt_counts$plate == p,]
    rt_df <- sub_df[,c("rows", "cols", "ReadCount")]
    
    data <- merge(plate, rt_df, by.x=c("Var1", "Var2"),
                  by.y=c("cols", "rows"), all.x=T)
    data$ReadCount[is.na(data$ReadCount)] <- 0
    sent_norm <- NULL
    if (!is.null(sent_barcs)) {
      if (any(sent_barcs %in% sub_df$V1)) {
        sent_norm <- mean(sub_df[sub_df$V1 %in% sent_barcs,]$V2)
      }
    }
    barn_norm <- NULL
    if (!is.null(barn_barcs)) {
      if (any(barn_barcs %in% sub_df$V1)) {
        barn_norm <- mean(sub_df[sub_df$V1 %in% barn_barcs,]$V2)
      }
    }
    
    if (!is.null(sent_norm)) {
      png(file = paste0("demux_dash/img/", lane, "_", p, ".rt_plate_sent_norm.png"), width = 6, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2, fill = log2(ReadCount/sent_norm)), data = data) +
              geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
              scale_fill_gradient2(name = "log2 norm value", low = "red", mid="white", high = "blue"))
      dev.off()
      
      sent_norm_pics <- c(sent_norm_pics, paste0("demux_dash/img/", lane, "_", p, ".rt_plate_sent_norm.png"))
      include_norm <- "true"
    } else {
      png(file = paste0("demux_dash/img/", lane, "_", p, ".rt_plate_sent_norm.png"), width = 5.5, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2), data = data) +
              geom_point(shape=21, size = 10) + geom_text(aes(x = 6.5, y = "D", label = "No Sentinel detected for this plate")) +
              theme_bw() + labs(x = "", y = "") )
      dev.off()
      
      sent_norm_pics <- c(sent_norm_pics, paste0("demux_dash/img/", lane, "_", p, ".rt_plate_sent_norm.png"))
    }
    
    
    if (!is.null(barn_norm)) {
      png(file = paste0("demux_dash/img/", lane, "_", p, ".rt_plate_barn_norm.png"), width = 6, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2, fill = log2(ReadCount/barn_norm)), data = data) +
              geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
              scale_fill_gradient2(name = "log2 norm value", low = "red", mid="white", high = "blue"))
      dev.off()
      barn_norm_pics <- c(barn_norm_pics, paste0("demux_dash/img/", lane, "_", p, ".rt_plate_barn_norm.png"))
      include_norm <- "true"
    } else {
      png(file = paste0("demux_dash/img/", lane, "_", p, ".rt_plate_barn_norm.png"), width = 5.5, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2), data = data) +
              geom_point(shape=21, size = 10) + geom_text(aes(x = 6.5, y = "D", label = "No barnyard detected for this plate")) +
              theme_bw() + labs(x = "", y = "") )
      dev.off()
      
      barn_norm_pics <- c(barn_norm_pics, paste0("demux_dash/img/", lane, "_", p, ".rt_plate_barn_norm.png"))
      
    }
    
    png(file = paste0("demux_dash/img/", lane, "_", p, ".rt_plate.png"), width = 6, height = 4, res = 200, units = "in")
    print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount), data = data) +
            geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
            scale_fill_gradient(low = "white", high = "blue"))
    dev.off()
  }
  
  if(args$level == "3") {
    rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
    cols <- 1:12
    plate <- data.frame(expand.grid(cols, rows))
    lig_map <- do.call(rbind, list(plate, plate, plate, plate))
    lig_map$lig <- paste0("LIG", 1:384)
    lig_map$plate <- c(rep("P1", 96), rep("P2", 96), rep("P3", 96), rep("P4", 96))
    levels(plate$Var2) <- rev(levels(plate$Var2))
    
    
    lig_counts <- read.csv(paste0(lane, ".lig_counts.csv"), header = F, stringsAsFactors = FALSE)
    lig_counts <- merge(lig_counts, lig_map, by.x="V1", by.y="lig", all.x=T)
    lig_counts$well <- paste0(lig_counts$Var2, lig_counts$Var1)
    for(p in unique(lig_counts$plate)) {
      sub_df <- lig_counts[lig_counts$plate == p,]
      lig_df <- data.frame(rows = substring(sub_df$well, first = 1, last = 1),
                           cols = as.numeric(substring(sub_df$well, first = 2,
                                                       last = nchar(sub_df$well))), ReadCount = sub_df$V2)
      
      
      data <- merge(plate, lig_df, by.x=c("Var1", "Var2"),
                    by.y=c("cols", "rows"), all.x=T)
      data$ReadCount[is.na(data$ReadCount)] <- 0
      
      png(file = paste0("demux_dash/img/", lane, "_", p, ".lig_plate.png"), width = 6, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount), data = data) +
              geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
              scale_fill_gradient(low = "white", high = "blue"))
      dev.off()
    }
  }
  if (nrow(bad_well) > 0) {
    barcode_list <- union(barcode_list, as.character(bad_well$V1))
    temp <- as.list(bad_well$V2)
    names(temp) <- bad_well$V1
    bad_well_list[[length(bad_well_list) + 1]] <- temp
    names(bad_well_list)[[length(bad_well_list)]] <- gsub("L00", "Lane ", lane)
  }
 
  pcr_counts <- read.csv(paste0(lane, ".pcr_counts.csv"), header = F, stringsAsFactors = FALSE)
  if(p5_cols[1] != "none") {
    pcr_counts$p5_row <- substring(pcr_counts$V1, first = 1, last = 1)
    pcr_counts$p5_col <- as.numeric(substring(pcr_counts$V1, first = 2,
                                              last = nchar(pcr_counts$V1)))
  }
  
  if (p7_rows[1] != "none") {
    pcr_counts$p7_row <- substring(pcr_counts$V2, first = 1, last = 1)
    pcr_counts$p7_col <- as.numeric(substring(pcr_counts$V2, first = 2,
                                              last = nchar(pcr_counts$V2)))
  }
  if (p7_rows[1] == "none" || p5_cols[1] == "none") {
    pcr_plate_list <- c()
    rel_barc <- ifelse(p7_rows[1] == "none", "p5", "p7")
    if(p7_rows[1] == "none") {
      pcr_counts$rows <- substring(pcr_counts$V1, first = 1, last = 1)
      pcr_counts$cols <- as.numeric(substring(pcr_counts$V1, first = 2,
                                              last = nchar(pcr_counts$V1)))
    } else {
      pcr_counts$rows <- substring(pcr_counts$V2, first = 1, last = 1)
      pcr_counts$cols <- as.numeric(substring(pcr_counts$V2, first = 2,
                                              last = nchar(pcr_counts$V2)))
    }
    p5_df <- data.frame(rows = pcr_counts$rows,
                        cols = pcr_counts$cols, ReadCount = pcr_counts$V3)
    data <- merge(plate, p5_df, by.x=c("Var1", "Var2"),
                  by.y=c("cols", "rows"), all.x=T)
    data$ReadCount[is.na(data$ReadCount)] <- 0
    
    data$outlier <- data$ReadCount < 0.01 * median(data$ReadCount)
    pcr_plate_list <- c(pcr_plate_list, paste0(rel_barc))
    png(file = paste0("demux_dash/img/", lane,"_",rel_barc, ".pcr_plate.png"), width = 6, height = 4, res = 200, units = "in")
    print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount, color=outlier, stroke=outlier), data = data) +
            geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
            scale_fill_gradient(low = "white", high = "blue") + scale_color_manual(values=c("black", "red"), guide=FALSE) +
            scale_discrete_manual(aesthetics = "stroke", values = c(0.5,1), guide=FALSE
            ))
    dev.off()
    well_df_lane <- data.frame()
  } else if (!well_wise) {
    pcr_plate_list <- c()
    for (i in  1:length(p7_rows)) {
      sub <- subset(pcr_counts, p7_row == p7_rows[i] & p5_col == as.numeric(p5_cols[i]))
      p5_df <- data.frame(rows = sub$p5_row,
                          cols = sub$p7_col, ReadCount = sub$V3)
      data <- merge(plate, p5_df, by.x=c("Var1", "Var2"),
                    by.y=c("cols", "rows"), all.x=T)
      data$ReadCount[is.na(data$ReadCount)] <- 0
      
      data$outlier <- data$ReadCount < 0.01 * median(data$ReadCount)
      pcr_plate_list <- c(pcr_plate_list, paste0(p7_rows[i], p5_cols[i]))
      png(file = paste0("demux_dash/img/", lane,"_", p7_rows[i], p5_cols[i], ".pcr_plate.png"), width = 6, height = 4, res = 200, units = "in")
      print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount, color=outlier, stroke=outlier), data = data) +
              geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
              scale_fill_gradient(low = "white", high = "blue") + scale_color_manual(values=c("black", "red"), guide=FALSE) +
              scale_discrete_manual(aesthetics = "stroke", values = c(0.5,1), guide=FALSE
              ))
      dev.off()
      well_df_lane <- data.frame()
    }
    
  } else {
    pcr_plate_list <- c()
    well_df_lane <- merge(well_df, pcr_counts[,c("V1", "V2", "V3")], by.x=c("p5", "p7"), by.y=c("V1", "V2"))
    names(well_df_lane)[3] <- gsub("L00", "Lane ", lane)
  }
  well_df_list[[lane]] <- well_df_lane
}

outtab <- lapply(lane_list, function(x) {
  sumstats <- jsonlite::fromJSON(paste0(x, ".stats.json"))
  list("Lane" = gsub("L00", "Lane ", x),
       "tot_inp_reads" = round(sumstats$total_input_reads),
       "tot_pass_reads" = round(sumstats$total_passed_reads),
       "pass_perc" = round(sumstats$fraction_passed_reads * 100, digits = 2),
       "perc_uncorr"= round(sumstats$fraction_uncorrected_reads * 100, digits = 2),
       "perc_inval_rt"= round(sumstats$fraction_invalid_rt_well * 100, digits = 2),
       "perc_pcr_mismatch"= round(sumstats$fraction_pcr_mismatch * 100, digits = 2)
  )
})

if(well_wise) {
  well_df_out <- tryCatch(
      error = function(condition) {
        message("generate_html.R: unable to make pcr_well_info for the demux_dash:\n  ", condition, appendLF=FALSE)
        message("This does not affect the analysis result.")
        return(well_df_out=data.frame())
      },
      Reduce(function(df1, df2) merge(df1, df2, by = c("p5", "p7")), well_df_list)
  )
} else {
  well_df_out <- data.frame()
}

json_info <- list("run_name" = args$project_name,
                  "lane_list" = lane_nums,
                  "plate_list" = unique(rt_counts$plate),
                  "pcr_combo_list" = unique(pcr_plate_list),
                  "lig_combo_list" = c("P1", "P2", "P3", "P4"),
                  "level" = args$level,
                  "lane_stats" = outtab,
                  "include_norm" = include_norm,
                  "bad_wells_barcodes" = barcode_list, 
                  "bad_wells" = bad_well_list,
                  "pcr_well_info" = well_df_out)

fileConn<-file("demux_dash/js/run_data.js")
writeLines(c("const run_data =", toJSON(json_info, pretty=TRUE)), fileConn)
close(fileConn)

