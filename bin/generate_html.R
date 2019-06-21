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
parser$add_argument('--level', required=TRUE, help='2 or 3 level')
args = parser$parse_args()


lane_list <- gsub(".stats.json", "", list.files(args$input_folder, pattern = ".json"))
lane_names <- gsub("L00", "Lane ", lane_list)
lane_nums <- gsub("L00", "", lane_list)

p7_rows <- unlist(stringr::str_split(args$p7_rows, " "))
p5_cols <- unlist(stringr::str_split(args$p5_cols, " "))


for (lane in lane_list) {
  rows <- c("A", "B", "C", "D", "E", "F", "G", "H")
  cols <- 1:12
  plate <- data.frame(expand.grid(cols, rows))
  levels(plate$Var2) <- rev(levels(plate$Var2))
  
  rt_counts <- read.csv(paste0(lane, ".rt_counts.csv"), header = F, stringsAsFactors = FALSE)
  rt_counts$plate <- stringr::str_split_fixed(rt_counts$V1, "-", 2)[,1]
  rt_counts$well <- stringr::str_split_fixed(rt_counts$V1, "-", 2)[,2]
  for(p in unique(rt_counts$plate)) {
    sub_df <- rt_counts[rt_counts$plate == p,]
    rt_df <- data.frame(rows = substring(sub_df$well, first = 1, last = 1),
                        cols = as.numeric(substring(sub_df$well, first = 2,
                                                    last = nchar(sub_df$well))), ReadCount = sub_df$V2)
    
    
    data <- merge(plate, rt_df, by.x=c("Var1", "Var2"),
                  by.y=c("cols", "rows"), all.x=T)
    data$ReadCount[is.na(data$ReadCount)] <- 0
    
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
  pcr_counts <- read.csv(paste0(lane, ".pcr_counts.csv"), header = F, stringsAsFactors = FALSE)
  pcr_counts$p5_row <- substring(pcr_counts$V1, first = 1, last = 1)
  pcr_counts$p5_col <- as.numeric(substring(pcr_counts$V1, first = 2,
                                            last = nchar(pcr_counts$V1)))
  
  pcr_counts$p7_row <- substring(pcr_counts$V2, first = 1, last = 1)
  pcr_counts$p7_col <- as.numeric(substring(pcr_counts$V2, first = 2,
                                            last = nchar(pcr_counts$V2)))
  for (i in  1:length(unlist(stringr::str_split(args$p7_rows, " ")))) {
    sub <- subset(pcr_counts, p7_row == p7_rows[i] & p5_col == as.numeric(p5_cols[i]))
    p5_df <- data.frame(rows = sub$p5_row,
                        cols = sub$p7_col, ReadCount = sub$V3)
    data <- merge(plate, p5_df, by.x=c("Var1", "Var2"),
                  by.y=c("cols", "rows"), all.x=T)
    data$ReadCount[is.na(data$ReadCount)] <- 0
    
    png(file = paste0("demux_dash/img/", lane,"_", p7_rows[i], p5_cols[i], ".pcr_plate.png"), width = 6, height = 4, res = 200, units = "in")
    print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount), data = data) +
            geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
            scale_fill_gradient(low = "white", high = "blue"))
    dev.off()
  }
  p5_agg <- aggregate(pcr_counts$V3, by = list(orig=pcr_counts$V1), sum)
  p5_df <- data.frame(rows = substring(p5_agg$orig, first = 1, last = 1),
                      cols = as.numeric(substring(p5_agg$orig, first = 2,
                                                  last = nchar(p5_agg$orig))), ReadCount = p5_agg$x)
  
  data <- merge(plate, p5_df, by.x=c("Var1", "Var2"),
                by.y=c("cols", "rows"), all.x=T)
  data$ReadCount[is.na(data$ReadCount)] <- 0
  
  
  
  
  
  p7_agg <- aggregate(pcr_counts$V3, by = list(orig=pcr_counts$V2), sum)
  p7_df <- data.frame(rows = substring(p7_agg$orig, first = 1, last = 1),
                      cols = as.numeric(substring(p7_agg$orig, first = 2,
                                                  last = nchar(p7_agg$orig))), ReadCount = p7_agg$x)
  
  data <- merge(plate, p7_df, by.x=c("Var1", "Var2"),
                by.y=c("cols", "rows"), all.x=T)
  data$ReadCount[is.na(data$ReadCount)] <- 0
  
  png(file = paste0("demux_dash/img/", lane, ".p7_plate.png"), width = 6, height = 4, res = 200, units = "in")
  print(ggplot(aes(as.factor(Var1), Var2, fill = ReadCount), data = data) +
          geom_point(shape=21, size = 10) + theme_bw() + labs(x = "", y = "") +
          scale_fill_gradient(low = "white", high = "blue") )
  dev.off()
  
}




outtab <- lapply(lane_list, function(x) {
  sumstats <- jsonlite::fromJSON(paste0(x, ".stats.json"))
  data.frame(`Total input reads` = round(sumstats$total_input_reads),
             `Total passed reads` = round(sumstats$total_passed_reads),
             `Pass percentage` = round(sumstats$fraction_passed_reads * 100, digits = 2),
             `Percent uncorrected`= round(sumstats$fraction_uncorrected_reads * 100, digits = 2),
             `Percent invalid RT`= round(sumstats$fraction_invalid_rt_well * 100, digits = 2),
             `Percent PCR mismatch`= round(sumstats$fraction_pcr_mismatch * 100, digits = 2)
  )
})
tab <- as.data.frame(do.call(rbind, outtab))
row.names(tab) <- lane_names

top <- HTML('<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
    <style>
        .navbar-brand
        {
          font-size: xx-large !important;
          display: flex;
          align-items: center;
          text-align: center;
        }
        .sidebar {
          position: fixed;
          top: 0;
          bottom: 0;
          left: 10px;
          z-index: 100; /* Behind the navbar */
          padding: 100px 0 0; /* Height of navbar */
          box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
        }
        </style>
  </head>')
lig_code <- ""
lig_head <- ""
if(args$level == "3") {
  lig_head <- HTML('<li class="nav-item">
                <a class="nav-link" href="#lig">
                  <span data-feather="shopping-cart"></span>
                  Ligation Barcodes
                </a>
              </li>')
  lig_plates <- unique(lig_counts$plate)
  lig_code <- list(
  HTML('          <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                <h1 class="h3" id="lig">Ligation Barcodes</h1>
            </div>
          <nav>
              <div class="nav nav-tabs" id="navlig-tab" role="tablist">'),
  lapply(lane_nums, function(num) {
    tags$a(class="nav-item nav-link", id=paste0("navlig-lane", num, "-tab"),
           `data-toggle`="tab", href=paste0("#navlig-lane", num), role="tab",
           `aria-controls`=paste0("navlig-lane", num),`aria-selected`="false", paste("Lane", num))
  }),
  HTML('              </div>
          </nav>
       <div class="tab-content" id="nav-tabContent">'),
  lapply(lane_nums, function(num) {
    tags$div(class="tab-pane fade", id=paste0("navlig-lane", num),
             role="tabpanel", `aria-labelledby`=paste0("navlig-lane", num, "-tab"),
             list(lapply(lig_plates, function(i) {
               list(tags$h4(paste0("Ligation plate ",i)),
                    tags$img(src=paste0("img/L00", num, "_", i, ".lig_plate.png"),
                             width = "50%",
                             class="rounded mx-auto d-block",
                             alt="..."))
             })))
    
  }))
}


body <- tags$body(
  list(
    HTML('    <nav class="navbar navbar-expand-md sticky-top navbar-light" style="background-color: #e3f2fd;">
        <div class="navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2">
            <ul class="navbar-nav mr-auto">
                <img src="img/bbi_icon.png" height="70" class="d-inline-block align-top" alt="">
            </ul>
        </div>
        <div class="mx-auto order-0">
            <a class="navbar-brand mx-auto" href="#">Demultiplexing QC Dashboard</a>
        </div>
        <div class="navbar-collapse collapse w-100 order-3 dual-collapse2">
        </div>
    </nav>
    <div class="container-fluid">
      <div class="row">
        <nav class="col-md-2 d-none d-md-block bg-light sidebar">
          <div class="sidebar-sticky">
            <ul class="nav flex-column">
              <li class="nav-item">
                <a class="nav-link active" href="#summary">
                  <span data-feather="home"></span>
                  Summary Statistics <span class="sr-only">(current)</span>
                </a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="#rt">
                  <span data-feather="file"></span>
                  RT Barcodes
                </a>
              </li>
              <li class="nav-item">
                <a class="nav-link" href="#pcr">
                  <span data-feather="shopping-cart"></span>
                  PCR Barcodes
                </a>
              </li>',
         lig_head,
            '</ul>
          </div>
        </nav>
       <main role="main" class="col-md-9 ml-sm-auto col-lg-10 px-4" style="padding-top: 15px;">'),
    
    tags$div(class="tab-content", id="nav-tabContent",
             list(HTML('<div class="tab-content" id="nav-tabContent">
              <div class="tab-pane fade show active" id="navstat-lane1" role="tabpanel" aria-labelledby="navstat-lane1-tab">
                  <h3 class="h3" id="summary">Summary statistics</h3>
                  <table class="table table-hover">
                      <thead>
                        <tr>
                          <th scope="col"></th>'), 
                  lapply(lane_names, function(lane) {
                    tags$th(scope="col", lane)
                  }),
                  HTML('</tr>
                      </thead>
                      <tbody>
                        <tr>
                          <th scope="row">Total input reads</th>'),
                  lapply(tab$Total.input.reads, function(lane) {
                    tags$td(lane)
                  }),
                  HTML('</tr>
                        <tr>
                          <th scope="row">Total passed reads</th>'),
                  lapply(tab$Total.passed.reads, function(lane) {
                    tags$td(lane)
                  }),   
                  HTML('                        </tr>
                        <tr>
                          <th scope="row">Pass percentage</th>'),
                  lapply(tab$Pass.percentage, function(lane) {
                    tags$td(lane)
                  }),
                  HTML('                        </tr>
                        <tr>
                          <th scope="row">Percent uncorrected</th>'),
                  lapply(tab$Percent.uncorrected, function(lane) {
                    tags$td(lane)
                  }),
                  HTML('                        </tr>
                        <tr>
                          <th scope="row">Percent invalid RT</th>'),
                  lapply(tab$Percent.invalid.RT, function(lane) {
                    tags$td(lane)
                  }),
                  HTML('                        </tr>
                        <tr>
                          <th scope="row">Percent PCR mismatch</th>'),
                  lapply(tab$Percent.PCR.mismatch, function(lane) {
                    tags$td(lane)
                  }),
                  HTML('                       </tr>
                      </tbody>
                    </table>'))),
    HTML('          <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                <h1 class="h3" id="pcr">RT Barcodes</h1>
            </div>
          <nav>
              <div class="nav nav-tabs" id="navrt-tab" role="tablist">'),
    lapply(lane_nums, function(num) {
      tags$a(class="nav-item nav-link", id=paste0("navrt-lane", num, "-tab"), 
             `data-toggle`="tab", href=paste0("#navrt-lane", num), role="tab", 
             `aria-controls`=paste0("navrt-lane", num),`aria-selected`="false", paste("Lane", num))
    }),
    HTML('              </div>
          </nav>
       <div class="tab-content" id="nav-tabContent">'),
    lapply(lane_nums, function(num) {
      tags$div(class="tab-pane fade", id=paste0("navrt-lane", num), 
               role="tabpanel", `aria-labelledby`=paste0("navrt-lane", num, "-tab"),
               list(lapply(unique(rt_counts$plate), function(p) {
                 list(tags$h4(paste0("Plate ",p)),
                      tags$img(src=paste0("img/L00", num, "_", p, ".rt_plate.png"), 
                               width = "50%",
                               class="rounded mx-auto d-block",
                               alt="..."))
               })))
      
    }),
    
    # PCR barcodes
    
    HTML('          <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom">
                <h1 class="h3" id="pcr">PCR Barcodes</h1>
            </div>
          <nav>
              <div class="nav nav-tabs" id="navrt-tab" role="tablist">'),
    lapply(lane_nums, function(num) {
      tags$a(class="nav-item nav-link", id=paste0("navpcr-lane", num, "-tab"), 
             `data-toggle`="tab", href=paste0("#navpcr-lane", num), role="tab", 
             `aria-controls`=paste0("navpcr-lane", num),`aria-selected`="false", paste("Lane", num))
    }),
    HTML('              </div>
          </nav>
       <div class="tab-content" id="nav-tabContent">'),
    lapply(lane_nums, function(num) {
      tags$div(class="tab-pane fade", id=paste0("navpcr-lane", num), 
               role="tabpanel", `aria-labelledby`=paste0("navpcr-lane", num, "-tab"),
               list(lapply(1:length(unlist(stringr::str_split(args$p7_rows, " "))), function(i) {
                 list(tags$h4(paste0("PCR Combo ",p7_rows[i], p5_cols[i])),
                      tags$img(src=paste0("img/L00", num, "_", p7_rows[i], p5_cols[i], ".pcr_plate.png"), 
                               width = "50%",
                               class="rounded mx-auto d-block",
                               alt="..."))
               })))
      
    }),
    
    lig_code,
    
    HTML('
    </main>
      </div>
      </div>
      
      
      <!-- Optional JavaScript -->
      <!-- jQuery first, then Popper.js, then Bootstrap JS -->
      <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js" integrity="sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1" crossorigin="anonymous"></script>
      <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js" integrity="sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM" crossorigin="anonymous"></script>
      ')
  )
)

end <- HTML('</html>')


fileConn<-file("demux_dash/demux_dash.html")
writeLines(c(as.character(top),as.character(body), as.character(end)), fileConn)
close(fileConn)


    
