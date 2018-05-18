drawHeatmap <- function(mirna_data, cutoff, surv = FALSE){
  
  library(ComplexHeatmap)
  library(circlize)
  
  if(surv){
    # Nuskaitomas klinikinių duomenų failas.
    clinical_data <- fread(
      "https://raw.githubusercontent.com/dovydask/gbm-bak/master/gliovis_clinical_data.txt"
    )
    clinical_data <- subset(clinical_data, Recurrence != "Recurrent")

    # Pasirinkta grupių skirstymo riba
    survival_cutoff <- cutoff
    
    # Toliau vykdoma įvairi tarpinė duomenų transformacija.
    long_term_survival <- clinical_data[clinical_data$survival >= 
                                          survival_cutoff,
                                        c("Sample", "CIMP_status",
                                          "IDH1_status", "MGMT_status",
                                          "Therapy_Class", "survival")]
    long_term_survival <- mirna_data[mirna_data$Sample %in% 
                                       long_term_survival$Sample, ]
    long_term_survival["CIMP"] <- clinical_data[match(long_term_survival$Sample,
                                                      clinical_data$Sample),
                                                "CIMP_status"]
    long_term_survival["IDH1"] <- clinical_data[match(long_term_survival$Sample,
                                                      clinical_data$Sample),
                                                "IDH1_status"]
    long_term_survival["MGMT"] <- clinical_data[match(long_term_survival$Sample,
                                                      clinical_data$Sample),
                                                "MGMT_status"]
    long_term_survival$Survived_12_Months <- rep("Yes", nrow(long_term_survival))
    
    short_term_survival <- clinical_data[clinical_data$survival < 
                                           survival_cutoff,
                                         c("Sample", "CIMP_status",
                                           "IDH1_status", "MGMT_status",
                                           "Therapy_Class", "survival")]
    short_term_survival <- mirna_data[mirna_data$Sample %in%
                                        short_term_survival$Sample, ]
    short_term_survival["CIMP"] <- clinical_data[match(short_term_survival$Sample,
                                                       clinical_data$Sample),
                                                 "CIMP_status"]
    short_term_survival["IDH1"] <- clinical_data[match(short_term_survival$Sample,
                                                       clinical_data$Sample),
                                                 "IDH1_status"]
    short_term_survival["MGMT"] <- clinical_data[match(short_term_survival$Sample,
                                                       clinical_data$Sample),
                                                 "MGMT_status"]
    short_term_survival$Survived_12_Months <- rep("No", nrow(short_term_survival))
    
    cat("\n1 grupė: ", nrow(short_term_survival), "pacientai.\n")
    cat("2 grupė: ", nrow(long_term_survival), "pacientai.\n")
    
    mirna_data_with_survival <- rbind(long_term_survival, short_term_survival)
    mirna_data_with_survival <- mirna_data_with_survival[order(
      mirna_data_with_survival$Subtype,
      mirna_data_with_survival$Survived_12_Months), ]
    
    mirna_data_with_survival[, 2:(ncol(mirna_data_with_survival)-5)] <- 
      lapply(mirna_data_with_survival[, 2:(ncol(mirna_data_with_survival)-5)],
             function(x) as.numeric(x))
    
    mirna_data_with_survival <- mirna_data_with_survival[order(
      mirna_data_with_survival$Survived_12_Months, 
      mirna_data_with_survival$Subtype, 
      mirna_data_with_survival$MGMT, 
      mirna_data_with_survival$IDH1, 
      mirna_data_with_survival$CIMP), ]
    
    data <- mirna_data_with_survival[order(mirna_data_with_survival$CIMP,
                                           mirna_data_with_survival$IDH1), -1]
    colnames(data) <- make.names(colnames(data))
    colnames(data)[ncol(data)-4] <- "Potipis"
    data$CIMP <- factor(data$CIMP)
    data$IDH1 <- factor(data$IDH1)
    data$Potipis <- factor(data$Potipis)
    levels(data$Potipis) <- c("Klasikinis", "Mezenchiminis",
                              "Nervinis", "Pronervinis")
    levels(data$CIMP) <- c("G-CIMP", "Ne G-CIMP")
    levels(data$IDH1) <- c("Mutavęs", "Laukinio tipo")
    plot_data <- scale(as.matrix(data[, 1:(ncol(data)-5)]))
    plot_data_t <- t(plot_data)
    
    new_labels <- colnames(data)[1:(ncol(data)-5)]
    for(x in new_labels){
      if(substr(x, nchar(x), nchar(x)) == "."){
        new <- x
        substr(new, nchar(new), nchar(new)) <- "*"
        new_labels[match(x, new_labels)] <- new
      }
    }
    new_labels <- unname(sapply(new_labels, function(x) gsub("\\.", "-", x)))
    
    ha <- HeatmapAnnotation(df = data.frame(Potipis = data$Potipis,
                                            CIMP = data$CIMP,
                                            IDH1 = data$IDH1),
                            col = list(
                              Potipis = c("Klasikinis" = "red",
                                          "Mezenchiminis" = "forestgreen",
                                          "Nervinis" = "blue",
                                          "Pronervinis" = "yellow"),
                                       CIMP = c("G-CIMP" = "firebrick",
                                                "Ne G-CIMP" = "blue4",
                                                "NA" = "white"),
                                       IDH1 = c("Laukinio tipo" = "green",
                                                "Mutavęs" = "deepskyblue")))
    
    ht <- Heatmap(plot_data_t, cluster_columns = F, 
                  show_column_dend = FALSE, show_column_names = FALSE,
                  show_heatmap_legend = FALSE, show_row_names = FALSE,
                  km = 1, top_annotation = ha, name = "miRNR ekspresija",
                  column_title = "miRNR ekspresija")
    
    max_text_width(1)
    ht_list <- ht + rowAnnotation(labels = anno_text(new_labels,
                                                     which = "row",
                                                     just = "right"),
                                  width = max(grobWidth(textGrob(new_labels))))
    
    draw(ht_list, gap = unit(2, "cm"), annotation_legend_side = "bottom",
         heatmap_legend_side = "left")
    
  }
  else{
    data <- mirna_data[, -1]
    data <- data[order(data$Subtype), ]
    colnames(data)[ncol(data)] <- "Potipis"
    data$Potipis <- factor(data$Potipis)
    levels(data$Potipis) <- c("Klasikinis", "Mezenchiminis",
                              "Nervinis", "Pronervinis")
    plot_data <- scale(as.matrix(data[, 1:(ncol(data)-1)]))
    plot_data_t <- t(plot_data)
    
    new_labels <- colnames(data)[1:(ncol(data)-1)]
    for(x in new_labels){
      if(substr(x, nchar(x), nchar(x)) == "."){
        new <- x
        substr(new, nchar(new), nchar(new)) <- "*"
        new_labels[match(x, new_labels)] <- new
      }
    }
    
    ha <- HeatmapAnnotation(df = data.frame(type = data$Potipis),
                            col = list(type = 
                                         c("Klasikinis" = "red",
                                           "Mezenchiminis" = "forestgreen",
                                           "Nervinis" = "blue",
                                           "Pronervinis" = "yellow")),
                            annotation_legend_param = list(title = "Potipis")
    )
    
    ht <- Heatmap(plot_data_t, cluster_columns = F, 
                   show_column_dend = FALSE, show_column_names = FALSE,
                  show_row_names = FALSE, show_heatmap_legend = FALSE,
                   km = 4, km_title = "%i Klasteris", top_annotation = ha,
                  name = "miRNR ekspresija", column_title = "miRNR ekspresija"
    )
    
    ht_list <- ht + rowAnnotation(labels = anno_text(new_labels,
                                                     which = "row",
                                                     just = "right"),
                                  width = max(grobWidth(textGrob(new_labels))))
    
    draw(ht_list, gap = unit(2, "cm"), annotation_legend_side = "bottom",
         heatmap_legend_side = "left")
    
  }
}