therapyAnalysis <- function(mirna_data, cutoff = 12){
  
  set.seed(123)
  
  # Reikalingos bibliotekos.
  library(data.table)
  library(survival)
  library(ggkm)
  
  mirna_data <- mirna_data[, c("Sample", lm_mirna_list, "Subtype")]
  
  # Klinikinių duomenų nuskaitymas.
  clinical_data <- read.delim("./gliovis_clinical_data.txt")
  clinical_data <- subset(clinical_data, Recurrence != "Recurrent")
  
  # Toliau - įvairios duomenų transformacijos ir grupių sudarymas.
  rad_treatment <- clinical_data[clinical_data$Therapy_Class %in% 
                                   c("Standard Radiation"),
                                 c("Sample", "CIMP_status", "IDH1_status",
                                   "MGMT_status", "Therapy_Class",
                                   "survival")]
  
  rad_treatment_survival <- clinical_data[clinical_data$Therapy_Class %in% 
                                            c("Standard Radiation"),
                                          c("Sample", "survival", "status")]
  
  rad_tmz_treatment <- 
    clinical_data[clinical_data$Therapy_Class == 
                    "Standard Radiation, TMZ Chemo",
                  c("Sample", "CIMP_status", "IDH1_status", "MGMT_status",
                    "Therapy_Class", "survival")]
  
  rad_tmz_treatment_survival <- clinical_data[clinical_data$Therapy_Class %in% 
                                                c("Standard Radiation,
                                                  TMZ Chemo"),
                                           c("Sample", "survival", "status")]
  
  tmz_treatment <- clinical_data[clinical_data$Therapy_Class %in% 
                                   c("TMZ Chemoradiation, TMZ Chemo"),
                              c("Sample", "CIMP_status", "IDH1_status",
                                "MGMT_status", "Therapy_Class", "survival")]
  
  tmz_treatment_survival <- clinical_data[clinical_data$Therapy_Class %in% 
                                            c("TMZ Chemoradiation, TMZ Chemo"),
                                       c("Sample", "survival", "status")]
  
  survival_cutoff <- cutoff
  
  rad_treatment_short <- subset(rad_treatment, survival < survival_cutoff)
  rad_treatment_short <- mirna_data[mirna_data$Sample %in% 
                                      rad_treatment_short$Sample, ]
  rad_treatment_long <- subset(rad_treatment, survival >= survival_cutoff)
  rad_treatment_long <- mirna_data[mirna_data$Sample %in% 
                                     rad_treatment_long$Sample, ]
  rad_treatment_long$Survived_12_Months <- rep("Yes", 
                                               nrow(rad_treatment_long))
  rad_treatment_short$Survived_12_Months <- rep("No", 
                                                nrow(rad_treatment_short))
  rad_treatment_surv <- rbind(rad_treatment_long, rad_treatment_short)
  rad_treatment_surv <- 
    rad_treatment_surv[order(rad_treatment_surv$Subtype, 
                             rad_treatment_surv$Survived_12_Months), ]
  rad_treatment_surv$Survived_12_Months <- 
    factor(rad_treatment_surv$Survived_12_Months)
  rad_treatment_surv$CIMP <- 
    clinical_data[match(rad_treatment_surv$Sample,
                        clinical_data$Sample), "CIMP_status"]
  rad_treatment_surv$IDH1 <- 
    clinical_data[match(rad_treatment_surv$Sample,
                        clinical_data$Sample), "IDH1_status"]
  rad_treatment_surv$MGMT <- 
    clinical_data[match(rad_treatment_surv$Sample,
                        clinical_data$Sample), "MGMT_status"]
  rad_treatment_surv <- 
    rad_treatment_surv[order(rad_treatment_surv$Survived_12_Months,
                             rad_treatment_surv$Subtype, 
                             rad_treatment_surv$MGMT, rad_treatment_surv$IDH1,
                             rad_treatment_surv$CIMP), ]
  
  rad_tmz_treatment_short <- subset(rad_tmz_treatment,
                                    survival < survival_cutoff)
  rad_tmz_treatment_short <- mirna_data[mirna_data$Sample %in% 
                                          rad_tmz_treatment_short$Sample, ]
  rad_tmz_treatment_long <- subset(rad_tmz_treatment, 
                                   survival >= survival_cutoff)
  rad_tmz_treatment_long <- mirna_data[mirna_data$Sample %in% 
                                         rad_tmz_treatment_long$Sample, ]
  rad_tmz_treatment_long$Survived_12_Months <- 
    rep("Yes", nrow(rad_tmz_treatment_long))
  rad_tmz_treatment_short$Survived_12_Months <- 
    rep("No", nrow(rad_tmz_treatment_short))
  rad_tmz_treatment_surv <- rbind(rad_tmz_treatment_long, 
                                  rad_tmz_treatment_short)
  rad_tmz_treatment_surv <- 
    rad_tmz_treatment_surv[order(rad_tmz_treatment_surv$Subtype, 
                                 rad_tmz_treatment_surv$Survived_12_Months), ]
  rad_tmz_treatment_surv$Survived_12_Months <- 
    factor(rad_tmz_treatment_surv$Survived_12_Months)
  rad_tmz_treatment_surv$CIMP <- 
    clinical_data[match(rad_tmz_treatment_surv$Sample, 
                        clinical_data$Sample), "CIMP_status"]
  rad_tmz_treatment_surv$IDH1 <- 
    clinical_data[match(rad_tmz_treatment_surv$Sample, 
                        clinical_data$Sample), "IDH1_status"]
  rad_tmz_treatment_surv$MGMT <- 
    clinical_data[match(rad_tmz_treatment_surv$Sample, 
                        clinical_data$Sample), "MGMT_status"]
  rad_tmz_treatment_surv <- 
    rad_tmz_treatment_surv[order(rad_tmz_treatment_surv$Survived_12_Months, 
                                 rad_tmz_treatment_surv$Subtype, 
                                 rad_tmz_treatment_surv$MGMT, 
                                 rad_tmz_treatment_surv$IDH1, 
                                 rad_tmz_treatment_surv$CIMP), ]
  
  tmz_treatment_short <- subset(tmz_treatment, survival < survival_cutoff)
  tmz_treatment_short <- mirna_data[mirna_data$Sample %in% 
                                      tmz_treatment_short$Sample, ]
  tmz_treatment_long <- subset(tmz_treatment, survival >= survival_cutoff)
  tmz_treatment_long <- mirna_data[mirna_data$Sample %in% 
                                     tmz_treatment_long$Sample, ]
  tmz_treatment_long$Survived_12_Months <- rep("Yes", 
                                               nrow(tmz_treatment_long))
  tmz_treatment_short$Survived_12_Months <- rep("No", 
                                                nrow(tmz_treatment_short))
  tmz_treatment_surv <- rbind(tmz_treatment_long, tmz_treatment_short)
  tmz_treatment_surv <- 
    tmz_treatment_surv[order(tmz_treatment_surv$Subtype,
                             tmz_treatment_surv$Survived_12_Months), ]
  tmz_treatment_surv$Survived_12_Months <- 
    factor(tmz_treatment_surv$Survived_12_Months)
  tmz_treatment_surv$CIMP <- clinical_data[match(tmz_treatment_surv$Sample,
                                                 clinical_data$Sample),
                                           "CIMP_status"]
  tmz_treatment_surv$IDH1 <- clinical_data[match(tmz_treatment_surv$Sample,
                                                 clinical_data$Sample),
                                           "IDH1_status"]
  tmz_treatment_surv$MGMT <- clinical_data[match(tmz_treatment_surv$Sample,
                                                 clinical_data$Sample),
                                           "MGMT_status"]
  tmz_treatment_surv <- 
    tmz_treatment_surv[order(tmz_treatment_surv$Survived_12_Months, 
                             tmz_treatment_surv$Subtype,
                             tmz_treatment_surv$MGMT, tmz_treatment_surv$IDH1,
                             tmz_treatment_surv$CIMP), ]
  
  rad_gbm <- mirna_data[mirna_data$Sample %in% rad_treatment$Sample, ]
  rad_tmz_gbm <- mirna_data[mirna_data$Sample %in% rad_tmz_treatment$Sample, ]
  tmz_gbm <- mirna_data[mirna_data$Sample %in% tmz_treatment$Sample, ]
  
  cat("Stand. radiacija:", nrow(rad_gbm), "pac.\n")
  cat("Stand. radiacija / TMZ:", nrow(rad_tmz_gbm), "pac.\n")
  cat("TMZ:", nrow(tmz_gbm), "pac.\n")
  
  # Išgyvenamumo duomenys.
  surv_data <- clinical_data[match(rad_treatment_surv$Sample, 
                                   clinical_data$Sample), 
                             c("survival", "status", "Therapy_Class")]
  surv_data <- rbind(surv_data, 
                     clinical_data[match(rad_tmz_treatment_surv$Sample,
                                         clinical_data$Sample), 
                                   c("survival", "status", "Therapy_Class")])
  surv_data <- rbind(surv_data, clinical_data[match(tmz_treatment_surv$Sample,
                                                    clinical_data$Sample),
                                              c("survival", "status", 
                                                "Therapy_Class")])
  surv_data$Therapy_Class <- factor(surv_data$Therapy_Class)
  
  # Kaplan-Meier grafikas.
  fit <- survfit(Surv(survival, status) ~ Therapy_Class, data = surv_data)
  ggkm(fit, 
       ystratalabs = c("Standartinė Radiacija", 
                       "Standartinė Radiacija/TMZ Chemoterapija", 
                       "TMZ Chemoradiacija/TMZ Chemoterapija"),
       ystrataname = "Terapijos grupė", legendposition = c(0.75, 0.8),
       ylab = "Išgyvenamumas (%)", xlab = "Laikas nuo tyrimo pradžios (mėn.)"
  )
  
  # Sudaromos populiacijos statistiniams testams.
  population_1 <- mirna_data[mirna_data$Sample %in% rad_gbm$Sample, ]
  population_2 <- mirna_data[mirna_data$Sample %in% rad_tmz_gbm$Sample, ]
  population_2 <- rbind(population_2, mirna_data[mirna_data$Sample %in% 
                                                   tmz_gbm$Sample, ])
  
  cat("\nStd. Rad. vs Rad/TMZ + TMZ\n")
  results_df <- bootstrap(population_1, population_2)
  
  population_1 <- mirna_data[mirna_data$Sample %in% rad_tmz_gbm$Sample, ]
  population_2 <- mirna_data[mirna_data$Sample %in% rad_gbm$Sample, ]
  population_2 <- rbind(population_2, mirna_data[mirna_data$Sample %in% 
                                                   tmz_gbm$Sample, ])
  
  cat("\nRad/TMZ vs. Std. Rad. + TMZ\n")
  results_df <- bootstrap(population_1, population_2)
  
  res <- sort(summary(results_df[results_df$Subtype == "Classical", miRNA]))
  res <- res[length(res):1]
  res <- res[1:20]
  
  population_2 <- mirna_data[mirna_data$Sample %in% tmz_gbm$Sample, ]
  population_1 <- mirna_data[mirna_data$Sample %in% rad_gbm$Sample, ]
  population_1 <- rbind(population_1, mirna_data[mirna_data$Sample %in% 
                                                   rad_tmz_gbm$Sample, ])
  
  cat("\nTMZ vs. Std. Rad. + Rad/TMZ\n")
  results_df <- bootstrap(population_1, population_2)
  
  x <- barplot(res, xaxt = "n")
  text(cex = 0.8, x=x, y=-80, sapply(paste(names(res)), 
                                     function(x) 
                                       substr(x, 5, nchar(x))), xpd=T, srt=90)
}

# Bootstrap funkcija.
bootstrap <- function(population_1, population_2){
  set.seed(123)
  cat("Bootstrapping...\n")
  stopifnot(colnames(population_1) == colnames(population_2))
  subtypes_t_test <- c("Classical", "Mesenchymal", "Neural", "Proneural")
  mirnas_t_test <- colnames(population_1)[2:(length(colnames(population_1))-1)]
  results_df <- data.table(miRNA=character(), Subtype=character(), 
                           pValue=numeric())
  for(subtype_t_test in subtypes_t_test){
    cat(subtype_t_test, "... ")
    pop_1 <- population_1[population_1$Subtype == subtype_t_test, ]
    pop_2 <- population_2[population_2$Subtype == subtype_t_test, ]
    stopifnot(nrow(pop_1) < nrow(pop_2))
    
    for(i in c(1:1000)){
      random_sample <- pop_2[sample(nrow(pop_2), nrow(pop_1)), ]
      for(mirna_t_test in mirnas_t_test){
        output <- wilcox.test(pop_1[, mirna_t_test],
                              random_sample[, mirna_t_test])
        
        if(output$p.value < 0.025){
          #cat(paste("Wilcox: ", round(output$p.value, digits=5), 
          #mirna_t_test, subtype_t_test, "\n", sep="\t"))
          results_df <- rbindlist(list(results_df, 
                                       list(mirna_t_test, subtype_t_test,
                                            round(output$p.value,
                                                  digits = 5))))
        }
      }
    }
    cat("Done.\n")
  }
  cat("\nResults:\n")
  
  results_df$miRNA <- factor(results_df$miRNA)
  
  res <- sort(summary(results_df[results_df$Subtype == "Classical", miRNA]))
  cat("\nClassical\n")
  print(res[res>=500])
  
  res <- sort(summary(results_df[results_df$Subtype == "Mesenchymal", miRNA]))
  cat("\nMesenchymal\n")
  print(res[res>=500])
  
  res <- sort(summary(results_df[results_df$Subtype == "Neural", miRNA]))
  cat("\nNeural\n")
  print(res[res>=500])
  
  res <- sort(summary(results_df[results_df$Subtype == "Proneural", miRNA]))
  cat("\nProneural\n")
  print(res[res>=500])
  
  return(results_df)
}
