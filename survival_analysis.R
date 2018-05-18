survivalAnalysis <- function(mirna_data, tp = 0.9, cutoff = 12, rf = FALSE){
  # Reikalingos bibliotekos.
  library(data.table)
  library(survival)
  library(ggkm)
  library(randomForest)
  library(e1071)
  
  cat("Parametrai:\n")
  cat("Treniruojama", tp*100, "proc. visų duomenų.\n")
  cat("Klasifikuojama", ncol(mirna_data)-2, "miRNR duomenimis.\n")
  
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
  
  survival_stats <- clinical_data[match(mirna_data$Sample,
                                        clinical_data$Sample), "survival"]
  survival_stats$survival
  survival_stats <- sort(survival_stats$survival)
  
  # Pacientų išgyvenamumo histograma ir tankio grafikas.
  hist(survival_stats, main = "Išgyvenamumas", xlab = "Mėnesiai",
       ylab = "Pacientų skaičius")
  
  plot(density(survival_stats), main = "Išgyvenamumas",xlab = "Mėnesiai",
       ylab = "Tankis")
  abline(v = density(survival_stats)$x[which.max(density(survival_stats)$y)],
         col = "red", lty = 2)
  abline(v = survival_cutoff, col = "blue", lty = 2)
  density(survival_stats)$x[which.max(density(survival_stats)$y)]
  
  legend(60, 0.04, legend=c(paste(round(
    density(survival_stats)$x[which.max(density(survival_stats)$y)],
    digits = 2), "mėn."), paste(survival_cutoff, "mėn.")),
    col=c("red", "blue"), lty=2, cex=1)
  
  surv_data <- clinical_data[match(
    mirna_data[mirna_data$Subtype == "Classical", "Sample"],
    clinical_data$Sample),
    c("survival", "status", "Subtype_Verhaak_2010")]
  surv_data <- rbind(surv_data, clinical_data[match(
    mirna_data[mirna_data$Subtype == "Mesenchymal", "Sample"],
    clinical_data$Sample),
    c("survival", "status", "Subtype_Verhaak_2010")])
  surv_data <- rbind(surv_data, clinical_data[match(
    mirna_data[mirna_data$Subtype == "Neural", "Sample"],
    clinical_data$Sample),
    c("survival", "status", "Subtype_Verhaak_2010")])
  surv_data <- rbind(surv_data, clinical_data[match(
    mirna_data[mirna_data$Subtype == "Proneural", "Sample"],
    clinical_data$Sample),
    c("survival", "status", "Subtype_Verhaak_2010")])
  
  # Kaplan-Meier potipių išgyvenamumo grafikas.
  fit <- survfit(Surv(survival, status) ~ Subtype_Verhaak_2010,
                 data = surv_data)
  print(ggkm(fit, 
       ystratalabs = c("Klasikinis", "Mezenchiminis",
                       "Nervinis", "Pronervinis"),
       ystrataname = "Molekulinis potipis", legendposition = c(0.75, 0.8),
       ylab = "Išgyvenamumas (%)",
       xlab = "Laikas nuo tyrimo pradžios (mėn.)"
  ))
  
  cat("Nupiešti 3 grafikai.\n")
  
  if(rf){
    # Klasifikavimas į išgyvenamumo grupes.
    cat("\nRandom Forest\n")
    train.prop <- tp
    total <- 0
    for(i in c(1:100)){
      data <- mirna_data_with_survival[sample(
        nrow(mirna_data_with_survival)), -1]
      colnames(data) <- make.names(colnames(data))
      train <- data[1:(train.prop * nrow(data)), ]
      test <- data[(train.prop * nrow(data)+1):nrow(data), ]
      
      train$Subtype <- NULL
      train$CIMP <- NULL
      train$IDH1 <- NULL
      train$MGMT <- NULL
  
      test$Subtype <- NULL
      test$CIMP <- NULL
      test$IDH1 <- NULL
      test$MGMT <- NULL
      
      train$Survived_12_Months <- factor(train$Survived_12_Months)
      test$Survived_12_Months <- factor(test$Survived_12_Months)
      
      model_rf <- randomForest(Survived_12_Months ~ ., data = train)
      preds <- predict(model_rf, test)
      cat(i, "bandymas:",
          sum((unname(preds) == test$Survived_12_Months) == TRUE) /
            length(test$Survived_12_Months), "\n")
      
      total <- total +
        sum((unname(preds) == test$Survived_12_Months) == TRUE) /
        length(test$Survived_12_Months)
    }
    cat("Random Forest tikslumo vidurkis: ", total/100, "\n")
  }
  # Naive Bayes klasifikavimas
  nbc_data <- mirna_data_with_survival[
    , c(colnames(mirna_data_with_survival)
        [2:(ncol(mirna_data_with_survival)-5)], "Survived_12_Months")]
  new_nbc_data <- nbc_data
  model_nbc <- naiveBayes(as.factor(Survived_12_Months) ~ ., data = nbc_data)
  
  cat("\nNaive Bayes apriori tikimybės:\n")
  print(model_nbc$apriori/nrow(nbc_data))
  
  preds_nbc <- predict(model_nbc, 
                       newdata = new_nbc_data[, -ncol(new_nbc_data)])
  
  cat("\nNaive Bayes klaidų matrica:\n")
  table(preds_nbc, new_nbc_data$Survived_12_Months)
  
  cat("\nNaive Bayes tikslumas: ",
      sum((unname(preds_nbc) == new_nbc_data$Survived_12_Months) == TRUE) /
        length(new_nbc_data$Survived_12_Months))
  
  cat("\n\nKolmogorov-Smirnov testai:\n")
  cat("Testas / p-reikšmė / miRNR / Potipis\n")
  
  population_1 <- mirna_data_with_survival[mirna_data_with_survival$Survived_12_Months == "Yes", ]
  population_2 <- mirna_data_with_survival[mirna_data_with_survival$Survived_12_Months == "No", ]

  stopifnot(colnames(population_1) == colnames(population_2))
  mirnas_t_test <- colnames(population_1)[2:(length(colnames(population_1))-5)]
  
  for(subtype in c("Classical", "Mesenchymal", "Neural", "Proneural")){
    for(mirna_t_test in mirnas_t_test){
      output <- ks.test(population_1[population_1$Subtype == subtype, mirna_t_test],
                          population_2[population_2$Subtype == subtype, mirna_t_test])
      
      if(output$p.value < 0.025){
        cat(paste("KS: ", round(output$p.value, digits=5), mirna_t_test, subtype, "\n", sep="\t"))
      }
    }
  }
}