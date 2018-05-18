getData <- function(){
  # TCGAbiolinks biblioteka reikalinga duomenų gavimui iš TCGA Legacy
  # duomenų bazės, o data.table funkcija "fread" leidžia nuskaityti failą
  # iš GitHub repozitorijos.
  library(TCGAbiolinks)
  library(data.table)
  
  # Nuskaitomas klinikinių duomenų failas (gautas iš GlioVis duomenų bazės) ir 
  # pašalinami pacientai su rekurentiniu auglio tipu.
  clinical_data <- fread(
    "https://raw.githubusercontent.com/dovydask/gbm-bak/master/gliovis_clinical_data.txt"
    )
  clinical_data <- subset(clinical_data, Recurrence != "Recurrent")
  
  # Formuojama užklausa TCGA Legacy duomenų bazei - reikalaujama
  # gauti visus miRNR ekspresijos duomenis iš multiforminės glioblastomos
  # tyrimo projekto TCGA-GBM.
  query <- GDCquery(project = "TCGA-GBM",
                    data.category = "Gene expression",
                    data.type = "miRNA gene quantification",
                    platform = "H-miRNA_8x15K",
                    legacy = TRUE)
  
  # Atsiunčiami duomenys pagal suformuotą užklausą ir sukuriamas R objektas.
  GDCdownload(query, method = "client")
  data <- GDCprepare(query, directory = "./GDCdata")

  # Duomenų transformacijos patogesnei tolimesnei analizei ir molekulinių
  # potipių iš klinikinių duomenų pridėjimas. Duomenys taip pat normalizuojami
  # ir centruojami.
  data_t <- t(data)
  colnames(data_t) <- data_t[1, ]
  data_t <- data_t[-1, -1]
  cases <- unname(sapply(rownames(data_t), function(x) substr(x, 0, 12)))
  samples <- unname(sapply(rownames(data_t), function(x) substr(x, 0, 28)))
  rownames(data_t) <- samples
  data_t <- cbind(Sample = make.names(cases), data_t)
  data_t <- as.data.frame(data_t)
  data_t <- data_t[order(data_t$Sample), ]
  data_t$Subtype <- clinical_data$Subtype_Verhaak_2010[match(
    data_t$Sample,clinical_data$Sample)]
  data_t <- data_t[complete.cases(data_t), ]
  data_t <- data_t[!duplicated(data_t$Sample), ]
  data_t$Subtype <- factor(data_t$Subtype)
  data_t$Sample <- as.character(data_t$Sample)
  data_t[, -c(1, ncol(data_t))] <- lapply(data_t[, -c(1, ncol(data_t))],
                                          function(x) as.character(x))
  data_t[, -c(1, ncol(data_t))] <- lapply(data_t[, -c(1, ncol(data_t))],
                                          function(x) as.numeric(x))
  data_t[, -c(1, ncol(data_t))] <- scale(data_t[, -c(1, ncol(data_t))])
  
  return(data_t)
}
