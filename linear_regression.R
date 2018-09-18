linearRegression <- function(mirna.data, cutoff = 0.24){
  # Biblioteka potipių perstatoms atlikti.
  library(combinat)
  library(data.table)
  
  # R žodynas potipių perstatoms ir jų vertimu į skaičius.
  dict <- vector(mode="list", length=4)
  names(dict) <- c("Classical", "Mesenchymal", "Neural", "Proneural")
  dict[[1]] <- 1
  dict[[2]] <- 2
  dict[[3]] <- 3
  dict[[4]] <- 4
  
  # Tiesinės regresijos analizė: vykdoma kiekvienai miRNR, visoms 
  # įmanomoms potipių perstatoms. Priimami visi modeliai su
  # nenuliniais koeficientais.
  mirna.lr <- data.table(miRNA=character(), pValue=numeric())
  i <- 1
  subtype_levels <- lapply(levels(factor(mirna.data$Subtype)),
                           function(x) dict[x])
  class_levels <- unname(unlist(subtype_levels))
  classes <- as.numeric(mirna.data$Subtype)
  for (miRNA in names(mirna.data)[2:(length(names(mirna.data))-1)]) {
    i <- i + 1
    if (miRNA != "Subtype") {
      for (permutation in permn(class_levels)) {
        mirna.expressions <- as.numeric(mirna.data[[miRNA]])
        df <- data.frame(classes, mirna.expressions)
        df$classes <- as.numeric(factor(df$classes, levels=permutation))
        df <- df[order(df$classes), ]
        model <- lm(df$mirna.expressions ~ df$classes)
        pvalue <- as.numeric(anova(model)$'Pr(>F)'[1])
        if(abs(model$coefficients[2]) > 0){
          mirna.lr <- rbindlist(list(mirna.lr,
                                     list(miRNA, abs(model$coefficients[2]))))
        }
      }
    }
  }

  # Pasirenkama modelio koeficiento absoliučios reikšmės riba.
  coefficient_cutoff <- cutoff
  
  # Analizės metu gautos miRNR surūšiuojamos ir išmetančios pasikartojančios.
  colnames(mirna.lr) <- c("miRNA", "pValue")
  mirna.ordered <- mirna.lr[order(mirna.lr$miRNA, -mirna.lr$pValue), ]
  mirna.unique <- subset(mirna.ordered, !duplicated(miRNA))
  mirna.filtered <- mirna.unique[mirna.unique$pValue > coefficient_cutoff, ]
  mirna.filtered <- mirna.filtered[order(mirna.filtered$pValue), ]
  mirna.subtyped <- subset(mirna.data,
                           select=c(as.vector(mirna.filtered$miRNA),
                                    "Subtype"))
  mirna.complete <- mirna.subtyped[complete.cases(mirna.subtyped), ]
  sort(colnames(mirna.complete))
  
  return(sort(mirna.filtered$miRNA[1:(nrow(mirna.filtered)-1)]))
}