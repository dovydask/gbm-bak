classificationAnalysis <- function(mirna_data, tp = 0.9){
  # Random Forest ir Naive Bayes klasifikatorių bibliotekos.
  library(randomForest)
  library(e1071)
  
  cat("Parametrai:\n")
  cat("Treniruojama", tp*100, "proc. visų duomenų.\n")
  cat("Klasifikuojama", ncol(mirna_data)-1, "miRNR duomenimis.\n")
  
  data <- mirna_data
  colnames(data) <- make.names(colnames(data), unique = T)

  # Random Forest klasifikavimas. Nurodoma visos duomenų aibės proporcija
  # treniravimo aibės sudarymui. Klasifikavimas vyksta 100 kartų (su 0.9
  # treniravimo proporcija, vykdoma 10-kartų kryžminė validacija)
  cat("\nRandom Forest\n")
  train.prop <- tp
  total <- 0
  for(i in c(1:100)){
  
    # Duomenys atsitiktinai sumaišomi ir sudaromos treniravimo ir
    # testavimo aibės
    data <- data[sample(nrow(data)), ]
    train <- data[1:(train.prop * nrow(data)), ]
    test <- data[(train.prop * nrow(data)+1):nrow(data), ]
    
    # Random Forest modelis treniruojamas ir testuojamas.
    model_rf <- randomForest(Subtype ~ ., data = train)
    preds <- predict(model_rf, test)

    cat(i, "bandymas:",
        sum((unname(preds) == test$Subtype) == TRUE)/length(test$Subtype),
        "\n")
    
    total <- total + sum((unname(preds) == test$Subtype) == TRUE) / 
      length(test$Subtype)
  }
  
  cat("Random Forest tikslumo vidurkis: ", total/100, "\n")
  
  # Sudaromas Naive Bayes klasifikatoriaus modelis.
  new_data <- data
  model_nbc <- naiveBayes(Subtype ~ ., data = data)
  
  cat("Naive Bayes apriori tikimybės:\n")
  print(model_nbc$apriori/nrow(data))
  
  # Naive Bayes modelis testuojamas ta pačia pilna duomenų aibe,
  # su kuria buvo ir treniruota.
  preds_nbc <- predict(model_nbc, newdata = new_data)
  cat("\nNaive Bayes klaidų matrica:\n")
  print(table(preds_nbc, new_data$Subtype))
  
  cat("\nNaive Bayes tikslumas: ",
      sum((unname(preds_nbc) == new_data$Subtype) == TRUE) /
        length(new_data$Subtype), "\n")
}