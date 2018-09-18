variableImportance <- function(rf_model, n){
  # Random Forest bibliotekoje yra ir kintamųjų svarbos funkcijos.
  library(randomForest)
  
  cat("Parametrai:\n")
  cat("Svarbiausių miRNR skaičius:", n, "\n")
  
  # Gaunami nurodyto modelio n kintamųjų svarbos reikšmės.
  importance <- as.data.frame(importance(rf_model))
  importance <- importance[order(-importance$MeanDecreaseGini), ,
                           drop=FALSE]
  
  # Kintamųjų vardai konvertuojami į standartinį pavidalą.
  new_labels <- rownames(importance)[1:n]
  for(x in new_labels){
    if(substr(x, nchar(x), nchar(x)) == "."){
      new <- x
      substr(new, nchar(new), nchar(new)) <- "*"
      new_labels[match(x, new_labels)] <- new
    }
  }
  new_labels <- unname(sapply(new_labels, function(x) gsub("\\.", "-", x)))
  
  # Piešiamas kintamųjų svarbos grafikas.
  varImpPlot(rf_model, n.var = n, main = "Kintamųjų svarba",
             labels = rev(new_labels))
  
  # Modelio rezultatai (tarp jų ir klaidų matrica) išspausdinami.
  print(rf_model)
  
  # Grąžinamas n svarbiausių kintamųjų vardų vektorius.
  return(new_labels)
}