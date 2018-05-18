hybridClassifier <- function(mirna_data, k = 10, tp = 0.9){
  library(RANN)
  library(randomForest)
  library(e1071)
  
  cat("Parametrai:\n")
  cat("Treniruojama", tp*100, "proc. visų duomenų.\n")
  cat("Klasifikuojama", ncol(mirna_data)-1, "miRNR duomenimis.\n")
  cat("Artimiausiems kaimynams rasti: k =", k, "\n")
  
  train.prop <- tp
  data.complete <- mirna_data
  colnames(data.complete) <- make.names(colnames(data.complete), unique = T)
  
  first_total <- 0
  final_total <- 0
  total_unchanged <- 0
  total_increases <- 0
  cat("\nRandomForest / Hibridinis klas.\n")
  for(n in 1:100){
    data.complete <- data.complete[sample(nrow(data.complete)), ]
    rf_train <- data.complete[1:(train.prop * nrow(data.complete)), ]
    rf_test <- data.complete[(train.prop *
                                nrow(data.complete)+1):nrow(data.complete), ]
    
    rf_model <- randomForest(Subtype ~ ., data = rf_train)
    rf_preds <- predict(rf_model, newdata = rf_test)
    first_result <- sum((unname(rf_preds) == rf_test$Subtype) == TRUE) /
      length(rf_test$Subtype)
    first_total <- first_total + first_result
    
    nb_test <- rf_test
    nb_test$Subtype <- unname(rf_preds)
    
    closest <- nn2(nb_test[, -ncol(nb_test)], searchtype = "standard", k = k)
    closest <- closest$nn.idx
    
    final_guesses <- c()
    for(i in 1:nrow(closest)){
      nearest_data <- nb_test[c(closest[i, ]), ]
      nb_model <- naiveBayes(Subtype ~ ., data = nearest_data)
      nb_preds <- predict(nb_model, newdata = nb_test[i, ])
      guess <- as.character(nb_preds)
      final_guesses <- c(final_guesses, guess)
    }
    
    final_result <- 
      sum((final_guesses == as.character(rf_test$Subtype)) == TRUE) /
      length(rf_test$Subtype)
    final_total <- final_total + final_result
    if(first_result < final_result){
      total_increases <- total_increases + 1
    }
    else if(first_result == final_result){
      total_unchanged <- total_unchanged + 1
    }
    
    cat(n, "bandymas: ", first_result, " ", final_result, "\n")
  }
  cat("\nPagerėjimas / Nėra pokyčio / Pablogėjimas klas. tikslume\n")
  cat(total_increases, " ", total_unchanged, " ",
      100-total_increases-total_unchanged, "\n")
  cat("\nGalutiniai vidurkiai:\n")
  cat(first_total/100, " ", final_total/100, "\n")
}
  
