roc <- function(mirna_data, tp = 0.9){
  # Reikalingos bibliotekos Naive Bayes ir Random Forest klasifikatoriams.
  library(ROCR)
  library(klaR)
  
  set.seed(123)
  
  z <- (1/(100 - (tp*100)))*100
  
  data <- mirna_data
  colnames(data) <- make.names(colnames(data), unique = T)
  
  subtypes = levels(data$Subtype)
  data <- data[sample(nrow(data)), ]
  
  # Sudaromos treniravimo ir testavimo aibės.
  test_cols = which(1:length(data[, 1]) %% z == 0)
  mirna.test = data[test_cols, ]
  mirna.train = data[-test_cols, ]
  subtype_colors = c("red", "green", "blue", "cyan")
  
  # Tuščias grafikas.
  aucs = c()
  plot(x = NA, y = NA, xlim = c(0, 1), ylim = c(0, 1),
       ylab = "Tikrų pozityvų dažnis (TPR)",
       xlab = "Netikrų pozityvų dažnis (FPR)", bty = "n")
  
  # ROC analizė.
  for(subtype_id in 1:length(subtypes)){
    subtype = as.factor(mirna.train$Subtype == subtypes[subtype_id])
  
    rf_model <- randomForest(subtype ~ .,
                             data = mirna.train[, -ncol(mirna.train)])
    rf_prediction <- predict(rf_model,
                             mirna.test[, -ncol(mirna.test)], type = "prob")
  
    score = rf_prediction[, 'TRUE']
    actual.class = mirna.test$Subtype == subtypes[subtype_id]
    
    pred = prediction(score, actual.class)
    perf = performance(pred, "tpr", "fpr")
    
    roc.x = unlist(perf@x.values)
    roc.y = unlist(perf@y.values)
    lines(roc.y ~ roc.x, col = subtype_colors[subtype_id], lwd = 2)
    
    auc = performance(pred, "auc")
    auc = unlist(slot(auc, "y.values"))
    aucs[subtype_id] = auc
  }
  
  lines(x=c(0,1), c(0,1))
  legend(0.65, 0.4,
         legend=c("Klasikinis", "Mezenchiminis", "Nervinis", "Pronervinis"),
         col=c("red", "green", "blue", "cyan"), lty=1, cex=0.6)
}