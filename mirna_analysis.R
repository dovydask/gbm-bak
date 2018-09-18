# Reikalingų bibliotekų sąrašas:
#
# Bioconductor
# TCGAbiolinks
# data.table
# combinat
# randomForest
# e1071
# ROCR
# klaR
# RANN
# survival
# ggkm
# ComplexHeatmap
# circlize
# githubinstall
# devtools

source("https://bioconductor.org/biocLite.R")
biocLite(ask = FALSE, suppressUpdates = TRUE)

if(!require(ComplexHeatmap)) { biocLite("ComplexHeatmap") }
if(!require(TCGAbiolinks)) { biocLite("TCGAbiolinks") }
if(!require(devtools)){ 
  install.packages("devtools", repos="https://cran.rstudio.com/") }
if(!require(data.table)){ 
  install.packages("data.table", repos="https://cran.rstudio.com/") }
if(!require(combinat)){ 
  install.packages("combinat", repos="https://cran.rstudio.com/") }
if(!require(randomForest)){ 
  install.packages("randomForest", repos="https://cran.rstudio.com/") }
if(!require(e1071)){ 
  install.packages("e1071", repos="https://cran.rstudio.com/") }
if(!require(ROCR)){ 
  install.packages("ROCR",  repos="https://cran.rstudio.com/") }
if(!require(klaR)){ 
  install.packages("klaR", repos="https://cran.rstudio.com/") }
if(!require(RANN)){ 
  install.packages("RANN", repos="https://cran.rstudio.com/") }
if(!require(survival)){ 
  install.packages("survival", repos="https://cran.rstudio.com/") }
if(!require(circlize)){ 
  install.packages("circlize", repos="https://cran.rstudio.com/") }
if(!require(githubinstall)){ 
  install.packages("githubinstall", repos="https://cran.rstudio.com/") }
if(!require(ggkm)){
  library(githubinstall); githubinstall("michealway/ggkm", ask = FALSE) }

source("data.R", encoding = "utf-8")
source("linear_regression.R", encoding = "utf-8")
source("classification_analysis.R", encoding = "utf-8")
source("variable_importance.R", encoding = "utf-8")
source("roc.R", encoding = "utf-8")
source("hybrid.R", encoding = "utf-8")
source("survival_analysis.R", encoding = "utf-8")
source("draw_heatmap.R", encoding = "utf-8")
source("therapy_analysis.R", encoding = "utf-8")

# Duomenų gavimas iš TCGA duomenų bazės. Grąžina miRNR ekspresijų objektą.
# Taip pat iš GitHub repozitorijos nuskaitomas klinikinių duomenų failas.
mirna_data <- getData()

summary(mirna_data$Subtype)

# Žemiau pateikti pavyzdinės miRNR (hsa-miR-139) ekspresijos lygiai
# skirtinguose potipiuose ir atitinkamos miRNR ekspresijų taškinis 
# grafikas kartu su tiesinės regresijos tiese.
model <- lm(as.numeric(mirna_data[["hsa-miR-139"]]) ~ mirna_data$Subtype)
boxplot(as.numeric(mirna_data[["hsa-miR-139"]]) ~ mirna_data$Subtype,
        main = "hsa-miR-139 ekspresija", xaxt = "n",
        xlab = "Potipiai", ylab = "Normalizuota ekspresija")
axis(1, at=1:4, labels=c("Klasikinis", "Mezenchiminis",
                         "Nervinis", "Pronervinis"))

model <- 
  lm(as.numeric(mirna_data[["hsa-miR-139"]]) ~ as.numeric(mirna_data$Subtype))
plot(as.numeric(mirna_data[["hsa-miR-139"]]) ~ as.numeric(mirna_data$Subtype),
     main = "hsa-miR-139 ekspresija", xaxt = "n",
     xlab = "Potipiai", ylab = "Normalizuota ekspresija")
axis(1, at=1:4, labels=c("Klasikinis", "Mezenchiminis",
                         "Nervinis", "Pronervinis"))
abline(model)

model <- lm(as.numeric(mirna_data[["hsa-miR-144"]]) ~ mirna_data$Subtype)
boxplot(as.numeric(mirna_data[["hsa-miR-144"]]) ~ mirna_data$Subtype,
        main = "hsa-miR-144 ekspresija", xaxt = "n",
        xlab = "Potipiai", ylab = "Normalizuota ekspresija")
axis(1, at=1:4, labels=c("Klasikinis", "Mezenchiminis",
                         "Nervinis", "Pronervinis"))

model <-
  lm(as.numeric(mirna_data[["hsa-miR-21"]]) ~ mirna_data$Subtype)
boxplot(as.numeric(mirna_data[["hsa-miR-21"]]) ~ mirna_data$Subtype,
        main = "hsa-miR-21 ekspresija", xaxt = "n",
        xlab = "Potipiai", ylab = "Normalizuota ekspresija")
axis(1, at=1:4, labels=c("Klasikinis", "Mezenchiminis",
                         "Nervinis", "Pronervinis"))

# Tiesinės regresijos analizė. Reikalauja miRNR ekspresijų objekto (gauto 
# iš getData()), tiesės koeficiento absoliučios reikšmės ribos (neprivaloma,
# numatyta reikšmė - 0,24)ir grąžina miRNR sąrašą, sudarytą pagal nustatytą 
# koeficiento ribą (su 0,24 gaunamas 93 miRNR sąrašas).
lm_mirna_list <- linearRegression(mirna_data)

# Klasifikavimo analizė. Reikalauja miRNR ekspresijos duomenų failo
# (pateikiamas miRNR ekspresijų objektas be "Sample" identifikatorių
# stulpelio, su norimų miRNR stulpeliais ir potipiais.) ir duomenų
# proporcijos treniravimo aibei (numatytoji - 0,9). Išveda tik 
# tekstinę informaciją.
classificationAnalysis(mirna_data[, c(lm_mirna_list, "Subtype")])

# Random Forest ROC kreivės piešimas. Reikalauja miRNR ekspresijų objekto
# (be "Sample" identifikatorių) ir duomenų aibės proporcijos treniruotei
# (numatytoji - 0.9). Piešia ROC grafiką.
roc(mirna_data[, c(lm_mirna_list, "Subtype")], 0.8)

# Aprašytas 73 proc. tikslumo Random Forest modelis. Naudojimui su
# variableImportance funkcija.
rf_model <- get(load("./rf_model.RData"))
rf_model

# Kintamųjų svarba iš pateikto Random Forest modelio. Reikalauja RF modelio
# ir norimo kintamųjų skaičiaus. Grąžina kintamųjų sąrašą ir piešia kintamųjų
# svarbos grafiką.
top_classifier_list <- variableImportance(rf_model, 20)

# Intensyvumo žemėlapis (*heatmap*), vaizduojantis iš kintamųjų svarbos gautojo 
# miRNR sąrašo ekspresijos intensyvumą pacientuose. Funkcija reikalauja miRNR 
# ekspresijos duomenų failo (pateikiamas miRNR ekspresijų objektas be "Sample" 
# identifikatorių stulpelio, su norimų miRNR stulpeliais ir potipiais), loginio 
# indikatoriaus, žyminčio, ar skirstyti pacientus pagal išgyvenamumą, kartu 
# vaizduojant ir pacientų CIMP, IDH1 būsenas ir dviejų išgyvenamumo grupių 
# skirstymo ribos (numatytoji - 12 mėn.), jei prieš tai minėtasis loginis 
# indikatorius yra *TRUE*.
drawHeatmap(mirna_data[, c("Sample", top_classifier_list, "Subtype")])

# Hibridinis klasifikatorius. Reikalauja Reikalauja miRNR ekspresijų objekto
# (be "Sample" identifikatorių), skaičiaus k (numatytasis - 10) artimiausiems
# kaimynams rasti ir duomenų proporcijos treniruotei (numatytoji - 0,9). 
# Išveda tik tekstinę informaciją.
hybridClassifier(mirna_data[, c(lm_mirna_list, "Subtype")], 12)

# Išgyvenamumo analizė. Reikalauja miRNR ekspresijų objekto (kartu su "Sample"
# identifikatoriais), duomenų proporcijos treniruotei (numatytoji - 0,9), 
# dviejų išgyvenamumo grupių skirstymo riba (numatytoji - 12 mėn.) ir loginis
# indikatorius, žymintis, ar naudoti Random Forest klasifikavimą.
# Iš GitHub repozitorijos nuskaito klinikinių duomenų failą.
# Piešiami 3 grafikai - išgyvenamumo histograma, tankio grafikas ir
# Kaplan-Mieier kreivės.
survivalAnalysis(mirna_data[, c("Sample", lm_mirna_list, "Subtype")])

# Terapijų grupių analizė. Funkcija reikalauja miRNR ekspresijų objekto 
# (kartu su "Sample" identifikatoriais), dviejų išgyvenamumo grupių skirstymo 
# ribos (numatytoji - 12 mėn.).
therapyAnalysis(mirna_data[, c("Sample", lm_mirna_list, "Subtype")])

# Pavaizduojamas miRNR intensyvumo žemėlapis kartu su pacientų išgyvenamumu ir 
# CIMP bei IDH1 būsenomis.
drawHeatmap(mirna_data[, c("Sample", top_classifier_list, "Subtype")],
            cutoff = 12, surv = TRUE)
