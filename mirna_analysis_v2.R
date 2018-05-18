# Reikalingų bibliotekų sąrašas:
#
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

setwd("C:/Users/Dovydas/Desktop/gbm/R/mirna_analysis")
source("data.R", encoding = "utf-8")
source("linear_regression.R", encoding = "utf-8")
source("classification_analysis.R", encoding = "utf-8")
source("variable_importance.R", encoding = "utf-8")
source("roc.R", encoding = "utf-8")
source("hybrid.R", encoding = "utf-8")
source("survival_analysis.R", encoding = "utf-8")
source("draw_heatmap.R", encoding = "utf-8")

# Duomenų gavimas iš TCGA duomenų bazės. Grąžina miRNR ekspresijų objektą.
# Taip pat iš GitHub repozitorijos nuskaitomas klinikinių duomenų failas.
mirna_data <- getData()

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

# Aprašytas 73 proc. tikslumo Random Forest modelis. Naudojimui su
# variableImportance funkcija.
rf_model <- get(load("./rf_model.RData"))

# Kintamųjų svarba iš pateikto Random Forest modelio. Reikalauja RF modelio
# ir norimo kintamųjų skaičiaus. Grąžina kintamųjų sąrašą ir piešia kintamųjų
# svarbos grafiką.
top_classifier_list <- variableImportance(rf_model, 20)

# Random Forest ROC kreivės piešimas. Reikalauja miRNR ekspresijų objekto
# (be "Sample" identifikatorių) ir duomenų aibės proporcijos treniruotei
# (numatytoji - 0.9). Piešia ROC grafiką.
roc(mirna_data[, c(lm_mirna_list, "Subtype")], 0.8)

# Hibridinis klasifikatorius. Reikalauja Reikalauja miRNR ekspresijų objekto
# (be "Sample" identifikatorių), skaičiaus k (numatytasis - 10) artimiausiems
# kaimynams rasti ir duomenų proporcijos treniruotei (numatytoji - 0,9). 
# Išveda tik tekstinę informaciją.
hybridClassifier(mirna_data[, c(lm_mirna_list, "Subtype")], 7)

# Išgyvenamumo analizė. Reikalauja miRNR ekspresijų objekto (kartu su "Sample"
# identifikatoriais), duomenų proporcijos treniruotei (numatytoji - 0,9), 
# dviejų išgyvenamumo grupių skirstymo riba (numatytoji - 12 mėn.) ir loginis
# indikatorius, žymintis, ar naudoti Random Forest klasifikavimą.
# Iš GitHub repozitorijos nuskaito klinikinių duomenų failą.
# Piešiami 3 grafikai - išgyvenamumo histograma, tankio grafikas ir
# Kaplan-Mieier kreivės.
survivalAnalysis(mirna_data[, c("Sample", lm_mirna_list, "Subtype")])
survivalAnalysis(mirna_data)
survivalAnalysis(mirna_data[, c("Sample", top_classifier_list, "Subtype")],
                 rf = TRUE)

drawHeatmap(mirna_data[, c("Sample", top_classifier_list, "Subtype")])
drawHeatmap(mirna_data[, c("Sample", top_classifier_list, "Subtype")],
            cutoff = 12, surv = TRUE)
