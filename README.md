# Bioinformatikos bakalauro darbas
## Autorius: Dovydas Kičiatovas

Ši GitHub repozitorija skirta mano bioinformatikos bakalauro darbo rezultatų 
generavimui. Nuoroda į darbą: https://bit.ly/2QDoAq3

Šiame archyve darbo rezultatus generuoja R Markdown failas rezultatai.Rmd. Jokių
kitų skriptų vykdyti nereikia. Kaip alternatyva Markdown yra pateiktas failas 
mirna_analysis.R, kuriame galima įvykdyti visas Markdown metu vykdomas 
operacijas.

Kodas parašytas naudojant 1.1.423 versijos RStudio su 3.4.3 versijos R. Visi
papildomi paketai, kurie yra išvardyti rezultatai.Rmd faile, naudojami
naujausių versijų. Išbandyta ant Windows 10 ir Ubuntu 16.04 operacinių sistemų.

SVARBU: jei norite programinį kodą paleisti Unix (Ubuntu) sistemoje, reikia
įrašyti šiuos paketus: libssl-dev, libcurl4-gnutls-dev, libcurl4-openssl-dev,
libmariadb-client-lgpl-dev. Užtenka įvykdyti šias Shell komandas:

sudo apt-get install libssl-dev
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libmariadb-client-lgpl-dev

Šie paketai susiję su iš "trečiosios šalies" (third-party) šaltinių duomenų 
siuntimu, pavyzdžiui, iš TCGA (GDC) Legacy duomenų bazės bei su kai kurių
R paketų įrašymu.
