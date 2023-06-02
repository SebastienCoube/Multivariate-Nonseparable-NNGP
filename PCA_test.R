calidad_aire = read.csv2("../IA-Respiratory/data/calidad_aire/calidad_aire-data_clean.csv")

pollutants = calidad_aire[,  c("PM10","PM2.5","O3","NO2","NO","NOX","SO2","CO"#,"Benzene","Etilbenz","oXylene","Toluene" 
                               )]
pollutants = apply(pollutants, 2, as.numeric)

apply(pollutants, 2, function(x)round(100*sum(is.na(x))/length(x)))

cor(pollutants, use=  "pair")

ACP = FactoMineR::PCA(pollutants, ncp = 6)
FactoMineR::plot.PCA(ACP, choix = "var", axes = c(1,2))
FactoMineR::plot.PCA(ACP, choix = "var", axes = c(3,4))
FactoMineR::plot.PCA(ACP, choix = "var", axes = c(5,6))
hcl = hclust(d = dist(as.matrix(ACP$var$coord)), method = "complete")
plot(hcl)
