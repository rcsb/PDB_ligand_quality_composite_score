library(ggplot2)
selectSample <- function(data, n){
    set.seed(1)
    v_indices <- sample(1:length(data[,1]), n)
    da <- data[v_indices,]
    return(da)
}

data <- read.table(file = 'PCA_PDB_references.csv', sep = ',', header = TRUE)
da <- selectSample(data, 200)
gp_fit_linear <- ggplot(data=da, aes(x=rsr, y=rscc)) + geom_point() + geom_smooth(method="lm")  #linear RSR and RSCC
gp_fit_linear <- gp_fit_linear + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RSR", y="RSCC")
gp_fit_linear

gp_fit_smooth <- ggplot(data=da, aes(x=rsr, y=rscc)) + geom_point() + geom_smooth() #polynomial RSR and RSCC
gp_fit_smooth <- gp_fit_smooth + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RSR", y="RSCC")
gp_fit_smooth

gp_geo_linear <- ggplot(data=da, aes(x=mogul_bonds_rmsz, y=mogul_angles_rmsz)) + geom_point() + geom_smooth(method="lm") #linear RMSZs
gp_geo_linear <- gp_geo_linear + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RMSZ-bond-length", y="RMSZ-bond-angle")
gp_geo_linear

gp_geo_smooth <- ggplot(data=da, aes(x=mogul_bonds_rmsz, y=mogul_angles_rmsz)) + geom_point() + geom_smooth() #polynomial RMSZs
gp_geo_smooth <- gp_geo_smooth + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RMSZ-bond-length", y="RMSZ-bond-angle")
gp_geo_smooth
