library(ggplot2)
selectSample <- function(data, n){
    set.seed(1)
    v_indices <- sample(1:length(data[,1]), n)
    da <- data[v_indices,]
    return(da)
}

data <- read.table(file = 'PCA_PDB_references.csv', sep = ',', header = TRUE)
da1 <- selectSample(data, 500)
da <- da1

gp_fit_linear <- ggplot(data=da, aes(x=rsr, y=rscc)) + geom_point() + geom_smooth(method="lm")
## gp <- gp_fit_linear + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x="RSR", y="RSCC")
gp_fit_linear <- gp_fit_linear + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RSR", y="RSCC")
gp_fit_linear
## ggsave("rsr_rscc_sample_linear.png")

gp_fit_smooth <- ggplot(data=da, aes(x=rsr, y=rscc)) + geom_point() + geom_smooth()
gp_fit_smooth <- gp_fit_smooth + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RSR", y="RSCC")
gp_fit_smooth
## ggsave("rsr_rscc_sample_smooth.png")

gp_geo_linear <- ggplot(data=da, aes(x=mogul_bonds_rmsz, y=mogul_angles_rmsz)) + geom_point() + geom_smooth(method="lm")
gp_geo_linear <- gp_geo_linear + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RMSZ-bond-length", y="RMSZ-bond-angle")
gp_geo_linear
## ggsave("RMSZ_bond_angel_sample_linear.png")

gp_geo_smooth <- ggplot(data=da, aes(x=mogul_bonds_rmsz, y=mogul_angles_rmsz)) + geom_point() + geom_smooth()
gp_geo_smooth <- gp_geo_smooth + theme(axis.text.x=element_text(size=15), axis.text.y = element_text(size=15)) + labs(x="RMSZ-bond-length", y="RMSZ-bond-angle")
gp_geo_smooth
## ggsave("RMSZ_bond_angel_sample_smooth.png")


