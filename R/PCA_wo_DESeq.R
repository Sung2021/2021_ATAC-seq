## Run prcomp for PCA
input.data <- atac.peaks.read
pca.tmp <- prcomp(input.data,scale. = T)
pca.tmp %>% str()
pca.tmp$x
pca.tmp$sdev
## And we can use sdev to compute variance explained by each Principal Component.
var_explained <- pca.tmp$sdev^2/sum(pca.tmp$sdev^2)
pca.tmp$rotation %>% as.data.frame() %>% 
  ggplot(aes(PC1,PC2, shape= substr(rownames(pca.tmp$rotation), 1,6))) + 
  geom_point(size=3) +labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                           y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))
## PC2 and PC3 version
pca.tmp$rotation %>% as.data.frame() %>% 
  ggplot(aes(PC2,PC3, shape= substr(rownames(pca.tmp$rotation), 1,6))) + 
  geom_point(size=3) +labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
                           y=paste0("PC3: ",round(var_explained[3]*100,1),"%"))

pca.tmp$rotation %>% as.data.frame() %>% 
  ggplot(aes(PC1,PC3, shape= substr(rownames(pca.tmp$rotation), 1,6))) + 
  geom_point(size=3) +labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
                           y=paste0("PC3: ",round(var_explained[3]*100,1),"%"))

## sdev gives the standard deviation of principal component. 

## 3D plot of the pca
plotly::plot_ly(pca.tmp$rotation %>% as.data.frame(), x=~PC1,y=~PC2,z=~PC3, 
                color=~substr(rownames(pca.tmp$rotation),1,6))

                
