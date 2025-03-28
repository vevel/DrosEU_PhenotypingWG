rm(list=ls(all=TRUE))
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggforce)
library(cowplot)

master <- read.csv("/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/all_models_line_meta_compound_random_coefs_wide_220916.csv")

# Focusing on analysing the sexes separately. For males this is easy, but there are still a few different ways of 
# doing this (at least for females):
# These will be the PCAs analysed:
# M9.       9 male traits:
#           "CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M"
# F9.       9 female traits corresponding to these male traits:
#           "CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F"
# Fmax.     All 12 possible female traits:
#           "CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F"
# Fmaxplus. All 12 possible female traits plus viability:
#           "CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "Via_NA"

# Make variables for Country and the sexes
Country <- c(rep("Turkey", 20),rep("Portugal", 17),rep("Spain", 15),rep("Germany", 20),rep("Austria", 20),rep("Ukraine", 19),rep("Denmark", 20),rep("Russia", 20),rep("Finland", 22))
SexM <-c(rep("M", 173))
SexF <-c(rep("F", 173))

####
#### Formatting the data for the different versions of the PCA
####
#################################
# M9
masterM9 <- master[,c("Population","Line","CCRT_M","CSM_M","DT_A_M","DW_M","HSM_M","LS_M","SR_M","TL_M","WA_L_M")]
data1 <- na.omit(cbind(Country,masterM9))
# F9 
masterF9 <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","DW_F","HSM_F","LS_F","SR_F","TL_F","WA_L_F")]
data2 <- na.omit(cbind(Country,masterF9))
# Fmax
masterFmax <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F")]
data3 <- na.omit(cbind(Country,masterFmax))
# Fmaxplus
masterFmaxP <- master[,c("Population","Line","CCRT_F","CSM_F","DT_A_F","Dia_F","DW_F","Fec_F","HSM_F","LS_F","Pgm_Total_F","SR_F","TL_F","WA_L_F", "Via_NA")]
data4 <- na.omit(cbind(Country,masterFmaxP))

#####################
# Colour palette
#    Country   Color
#1   Austria #E38800
#2   Germany #F6C200
#4    Russia #095888
#7   Finland #A00E00
#8   Ukraine #0086A8
#13  Denmark #D04E00
#16 Portugal #7BA354
#18   Turkey #132B69
#71    Spain #B82E00

palette =c("#E38800","#D04E00","#A00E00","#F6C200","#7BA354","#095888","#B82E00","#132B69","#0086A8")

#####################
# PCA 1 = M9 (columns 1-3 are qualitative)
pca1<-PCA(data1,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca1$eig
pca1$var$coord
# Here are the individual coordinates for PC1 output with the qualitative vars
PCA1_ind<-cbind(data1[,c(1:3)], pca1$ind$coord)
p1<-plot(pca1,choix="var",axes=c(1,2),cex=1.4)
q1<-plot(pca1,choix="var",axes=c(2,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca1, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca1, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
p1a<-fviz(pca1, title = "Male PCA - 9 traits (M9) PC1 vs PC2",
            element = "ind",
                     habillage =  data1$Country,
                     geom = c("point","text"),
                     label= "quali",
                     pointsize = 2,
                     pointshape = 18,
                     palette = palette,
                     addEllipses = TRUE, # Concentration ellipses
                     ellipse.type="confidence",
                     legend.title = "Treatment",invisible="quali"
)
# labels for countries have to be added manually. I can do this later for this figure if required
p1aa<-p1a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.5, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -1.3, y = -0.2, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 2.4, y = 0.3, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -1.7, y = -1.2, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 2.5, y = -0.5, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 1.7, y = -1.8, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.0, y = 2.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.2, y = 1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.5, y = 1.4, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p1,p1aa)

#####################
# PCA 2 = F9 (columns 1-3 are qualitative)
pca2<-PCA(data2,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca2$eig
pca2$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA2_ind<-cbind(data2[,c(1:3)], pca2$ind$coord)
p2<-plot(pca2,choix="var",axes=c(1,2),cex=1.4)
q2<-plot(pca2,choix="var",axes=c(1,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca2, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca2, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
p2a<-fviz(pca2, title = "Female PCA - 9 traits (F9) PC1 vs PC2",
          element = "ind",
          habillage =  data2$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p2aa<-p2a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.6, y = 0.1, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = 1.7, y = 0.4, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.9, y = 1.0, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 1.0, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = -1.5, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.1, y = 1.8, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 2.8, y = -1.2, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.6, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.2, label = "Spain", cex = 5, colour = "#B82E00") 
# labels can be added later

p2b<-fviz(pca2, title = "Female PCA - 9 traits (F9) PC1 vs PC3",
          element = "ind", axes = c(1, 3),
          habillage =  data2$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p2bb<-p2b+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.7, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.4, y = 0.9, label = "Ger.", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.8, y = -0.6, label = "Rus.", cex = 5, colour = "#095888") +
  annotate("text", x = -1.1, y = 1.3, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = 2.3, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.9, y = -1.3, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.0, y = -0.2, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -2.0, y = -0.8, label = "Spain", cex = 5, colour = "#B82E00") 
# labels can be added later

plot_grid(p2,p2aa, q2, p2bb)

#####################
# PCA 3 = Fmax  (columns 1-3 are qualitative)
pca3<-PCA(data3,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca3$eig
pca3$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA3_ind<-cbind(data3[,c(1:3)], pca3$ind$coord)
p3<-plot(pca3,choix="var",axes=c(1,2),cex=1.4)
q3<-plot(pca3,choix="var",axes=c(1,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca3, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca3, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
p3a<-fviz(pca3, title = "Female PCA - 12 traits (Fmax) PC1 vs PC2",
          element = "ind",
          habillage =  data3$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p3aa<-p3a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 1.8, y = 0.3, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = 0.2, label = "Ger.", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.4, y = 1.3, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.1, y = 1.2, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.4, y = -1.3, label = "Ukr.", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.2, y = 1.9, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -2.0, y = -0.2, label = "Spain", cex = 5, colour = "#B82E00") 
# labels can be added later

p3b<-fviz(pca3, title = "Female PCA - 12 traits (Fmax) PC1 vs PC3",
          element = "ind", axes = c(1, 3),
          habillage =  data3$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p3bb<-p3b+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.0, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.6, y =1.0, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 0.8, y = -1.2, label = "Rus.", cex = 5, colour = "#095888") +
  annotate("text", x = -1.8, y = 0.7, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = 2.7, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.1, y = -1.2, label = "Den.", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -2.2, y = -0.4, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.4, y = -2.0, label = "Spain", cex = 5, colour = "#B82E00") 
# labels can be added later

plot_grid(p3,p3aa, q3, p3bb)

#####################
# PCA = Fmaxplus (columns 1-3 are qualitative)
pca4<-PCA(data4,scale.unit=TRUE,graph=FALSE,quali.sup=c(1:3),axes=c(1,2))
pca4$eig
pca4$var$coord
# Here are the individual coordinates for PC2 output with the qualitative vars
PCA4_ind<-cbind(data4[,c(1:3)], pca4$ind$coord)
p4<-plot(pca4,choix="var",axes=c(1,2),cex=1.4)
q4<-plot(pca4,choix="var",axes=c(1,3),cex=1.4)
# simple ellipses with chosen colours
plotellipses(pca4, keepvar = c(1), axes = c(1, 2),label = "quali", level = 0.95, palette=palette)
plotellipses(pca4, keepvar = c(1), axes = c(1, 3),label = "quali", level = 0.95, palette=palette)
p4a<-fviz(pca4, title = "Female PCA - 13 traits (Fmax Plus) PC1 vs PC2",
          element = "ind",
          habillage =  data4$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p4aa<-p4a+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.2, y = 0.3, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.8, y = 0.2, label = "Ger.", cex = 5, colour = "#F6C200") +
  annotate("text", x = -0.2, y = -0.2, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.4, y = 1.2, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 1.4, y = 1.4, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = -0.2, y = 1.9, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.3, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = -1.8, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -1.9, y = -0.5, label = "Spain", cex = 5, colour = "#B82E00") 

p4b<-fviz(pca4, title = "Female PCA - 13 traits (Fmax Plus) PC1 vs PC3",
          element = "ind", axes = c(1,3),
          habillage =  data4$Country,
          geom = c("point","text"),
          label= "quali",
          pointsize = 2,
          pointshape = 18,
          palette = palette,
          addEllipses = TRUE, # Concentration ellipses
          ellipse.type="confidence",
          legend.title = "Treatment",invisible="quali"
)
p4bb<-p4b+theme(text = element_text(size = 16), 
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14))+
  annotate("text", x = 2.0, y = 0.7, label = "Austria", cex = 5, colour = "#E38800") + 
  annotate("text", x = -0.6, y =1.0, label = "Germany", cex = 5, colour = "#F6C200") +
  annotate("text", x = 1.7, y = -1.0, label = "Russia", cex = 5, colour = "#095888") +
  annotate("text", x = -2.4, y = -0.4, label = "Finland", cex = 5, colour = "#A00E00") +
  annotate("text", x = 0.9, y = 2.7, label = "Ukraine", cex = 5, colour = "#0086A8") +
  annotate("text", x = 0.6, y = -1.5, label = "Denmark", cex = 5, colour = "#D04E00") +
  annotate("text", x = 3.1, y = -0.5, label = "Portugal", cex = 5, colour = "#7BA354") +
  annotate("text", x = -1.8, y = 0.7, label = "Turkey", cex = 5, colour = "#132B69") +
  annotate("text", x = -0.9, y = -1.8, label = "Spain", cex = 5, colour = "#B82E00") 

plot_grid(p4,p4aa, q4, p4bb)

#######
save(pca1, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/M9_drosEU.RData")
save(pca2, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/F9_drosEU.RData")
save(pca3, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/Fmax_drosEU.RData")
save(pca4, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/FmaxP_drosEU.RData")

write.csv(PCA1_ind, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/M9_drosEU_PCcoords.csv", row.names = F)
write.csv(PCA2_ind, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/F9_drosEU_PCcoords.csv", row.names = F)
write.csv(PCA3_ind, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/Fmax_drosEU_PCcoords.csv", row.names = F)
write.csv(PCA4_ind, file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/FmaxP_drosEU_PCcoords.csv", row.names = F)

save.image(file = "/Users/ewanharney/Dropbox/Barcelona_IBE/DrosEU/PCA_drosEU_220927.RData")
