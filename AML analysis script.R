#install package
install.packages("dplyr")
library(dplyr)


#//////////////////////////////read excel file single drug//////////////////////////////////////////////
library(readxl)
GILT_Day14_MV411 <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/GILT_Day14_MV411.xlsx")
TRAM_Day14_MOLM13 <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/TRAM_Day14_MOLM13.xlsx")
VEN_Day14_MOLM13 <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/VEN_Day14_MOLM13.xlsx") 

# //////////////// How to manage VEN_Day14_MOLM13 frame

#Rename column  
colnames(VEN_Day14_MOLM13)[2] <- 'num_sig'
colnames(VEN_Day14_MOLM13)[13] <- 'median_lfc'
colnames(VEN_Day14_MOLM13)[9] <- 'median_pval'
colnames(VEN_Day14_MOLM13)[1] <- 'Gene'

#delete row
VEN_Day14_MOLM13 <- VEN_Day14_MOLM13[-1,]

#convert colum from charactor to numeric 
VEN_Day14_MOLM13$median_lfc <- as.numeric(VEN_Day14_MOLM13$median_lfc)
VEN_Day14_MOLM13$median_pval <- as.numeric(VEN_Day14_MOLM13$median_pval)

#select column 
GILT_Day14_MV411 <- select(GILT_Day14_MV411, Gene, num_sig, tiers, median_lfc, median_pval)
TRAM_Day14_MOLM13 <- select(TRAM_Day14_MOLM13 , Gene, num_sig, tiers, median_lfc, median_pval)
VEN_Day14_MOLM13 <- (select(VEN_Day14_MOLM13 , Gene, num_sig, median_lfc, median_pval))

#Fillter()
FilGILT_14D_MV411 <- GILT_Day14_MV411 %>% filter(median_pval< 0.05, abs(median_lfc) >= 1, tiers != "Unassigned" )
FilTRAM_14D_MOLM13 <- TRAM_Day14_MOLM13 %>% filter(median_pval< 0.05, abs(median_lfc) >= 1, tiers != "Unassigned" )
FilVEN_14D_MOLM13 <- VEN_Day14_MOLM13 %>% filter(median_pval< 0.05, abs(median_lfc) >= 1)

#####################################################################################################################################
#####################################################################################################################################
#/////////////////////////////////venndiagram//////////////////////////////////////////////////////
#////////////////////////////////singledrug///////////////////////////////////////////////

in12 <- intersect(FilGILT_14D_MV411$Gene, FilTRAM_14D_MOLM13$Gene)
in13 <- intersect(FilGILT_14D_MV411$Gene, FilVEN_14D_MOLM13$Gene)
in23 <- intersect(FilTRAM_14D_MOLM13$Gene, FilVEN_14D_MOLM13$Gene)
in123 <- intersect(in12, in23)

#show(in123) 7
#"ELAC1" "ERCC1"     "HIST1H2AH" "IGFL2"     "OAZ1"      "PGS1"      "PTPN2" 

#install package venndiagram
install.packages("VennDiagram") 

#sg drug
library("VennDiagram")  
grid.newpage()                                        
draw.triple.venn(area1 = 1018,                        
                 area2 = 1777,
                 area3 = 794,
                 n12 = 136,
                 n23 = 73,
                 n13 = 37,
                 n123 = 7,
                 fill = c("cornflowerblue", "darkorchid1", "chartreuse1"),
                 lty = "blank",
                 category = c("GILT", "TRAM", "VEN"))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#///////////////////////////////////resistance////////////////////////////////////////////////////////
#filter >or= 1 :resistance
resisGILT_14D_MV411 <- FilGILT_14D_MV411 %>% filter(median_lfc >= 1)
resisTRAM_14D_MOLM13 <- FilTRAM_14D_MOLM13 %>% filter(median_lfc >= 1)
resisVEN_14D_MOLM13 <- FilVEN_14D_MOLM13 %>% filter(median_lfc >= 1)

#/////////////////////////////////venndiagrame//////////////////////////////////////////////////////

r12 <- intersect(resisGILT_14D_MV411$Gene, resisTRAM_14D_MV411$Gene)
r13 <- intersect(resisGILT_14D_MV411$Gene, resisVEN_14D_MOLM13$Gene)
r23 <- intersect(resisTRAM_14D_MOLM13$Gene, resisVEN_14D_MOLM13$Gene)
r123 <- intersect(r12, r23)

library("VennDiagram")  
grid.newpage()                                        
draw.triple.venn(area1 = 228,                        
                 area2 = 971,
                 area3 = 743,
                 n12 = 136,
                 n23 = 40,
                 n13 = 13,
                 n123 = 0,
                 fill = c("cornflowerblue", "darkorchid1", "chartreuse1"),
                 lty = "blank",
                 category = c("GILT", "TRAM", "VEN"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#///////////////////////////////////////////////sensitive//////////////////////////////////////////////////

sensiGILT_14D_MV411 <- FilGILT_14D_MV411 %>% filter(median_lfc <= -1)
sensiTRAM_14D_MOLM13 <- FilTRAM_14D_MOLM13 %>% filter(median_lfc <= -1)
sensiVEN_14D_MOLM13 <- FilVEN_14D_MOLM13 %>% filter(median_lfc <= -1)

#/////////////////////////////////////////////venndiagrame/////////////////////////////////////////////////

s12 <- intersect(sensiGILT_14D_MV411$Gene, sensiTRAM_14D_MOLM13$Gene)
s13 <- intersect(sensiGILT_14D_MV411$Gene, sensiVEN_14D_MOLM13$Gene)
s23 <- intersect(sensiTRAM_14D_MOLM13$Gene, sensiVEN_14D_MOLM13$Gene)
s123 <- intersect(s12, s23)
#show s123 = ELAC1

library("VennDiagram")  
grid.newpage()                                        
draw.triple.venn(area1 = 730,                        
                 area2 = 806,
                 area3 = 51,
                 n12 = 54,
                 n23 = 3,
                 n13 = 4,
                 n123 = 1,
                 fill = c("cornflowerblue", "darkorchid1", "chartreuse1"),
                 lty = "blank",
                 category = c("GILT", "TRAM", "VEN"))


#####################################################################################################################################
#####################################################################################################################################
##//////////////////Stacked barchart////////////////////////
library(readxl)
genegroup <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/in123singledrug.xlsx")
genein123 <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/123singledrug.xlsx")
GILT <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/GILT.xlsx")
TRAM <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/TRAM.xlsx")
VEN <- read_excel("D:/AML Project/xlxs data/analysis/xlsx_data/VEN.xlsx")


#//stacked barchart
 #1 library
install.packages("ggeasy")
library(ggplot2)
library(RColorBrewer)
library(scales)

#GILTplot
GILTplot <- (ggplot(GILT, aes(x = gene, y = medianlfc, fill = gene )) + 
  geom_bar(stat="identity") +
  #geom_col(fill = "cornflowerblue") +
  geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
  theme_classic() + scale_fill_viridis_d() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  ggeasy::easy_rotate_labels(which = "x", angle = 90) +
  labs(x = "gene", y = "medianlfc") +
  labs(title = "GILTERITINIB") +
  labs(subtitle = "gene expression in Gilteritinib drug"))
  
#02
ggplot(GILT, aes(x = medianlfc, y = gene, fill = "gene" )) + 
               geom_bar(stat="identity") +
               #geom_col(fill = "cornflowerblue") +
               geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
               theme_classic() + scale_fill_viridis_d() +
               theme(legend.position = "none") + 
               theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
               ggeasy::easy_rotate_labels(which = "x", angle = 90) +
               labs(y = "gene", x = "medianlfc") +
               labs(title = "GILTERITINIB") +
               labs(subtitle = "gene expression in Gilteritinib drug")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#TRAMplot
TRAMplot <- (ggplot(TRAM, aes(x = gene, y = medianlfc, fill = gene )) + 
               geom_bar(stat="identity") +
               #geom_col(fill = "darkorchid1") +
               geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
               theme_classic() + scale_fill_viridis_d() +
               theme(legend.position = "none") + 
               theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
               ggeasy::easy_rotate_labels(which = "x", angle = 90) +
               labs(x = "gene", y = "medianlfc") +
               labs(title = "TRAMETINIB") +
               labs(subtitle = "gene expression in Trametinib drug"))

#02
ggplot(TRAM, aes(x = medianlfc, y = gene, fill = gene )) + 
  geom_bar(stat="identity") +
  #geom_col(fill = "darkorchid1") +
  geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
  theme_classic() + scale_fill_viridis_d() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  ggeasy::easy_rotate_labels(which = "x", angle = 90) +
  labs(y = "gene", x = "medianlfc") +
  labs(title = "TRAMETINIB") +
  labs(subtitle = "gene expression in Trametinib drug")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#VENplot
VENplot <- (ggplot(VEN, aes(x = gene, y = medianlfc, fill = gene )) + 
               geom_bar(stat="identity") +
               geom_col(fill = "aquamarine3") +
               geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
               theme_classic() + scale_fill_viridis_d() +
               theme(legend.position = "none") + 
               theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
               ggeasy::easy_rotate_labels(which = "x", angle = 90) +
               labs(x = "gene", y = "medianlfc") +
               labs(title = "VENETOCLAX") +
               labs(subtitle = "gene expression in Venetoclax drug"))

#02
ggplot(VEN, aes(x = medianlfc, y = gene, fill = gene )) + 
  geom_bar(stat="identity") +
  geom_col(fill = "aquamarine3") +
  geom_text(aes(label = medianlfc), size = 3, vjust = -0.3)+
  theme_classic() + scale_fill_viridis_d() +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  ggeasy::easy_rotate_labels(which = "x", angle = 90) +
  labs(y = "gene", x = "medianlfc") +
  labs(title = "VENETOCLAX") +
  labs(subtitle = "gene expression in Venetoclax drug")






