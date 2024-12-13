### Initial stats on ws23 microbiome

### setup

# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("indicspecies")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(indicspecies)

parseDistanceDF = function(phylip_file) {
    
    # Read the first line of the phylip file to find out how many sequences/samples it contains
    temp_connection = file(phylip_file, 'r')
    len = readLines(temp_connection, n=1)
    len = as.numeric(len)
    len = len +1
    close(temp_connection)
    
    
    phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
    colnames(phylip_data) <- row.names(phylip_data)
    return(phylip_data)
}

### read in data

otu <- read.table(file="../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.0.03.subsample.shared",
                  header=T, stringsAsFactors = FALSE, row.names = 2)
str(otu)
alpha <- read.table(file="../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.groups.ave-std.summary",
                    header=T, stringsAsFactors = FALSE)
str(alpha)

#### challenge Filter alpha to keep all columns but only rows that have average (ave) values

alpha <- filter(alpha,method == "ave")
# alpha <- filter(alpha, method != "std")

t
taxa <- read.table(textConnection(gsub("\\(.+?\\);", "\t", 
    readLines("../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy"))), 
    col.names=c("OTU", "Size", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), skip=1)
taxa <- taxa[taxa$OTU %in% names(otu),]

# read in environmental data, sample name linker
expdata <- read.table(file= "../may18ws.env.txt", header=T, stringsAsFactors = TRUE)
samples <- read.table(file = "../may18ws.sample.txt", header=T)

expdata <- left_join(expdata, samples, by="Sample")

# join expdata with alpha, selects just the samples with mothur data
alpha.expdata <- left_join(alpha, expdata, by="group")
expdata.alpha <- left_join(expdata, alpha, by="group")

## reorder factors so they display water, sed, soil
alpha.expdata$Type <- factor(alpha.expdata$Type, levels = c("water", "sed", "soil"))

## read in beta diversity
jc <- parseDistanceDF("../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.jest.0.03.lt.ave.dist")
bc <- parseDistanceDF("../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.braycurtis.0.03.lt.ave.dist")
tyc <- parseDistanceDF("../ws23.trim.contigs.good.unique.good.filter.precluster.denovo.vsearch.pick.opti_mcc.thetayc.0.03.lt.ave.dist")

# search help
?read.table


# manupilate data
#base R 
#square bracket object[row,column]
otu2 <- otu[,-1]
otu2 <- otu2[,-1]

#tidyverse
otu <- select(otu, -label, -numOtus)

maxab <- apply(otu, 2,max)
str(maxab)
# make vector all names of OTUs that are less than 1% of any sample (subsample to 7500)
n1 <- names(which(maxab<75))
str(n1)

#make otu table with just abundant OTUs
otu.ab <- otu[,-which(names(otu) %in% n1)]


type.col <- c(
    "water" = "#66c2a5",
    "sed" = "#8da0cb",
    "soil" = "#fc8d62"
)

### Alpha diversity!

ggplot(data=alpha.expdata, mapping=aes(x=Type, y=sobs))+
           geom_boxplot()
 
## Challenge make boxplots of samples grouped by Type for coverage, shannon, and invsimpson

ggplot(data=alpha.expdata, mapping=aes(x=Type, y=coverage))+
    geom_boxplot()

ggplot(data=alpha.expdata, mapping=aes(x=Type, y=shannon))+
    geom_boxplot()

ggplot(data=alpha.expdata, mapping=aes(x=Type, y=invsimpson))+
    geom_boxplot(varwidth = T, outlier.shape = NA)+
    geom_jitter(width=0.1, color="blue", size=4)+
    labs(y="Inverse Simpson")+
    ggtitle("Bacterial diversity by Sample Type")+
    theme_bw()

## run Anova to test for sig differences

res.aov <- aov(invsimpson~Type, data=alpha.expdata)
summary(res.aov)
TukeyHSD(res.aov)    
    
    
### Beta diversity

## jc=Jaccard presence absence

jc.nms <- metaMDS(as.dist(jc), k=2, try=50, trymax=500, wascores=F)

# use ordiplot for 3D ordinations
ordiplot(jc.nms, choices=c(1,2), type="points")

## plot ordination with ggplot

##### Challenge Improve Ordination plots, pick better colors, 
  #make points bigger/different shape, add plot title, fix legend title

jc.points <- data.frame(jc.nms$points)

jc.plot <- ggplot(jc.points, aes(x=MDS1, y=MDS2, label=rownames(jc)))

x <- max(jc.points$MDS1)/1.5
y <- max(jc.points$MDS2)


jc.plot +
    geom_point(aes(fill = factor(alpha.expdata$Type)), size = 3, shape = 23)+
    # scale_shape_manual(values = 23)+ # doesn't work because we're not mapping shape
    scale_fill_manual(values= type.col)+
    ggtitle("Jaccard")+
    labs(fill = "Sample Type")+
    annotate("text", x, y, label=paste("stress = ", round(jc.nms$stress, digits = 3)))
    
## bc = Bray Curtis

bc.nms <- metaMDS(as.dist(bc), k=2, try=50, trymax=500, wascores=F)
#
## plot ordination with ggplot

bc.points <- data.frame(bc.nms$points)

bc.plot <- ggplot(bc.points, aes(x=MDS1, y=MDS2, label=rownames(bc)))

x <- max(bc.points$MDS1)/2
y <- min(bc.points$MDS2)


bc.plot +
    geom_point(aes(fill = factor(alpha.expdata$Type)), size = 3, shape = 23)+
    # geom_text()+ # labels each sample with row.name
    # scale_shape_manual(values = 23)+ # doesn't work because we're not mapping shape
    scale_fill_manual(values= type.col)+
    ggtitle("Beta Diversity Bray Curtis")+
    labs(fill = "Sample Type")+
    annotate("text", x, y, label=paste("stress = ", round(bc.nms$stress, digits = 3)))

## tyc = Theta YC

tyc.nms <- metaMDS(as.dist(tyc), k=2, try=50, trymax=500, wascores=F)

## plot ordination with ggplot

tyc.points <- data.frame(tyc.nms$points)

tyc.plot <- ggplot(tyc.points, aes(x=MDS1, y=MDS2, label=rownames(tyc)))

x <- max(tyc.points$MDS1)/1.5
y <- max(tyc.points$MDS2)


tyc.plot +
    geom_point(aes(color = factor(alpha.expdata$Type)))+
    annotate("text", x, y, label=paste("stress = ", round(tyc.nms$stress, digits = 3)))


### hypothesis testing Permanova 

permanova <- adonis2(as.dist(jc)~alpha.expdata$Type, perm=99, rm.na=TRUE)
permanova
betadisp <- betadisper(as.dist(jc), alpha.expdata$Type)
anova(betadisp) ## p=0.56, no sig diff between dispersion of samples

permanova <- adonis2(as.dist(bc)~alpha.expdata$Type, perm=99, rm.na=TRUE)
permanova
betadisp <- betadisper(as.dist(bc), alpha.expdata$Type)
anova(betadisp) # no sig diff

permanova <- adonis2(as.dist(tyc)~alpha.expdata$Type, perm=99, rm.na=TRUE)
permanova
betadisp <- betadisper(as.dist(tyc), alpha.expdata$Type)
anova(betadisp) #no sig diff


### Indicator Species Analysis

indic <- multipatt(otu, alpha.expdata$Type, control = how(nperm=99))
summary(indic)

write.csv(file="indicator.species.csv", indic$sign %>%
              rownames_to_column(var = "OTU") %>%
              mutate(p.fdr = round(p.adjust(p.value, "fdr"),3)) %>%
              right_join(taxa, by = "OTU") %>%
              # filter(p.fdr <0.5) %>%
              filter(p.value <0.05)%>%
              arrange(index))


#repeat indicator analysis on just abundant otu

indic <- multipatt(otu.ab, alpha.expdata$Type, control = how(nperm=999))
summary(indic)

write.csv(file="indicator.species.999.csv", indic$sign %>%
              rownames_to_column(var = "OTU") %>%
              mutate(p.fdr = round(p.adjust(p.value, "fdr"),3)) %>%
              right_join(taxa, by = "OTU") %>%
              # filter(p.fdr <0.5) %>%
              filter(p.value <0.05)%>%
              arrange(index))


## probably shouldn't use indic for microbiome. should try ALDEx2 and ANCOM-II 