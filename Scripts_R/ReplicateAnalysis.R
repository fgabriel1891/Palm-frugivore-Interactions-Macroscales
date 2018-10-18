#### R Script to replicate the results presented at:
#### A synthesis of animal-mediated seed dispersal of palms reveals distinct biogeographic differences in species interactions 
#### The Journal of Biogeography. 2018.(in revision)
#### Gabriel Muñoz1,*, Kristian Trøjelsgaard2 & W. Daniel Kissling1

## 1Institute for Biodiversity and Ecosystem Dynamics (IBED), University of Amsterdam, P.O. Box 94248, 1090 GE Amsterdam, The Netherlands

## 2Faculty of Engineering and Science, Department of Chemistry and Bioscience, Section of Biology and Environmental Science, University of Aalborg. Denmark

## *Present address: NASUA, Andrade Marin E7-76 & Diego de Almagro, 170135 Quito, Ecuador.

## Comments or questions about this script: 
## Gabriel Muñoz mailto::fgabriel1891@gmail.com | GitHub: /fgabriel1891

#--------------------------------------------------
#### Load dependencies and pre-define settings: 

# Loading  libraries
library(dplyr)
library(vegan)
library(rgdal)
library(bipartite)
source("Scripts_R/CustomFunctions.R")

# Load palm distribution data at TDWG "Botanical Country" Level 3 

world.checklist <- read.csv("DATA/world.checklist.csv",
                            header=T, stringsAsFactors = F)

# Load shapefile of borders of botanical countries 

BotCountries <-  readOGR("DATA/BotanicalCountries/TDWG_level3_Coordinates.shp", 
                                layer = "TDWG_level3_Coordinates")
#--------------------------------------------------
#### Load dataset

Dataset = read.csv("DATA/PalmFrugDatasetOCT2018.csv", header = T)

#--------------------------------------------------
#### Data compilation

## No. records  of assembled species interactions dataset
dim(unique(Dataset))
## List variables included in dataset
names(Dataset)
## Number of unique palm-frugivore interactions 
length(unique(Dataset$InteractionID))
## Number of frugivore species 
length(unique(Dataset$FRUGIVORE))
## Number of palm species 
length(unique(Dataset$PALM))

#### Descriptive statistics of the dataset for the Neotropics
## Number of interactions recorded for the Neotropics 
dim(Dataset[Dataset$biogeographicRegion == "Neotropics",])
## Number of unique interactions for the Neotropics 
length(unique(Dataset$InteractionID[Dataset$biogeographicRegion == "Neotropics"]))
## Number of frugivore species for the Neotropics 
length(unique(Dataset$FRUGIVORE[Dataset$biogeographicRegion == "Neotropics"]))
## Number of palm species for the Neotropics 
length(unique(Dataset$PALM[Dataset$biogeographicRegion == "Neotropics"]))

#### Descriptive statistics of the dataset for the Afrotropics
## Number of interactions recorded for the Afrotropics 
dim(Dataset[Dataset$biogeographicRegion == "Afrotropics",])
## Number of unique interactions for the Afrotropics 
length(unique(Dataset$InteractionID[Dataset$biogeographicRegion == "Afrotropics"]))
## Number of frugivore species for the Afrotropics 
length(unique(Dataset$FRUGIVORE[Dataset$biogeographicRegion == "Afrotropics"]))
## Number of palm species for the Afrotropics 
length(unique(Dataset$PALM[Dataset$biogeographicRegion == "Afrotropics"]))

#--------------------------------------------------
#### Quantification of knowledge gaps

### Estimating sampling completeness at the dataset level 

# Subset dataset for the Neotropics and Afrotropics
neo = Dataset[Dataset$biogeographicRegion == "Neotropics",]
neo = droplevels(neo) # Drop extra "PALM" levels 
afr = Dataset[Dataset$biogeographicRegion == "Afrotropics",]
afr = droplevels(afr) # Drop extra "PALM" levels 
# Create binary matrix of frugivore occurrences per palm species 
t1 = table(neo$PALM, neo$FRUGIVORE) # Neotropics
t1[t1 >= 1] = 1 # Turn into presence/absence matrix
t2 = table(afr$PALM, afr$FRUGIVORE) # Afrotropics
t2[t2 >= 1] = 1  # Turn into presence/absence matrix

# Calculate accumulation curves
spAcum1 = vegan::specaccum(t1, "random",permutations = 100) # Neotropics
spAcum2 = vegan::specaccum(t2, "random", permutations  = 100) # Afrotropics 

# Calculate estimated asymptotes (richness of frugivores in function of observed palm species)
est1 = vegan::specpool(t1)
est2 = vegan::specpool(t2)

# Calculate sampling completness of the interaction dataset

SC1 = est1$Species/est1$chao # Neotropics 
SC2 = est2$Species/est2$chao # Afrotropics 
## Estimating individual sampling completness for all palm species 

# Erasing those records from Zona and Henderson
DatasetFC = Dataset[!Dataset$referenceKey == "FAIRCHILD",]
# Create individual datasets for each region 
neo1 = DatasetFC[DatasetFC$biogeographicRegion == "Neotropics",] 
neo1 = droplevels(neo1) # Drop extra "PALM" levels 
afr1 = DatasetFC[DatasetFC$biogeographicRegion == "Afrotropics",] 
afr1 = droplevels(afr1) # Drop extra "PALM" levels

# apply custom function to estimate Sampling Completeness per individual palm species 

SCTableNeo = data.frame(t(sapply(as.character(unique(neo1$PALM)), 
                                  function(x) makeSC(neo1, x)))) # Neotropics
SCTableNeo[SCTableNeo$NoStudies > 1 & SCTableNeo$NoInteractions > 1,] # Show only those with more than 1 study and 1 interaction

SCTableAfr = data.frame(t(sapply(as.character(unique(afr1$PALM)), 
                                 function(x) makeSC(afr1, x)))) # Afrotropics
SCTableAfr[SCTableAfr$NoStudies > 1 & SCTableAfr$NoInteractions > 1,] # Show only those with more than 1 study and 1 interaction



## Calculate Pearson correlations

PearsonSCTable = list("Neotropics" = cor(SCTableNeo),
                      "Afrotropics" = cor(SCTableAfr))
## Calculate linear regressions

# Regression as Interactions in function of No. studies 
LmSCTable = list("Neotropics" = lm(NoInteractions~NoStudies,SCTableNeo), 
                 "Afrotropics" = lm(NoInteractions~NoStudies, SCTableAfr))

# Regression as Interactions in function of  Sampling completeness 
# (erasing those of SC = 100 because is most likely that they have only one frugivore and one study )
LmSC_SC = list("Neotropics" = lm(SC~NoStudies,SCTableNeo[!SCTableNeo$NoStudies == 1,]), 
               "Afrotropics" = lm(SC~NoStudies, SCTableAfr[!SCTableAfr$NoStudies == 1,]))

# Order dataset based on the most representative species (i.e. most number of interactions)
SCTableNeo = SCTableNeo[order(SCTableNeo$NoInteractions, decreasing = T),]
SCTableAfr = SCTableAfr[order(SCTableAfr$NoInteractions, decreasing = T),]

### Knowledge gaps expresed in relation of Botanical Countries.

# Subset dataset so it contains only records for which BotCountry data is available.
BotCountry <- Dataset[!is.na(Dataset$botanicalCountry),]
# Create binary matrix of unique interactions in function of bot. countries
BotCountryMat <- table(BotCountry$InteractionID,BotCountry$botanicalCountry)
BotCountryMat[BotCountryMat>=1] <- 1
# Calculate number of interactions per botanical country
BotCountryDat <- data.frame("num"=colSums(BotCountryMat))
# create an id as column 
BotCountryDat$id <- rownames(BotCountryDat)

#Generate table to estimate ratio of species recorded/existent 
bot <- data.frame("bot"= unique(BotCountry$botanicalCountry))
# Calculate number of unique palms per bot. country (recorded)
for(i in 1:length(bot$bot)){
  bot$number[i]<-length(
    unique(BotCountry$PALM[BotCountry$botanicalCountry==bot$bot[i]]))
}

# Generate table to get total species (existent) per botanical country

wlist <- data.frame("bot"=unique(world.checklist$BOT_Country))

for(i in 1:length(wlist$bot)){
  wlist$number[i] <- length(unique(
    world.checklist$SpecName[world.checklist$BOT_Country==wlist$bot[i]]))
}

# Match total with recorded species

bot$total <- wlist$number[match(bot$bot,wlist$bot)]
bot$int <- BotCountryDat$num[match(bot$bot, BotCountryDat$id)]

# Calculate ratio of species per botanical country 

bot$ratio <- bot$number/bot$total
bot$ratio[which(bot$ratio > 1)] <- 1 # Set maximum values no more than 1
BotCountriesData <- subset(BotCountries, BotCountries$LEVEL_3_CO %in% bot$bot[complete.cases(bot$bot)])
BotCountriesTable <- data.frame("id" = droplevels(BotCountriesData$LEVEL_3_CO), 
                                "lat"= BotCountriesData$Lat, 
                                "lon" = BotCountriesData$Long,
                                "name" = droplevels(BotCountriesData$LEVEL_NAME))
BotCountriesTable$ratio <- bot$ratio[match( BotCountriesTable$id, bot$bot)]
BotCountriesTable$total <- bot$int[match( BotCountriesTable$id, bot$bot)]

  
ma1 = table(BotCountry$botanicalCountry, BotCountry$PALM)
ma1[ma1>=1]=1


BotCountriesTable$PalmPerC = rowSums(ma1)
BotCountriesTable$IntPPalm = BotCountriesTable$total/BotCountriesTable$PalmPerC

SCALL = rbind(SCTableAfr[SCTableAfr$NoStudies > 1 & SCTableAfr$NoInteractions > 1,] , 
              SCTableNeo[-which(rownames(SCTableNeo) == "Elaeis guineensis"),] # Because its found in both datasets
              [SCTableNeo$NoStudies > 1 & SCTableNeo$NoInteractions > 1,]) # We elimiate those with only one study and one interaction because it inflates the estimations of SC.

world.checklist$SC = SCALL$SC[match( world.checklist$SpecName,rownames(SCALL))]

world.checklist$SC[which(is.na(world.checklist$SC))] = 0 # For those species without information we assume sampling completness is 0 
aggBR = aggregate( world.checklist$SC, by = list("BC" = world.checklist$BOT_Country) , FUN = mean)# Average SC per botanical country (includei)

BotCountriesTable$SCm = aggBR$x[match(BotCountriesTable$id,aggBR$BC)]  # Match averages 
BotCountriesTable$rich = bot$total[match(BotCountriesTable$id,bot$bot)] # Add total palm richness 
BotCountriesTable$KG = BotCountriesTable$rich/ BotCountriesTable$SCm# Calculate KG as the inverse of the ratio SC / richness
BotCountriesTable$KG = BotCountriesTable$KG/0.0967258420 # relative to Colombia
BotCountriesTable$KG2 = BotCountriesTable$IntPPalm/BotCountriesTable$rich


#  write.csv(x = BotCountriesTable, "BotCount.csv")

#--------------------------------------------------
#### Comparison of networks

## Calculate and compare networks metrics at the Metanetwork level (Time consuming step!)
# Uncomment line below to replicate 
# MetaNetIndex = MetaNetworkAnalisis(net1, net2, 1)

## Calculate and comparte networks by subsampling the Neotropics networks to the same number of palm species as the Afrotropics 

# Subset to a dataset containing only Neotropics and Afrotropics pair-wise interactions

net1 = Dataset[Dataset$biogeographicRegion == "Neotropics",][c("PALM","FRUGIVORE")]
net1 = droplevels(net1)
net2 = Dataset[Dataset$biogeographicRegion == "Afrotropics",][c("PALM","FRUGIVORE")]  
net2 = droplevels(net1)

# Perform replicates and calculations 
# Time consuming step!, uncomment line below to replicate
#sensitivity = sensAnNet(net1, net2, NReps = 100, SpOrIn = "Sp", colToSample = 1, repNull = 999)

#--------------------------------------------------
#### Functional trait matching 

## Extracting trait information per interaction partners 

frugiv.am <- neo[c("PALM",
                   "FRUGIVORE", 
                   "frugClass",
                   "frugBodyMASS", 
                   "FruitLength_cm")][neo$frugClass == "MAMMAL" | neo$frugClass =="BIRDS",]


frugiv.af <- afr[c("PALM",
                   "FRUGIVORE", 
                   "frugClass", 
                   "frugBodyMASS",
                   "FruitLength_cm")][afr$frugClass == "MAMMAL" | afr$frugClass =="BIRDS",]

## Median matching traits between palms and frugivores (by medians)

# Aggregate frugivore mass trait data by medians
matchTraitsNeo <- aggregate( frugBodyMASS ~ PALM, 
                             frugiv.am, 
                             FUN = median )
matchTraitsAfr <- aggregate( frugBodyMASS ~ PALM, 
                             frugiv.af,
                             FUN = median )
# Aggregate fruit lenght trait data by medians
matchTraitsNeo$FruitLength_cm <- frugiv.am$FruitLength_cm[match(matchTraitsNeo$PALM, frugiv.am$PALM)] 
matchTraitsAfr$FruitLength_cm <- frugiv.af$FruitLength_cm[match(matchTraitsAfr$PALM, frugiv.af$PALM)] 
## Testing normality assumptions and transform if necessary 
shapiro.test(matchTraitsNeo$frugBodyMASS)
shapiro.test(matchTraitsAfr$frugBodyMASS)
## Log transform data to fulfill normality assumptions
matchTraitsNeo[,c("frugBodyMASS", "FruitLength_cm")] <- log1p(matchTraitsNeo[,c("frugBodyMASS", "FruitLength_cm")]) # Neotropics
matchTraitsAfr[,c("frugBodyMASS", "FruitLength_cm")] <- log1p(matchTraitsAfr[,c("frugBodyMASS", "FruitLength_cm")]) # Afrotropics
# Finding dispersal syndrome by calculating the mode 
mode.dispersal <- aggregate(frugClass ~ PALM, Dataset, FUN = Mode)
# Agreggate dispersal syndrome info 
matchTraitsNeo$frugClass = droplevels(mode.dispersal$frugClass[match(matchTraitsNeo$PALM,mode.dispersal$PALM )])
matchTraitsAfr$frugClass = droplevels(mode.dispersal$frugClass[match(matchTraitsAfr$PALM,mode.dispersal$PALM )])

## Calculate linear regressions between matching traits 

## All species
LmAllNeo <- lm(matchTraitsNeo$frugBodyMASS~matchTraitsNeo$FruitLength_cm) # Neotropics 
LmAllAfr <- lm(matchTraitsAfr$frugBodyMASS~matchTraitsAfr$FruitLength_cm) # Afrotropics 

## Bird dispersed species 
LmBIRDNeo <- lm(frugBodyMASS~FruitLength_cm, matchTraitsNeo[matchTraitsNeo$frugClass == "BIRDS",]) # Neotropics 
LmBIRDAfr <- lm(frugBodyMASS~FruitLength_cm, matchTraitsAfr[matchTraitsAfr$frugClass == "BIRDS",]) # Afrotropics 

## Mammal dispersed species
LmMAMMALNeo <- lm(frugBodyMASS~FruitLength_cm, matchTraitsNeo[matchTraitsNeo$frugClass == "MAMMAL",]) # Neotropics 
LmMAMMALAfr <- lm(frugBodyMASS~FruitLength_cm, matchTraitsAfr[matchTraitsAfr$frugClass == "MAMMAL",]) # Afrotropics 

#### END OF SCRIPT--------------------------------------------------




