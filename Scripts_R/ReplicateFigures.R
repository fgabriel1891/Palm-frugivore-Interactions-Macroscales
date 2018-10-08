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
#### Load analysis, dependencies and pre-define settings 
source("Scripts_R/ReplicateAnalisis.R")
library(plotrix)
library(circlize)
library(scales)



## Colors for Neotropics and Afrotropics 
col.afro <- "#f1a340"
col.neo  <- "#998ec3"


#--------------------------------------------------
#### FIGURE 1: WORKFLOW

## No script needed 


#--------------------------------------------------

#### FIGURE 2: SAMPLING COMPLETENESS 

## a) Accumulation curves representing sampling completeness (SC) of palm-frugivore interaction datasets in each biogeographic region (Neotropics, Afrotropics). 

setEPS() # Open EPS file
postscript("figs/SamplingCompletenessA.eps") # Set directory to save image
par(las = 1, mar =  c(6,6,2,2) ) # plot settings 
plot(spAcum1, 
     ylim = c(0,700),
     xlim = c(0,800),
     ylab = "",
     xlab = "Palm species",
     col = col.neo, 
     cex.axis = 1.5,
     cex.lab = 1.5
)
title(ylab="Frugivore species", line = 4, 
      cex.lab=1.5)
plotrix::ablineclip(h = est1$chao, 
                    x1 = -50,
                    x2 = 750,
                    v = 750, 
                    y1 = -50,
                    y2 = est1$chao,
                    lty = 1, 
                    col = col.neo,
                    lwd = 5)
plotrix::ablineclip(h = est1$chao-est1$chao.se, 
                    x1 = -50,
                    x2 = 750,
                    lty = 2, 
                    col = col.neo,
                    lwd = 2)
plotrix::ablineclip(h = est1$chao+est1$chao.se, 
                    x1 = -50,
                    x2 = 750,
                    v = 750, 
                    y1 = -50,
                    y2 = est1$chao+est1$chao.se,
                    lty = 2, 
                    col = col.neo,
                    lwd = 2)

legend(0,550,fill = c(col.neo, col.afro),
       legend = c("Neotropics SC = 47%",
                  "Afrotropics SC = 37%"),
       bty = "n", cex = 1.5)
plotrix::ablineclip(h = est2$chao, 
                    x1 = -50,
                    x2 = 220,
                    v = 220, 
                    y1 = -50,
                    y2 = est2$chao,
                    lty = 1, 
                    col = col.afro,
                    lwd = 5)
plotrix::ablineclip(h = est2$chao-est2$chao.se, 
                    x1 = -50,
                    x2 = 220,
                    y2 = est2$chao,
                    lty = 2, 
                    col = col.afro,
                    lwd = 2)
plotrix::ablineclip(h = est2$chao+est2$chao.se, 
                    x1 = -50,
                    x2 = 220,
                    v = 220, 
                    y1 = -50,
                    y2 = est2$chao+est2$chao.se,
                    lty = 2, 
                    col = col.afro,
                    lwd = 2)
plot(spAcum2, add =T,
     col = col.afro)

dev.off()

setEPS() # Open EPS file
postscript("figs/SamplingCompletenessB.eps") # Set directory to save image
## b) individual palm SC (expected/observed) in function of the number of studies and palm degree. 

par(las = 1,  mar = c(6,6,2,2))
plot(SC~NoStudies, 
     data = SCTableNeo[!SCTableNeo$NoStudies == 1,],
     col = col.neo, 
     pch = 16, 
     cex = log1p(SCTableNeo$NoInteractions[!SCTableNeo$NoStudies == 1]),
     xlab = "No. studies",
     ylab = "SC %(observed degree/expected degree) ",
     ylim = c(1,100),
     xlim = c(0,25),
     cex.lab = 1.5, 
     cex.axis = 1.5)
points(SC~NoStudies,
       data = SCTableAfr[!SCTableAfr$NoStudies == 1,],
       col = col.afro, 
       cex = log1p(SCTableAfr$NoInteractions[!SCTableNeo$NoStudies == 1]),
       pch = 16)
# plotrix::ablineclip(LmSC_SC$Neotropics, 
#                     x1 =2,
#                     x2 = max(SCTableNeo$NoStudies),
#                     lwd = 3,
#                     lty = 2,
#                     col = col.neo)
# plotrix::ablineclip(LmSC_SC$Afrotropics,
#                     x1 = 2,
#                     x2 = max(SCTableAfr$NoStudies),
#                     lwd = 3,
#                     lty = 2,
#                     col = col.afro)
legend("bottomright", 
       col = c(col.neo, col.afro),
       pch = c(16,16),
       c("Neotropics" ,"Afrotropics"), 
       bty = "n")
dev.off()

#--------------------------------------------------

#### FIGURE 3: GEOGRAPHICAL VARIATION OF PALM-FRUGIVORE INTERACTION RECORDS

### FIGURE 3 MADE WITH ARCGIS 

#--------------------------------------------------

#### FIGURE 4: REGIONAL META-NETWORKS

## Neotropics

## Chord Diagram for the Neotropics
# Build interactions matrix (at the family level to visualize better)
cex = 0.3 # set label size
int.fam.am <- table(neo$frugFamily,neo$PALM) # Create matrix
int.fam.am[int.fam.am >= 1] <- 1 # Make matrix binary
fam.matrix.am <- unclass(int.fam.am)
# Set the order of tracks 

order.frug.am <- unique(neo$frugFamily[order(neo$frugClass,
                                             neo$frugBodyMASS, 
                                             decreasing = F)]) # Set order of tracks frug
order.palms  <- unique(neo$PALM[order(neo$FruitLength_cm, 
                                               decreasing = T)]) # Set order of tracks palms

order.am <- c(as.character(order.palms), as.character(order.frug.am)) # Join order of tracks

# Determine colors for grids and links 

colNeo <- neo[,c("PALM","frugFamily","frugClass")]
colNeo$col <- ifelse(colNeo$frugClass =="BIRDS", "#efd14a",
                   ifelse(colNeo$frugClass =="MAMMAL","#7c4133",
                          ifelse(colNeo$frugClass =="FISH","#9cc0ec", 
                                 ifelse(colNeo$frugClass =="COLEOPTERA","#1a2042",
                                        ifelse(colNeo$frugClass =="REPTILES", "#c44031","#ff4ae2"))))) # Set colors by taxonomic class
# Grid colors
pallete.am <- c(colNeo$col,rep("#a3ae2c",length(colnames(fam.matrix.am))))
names(pallete.am) <- c(as.character(colNeo$frugFamily), colnames(fam.matrix.am))
# Link colors 
rn.am <- unique(colNeo)
link.col.am <- rn.am[,c(2,1,4)] 



circos.clear() # To ensure circos environment is clear
setEPS() ## Open EPS space
postscript("figs/ChordDiagramNeo.eps") # Define directory to save the figure
circos.par(start.degree=180, track.margin=c(0,0.01)) # initial settings

chordDiagram(fam.matrix.am,
             annotationTrack = "grid", 
             symmetric=FALSE,
             order = order.am, 
             grid.col = pallete.am, 
             col = link.col.am, 
             # preAllocateTracks = list(
             #   list(track.height = 0.3), list(track.height = 0.05)
             preAllocateTracks = 1)

# Add species labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] , 
              cex = cex, sector.name, 
              facing = "clockwise", 
              niceFacing = TRUE, 
              adj = c(0, 0))
  }, bg.border = NA)

# highlight.sector(unique(plot.america$PALM), track.index = 2, col = "#a3ae2c" )
# highlight.sector(unique(plot.america$Frug_Family[plot.america$CLASS == "MAMMAL"]), track.index = 2, col = "#7c4133" )
# highlight.sector(unique(plot.america$Frug_Family[plot.america$CLASS == "BIRDS"]), track.index = 2, col = "#efd14a" )
# highlight.sector(unique(plot.america$Frug_Family[plot.america$CLASS == "FISH"]), track.index = 2, col = "#9cc0ec" )
# highlight.sector(unique(plot.america$Frug_Family[plot.america$CLASS == "COLEOPTERA"]), track.index = 2, col = "#1a2042" )
# highlight.sector(unique(plot.america$Frug_Family[plot.america$CLASS == "REPTILES"]), track.index = 2, col = "#c44031" )
#legend("topleft", "a) Neotropics", bty="n")
legend("bottomleft", fill=c(unique(colNeo$col),"#a3ae2c"), 
       legend=c("MAMMALS", "FISH","BIRDS", "REPTILES","BEETLES","CRABS", "PALMS"),
       bty="n", 
       cex=0.3, border = c(unique(colNeo$col),"#a3ae2c"))
circos.clear()
dev.off()


# build the interaction matrix 
int.fam.af <- table(afr$frugFamily,afr$PALM)
int.fam.af[int.fam.af >= 1] <- 1
fam.matrix.af <- unclass(int.fam.af)

# Set the order of tracks

order.frug.af <- unique(afr$frugFamily[order( afr$frugClass, afr$frugBodyMASS, 
                                             decreasing = F)])
order.palms  <- unique(afr$PALM[order(afr$FruitLength_cm,
                                              decreasing = T)])
order.af <- c(as.character(order.palms), as.character(order.frug.af))


# Determine colors for grids and links 

col.af<-afr[,c("PALM","frugFamily","frugClass")]
col.af$col<-ifelse(col.af$frugClass=="BIRDS", "#efd14a",
                   ifelse(col.af$frugClass=="MAMMAL","#7c4133",
                          ifelse(col.af$frugClass=="FISH","#9cc0ec", 
                                 ifelse(col.af$frugClass=="COLEOPTERA","#1a2042",
                                        ifelse(col.af$frugClass=="REPTILES", "#006e10","#a8324f")))))

# Grid colors and set colors by taxonomical classes
pallete.af<-c(col.af$col,rep("#a3ae2c",length(colnames(fam.matrix.af))))
names(pallete.af)<-c(as.character(col.af$frugFamily), colnames(fam.matrix.af))


# Link colors 
rn<-unique(col.af)
link.col <- rn[,c(2,1,4)]
cex = 0.3 # set label size
circos.clear()
## Afrotropics
setEPS() ## Open EPS space
postscript("figs/ChordDiagramAfr.eps") # Define directory to save the figure
circos.par(start.degree=180, track.margin=c(0,0.01), gap.degree= 1)
chordDiagram(fam.matrix.af, 
             annotationTrack = "grid",
             symmetric=FALSE,
             order=order.af,
             grid.col = pallete.af, 
             col=link.col, 
             # preAllocateTracks = list(
             #   list(track.height = 0.3), list(track.height = 0.05))
             preAllocateTracks = 1)

# Add species labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] , cex = cex,sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0))
}, bg.border = NA)

# highlight.sector(unique(plot.africa$PALM), track.index = 2, col = "#a3ae2c" )
# highlight.sector(unique(plot.africa$Frug_Family[plot.africa$CLASS == "MAMMAL"]), track.index = 2, col = "#7c4133" )
# highlight.sector(unique(plot.africa$Frug_Family[plot.africa$CLASS == "BIRDS"]), track.index = 2, col = "#efd14a" )

legend("bottomleft", fill=c(unique(col.af$col),"#a3ae2c"), legend=c("MAMMALS","BIRDS", "PALMS"), bty="n", cex=cex, border = c(unique(col.af$col),"#a3ae2c"))
circos.clear()
dev.off()

#--------------------------------------------------

#### FIGURE 4: TRAIT MATCHING RELATIONSHIPS
setEPS() ## Open EPS space
postscript("figs/TraitMatching.eps")

## Neotropics

cex = 1
layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
par(las = 1)
# Plots all dataset 
plot(log1p(neo$frugBodyMASS)~log1p(neo$FruitLength_cm),
     xlab="Log[ fruit length (cm) ]", 
     ylab="Log[ median body mass (g) ]", 
     pch = 16,
     col = "grey",
     cex = cex,
     xlim = c(0,4), 
     ylim = c(0,15),
     main = "All species", cex.lab = cex)
points(matchTraitsNeo$frugBodyMASS~matchTraitsNeo$FruitLength_cm, 
       pch = 16, 
       col =  col.neo,
       xlim = c(0,4), 
       ylim = c(0,15),
       cex = 1.5)
# ablineclip(LmAllNeo, col = alpha(col.neo, 0.9), 
#            x1 = min(matchTraitsNeo$FruitLength_cm),
#            x2 = max(matchTraitsNeo$FruitLength_cm), 
#            lwd = 4, 
#            lty = 1)
#legend("bottomright", "p = 0.08", cex = cex)
#legend(-0.5,16, "a)", cex = cex, bty = "n")

# Plots bird dataset
bird.data <- frugiv.am[frugiv.am$PALM %in% matchTraitsNeo$PALM[matchTraitsNeo$frugClass == "BIRDS"],]


plot(log1p(frugBodyMASS)~log1p(FruitLength_cm),
     data = bird.data[bird.data$frugClass == "BIRDS",], 
     xlab="Log[ fruit length (cm) ]", 
     ylab="Log[ median body mass (g) ]", 
     pch = 16, 
     col = "grey",
     cex = cex, 
     xlim = c(0,4), 
     ylim = c(0,15),
     main = "Birds",cex.lab = cex)
points(frugBodyMASS~FruitLength_cm, 
       data = matchTraitsNeo[matchTraitsNeo$frugClass == "BIRDS",],
       pch = 16, 
       col =  col.neo,
       xlim = c(0,4), ylim = c(0,15),
       cex = 1.5)
#ablineclip(birds.neo.lm.palm, col = alpha(col.neo, 0.9), 
#x1 = 0.2,x2 = 2.6, lwd = 4, lty = 2)
#legend("bottomright", "p = 0.38", cex = cex)
#legend(-0.5,16, "b)", cex = cex, bty = "n")

# Plots mammal dataset 

mammal.data <- frugiv.am[frugiv.am$PALM %in% matchTraitsNeo$PALM[matchTraitsNeo$frugClass == "MAMMAL"],]

plot(log1p(frugBodyMASS) ~ log1p(FruitLength_cm),
     data = mammal.data[mammal.data$frugClass == "MAMMAL",], 
     xlab="Log[ fruit length (cm) ]", 
     ylab="Log[ median body mass (g) ]", 
     pch = 16, 
     col = "grey",
     cex = cex,
     xlim = c(0,4), 
     ylim = c(0,15),
     main = "Mammals",
     cex.lab = cex)
points(frugBodyMASS ~ FruitLength_cm , 
       data = matchTraitsNeo[matchTraitsNeo$frugClass == "MAMMAL",],
       pch = 16, 
       col =  col.neo,
       xlim = c(0,4), ylim = c(0,15),
       cex = 1.5)
#ablineclip(mammal.neo.lm.palm, col = alpha(col.neo, 0.9), 
# x1 = 0.2, x2 = 3, lwd = 4, lty = 4)
#legend("bottomright", "p = 0.89", cex = cex)
#legend(-0.5,16, "c)", cex = cex, bty = "n")

## Afrotropics 

# Plots all dataset 

plot(log1p(frugBodyMASS)~log1p(FruitLength_cm),
     data =  afr, 
     xlab="Log[ fruit length (cm) ]", 
     ylab="Log[ median body mass (g) ]", 
     pch = 16,
     col = "grey",
     cex = cex, xlim = c(0,4),
     ylim = c(0,15),
     main = "All species",
     cex.lab = cex)
points(matchTraitsAfr$frugBodyMASS~ matchTraitsAfr$FruitLength_cm, 
       pch = 16, 
       col =  col.afro,
       xlim = c(0,4),
       ylim = c(0,15),
       cex = 1.5)
ablineclip(LmAllAfr, 
           col = col.afro,
           x1 = 0.5,
           x2 = 2.6, 
           lwd = 4)
#legend("bottomright", "p = 0.02", cex = cex)
#legend(-0.5,16, "d)", cex = cex, bty = "n")

# Plots bird dataset

bird.data <- frugiv.af[frugiv.af$PALM %in% matchTraitsAfr$PALM[matchTraitsAfr$frugClass == "BIRDS"],]


plot(log1p(frugBodyMASS)~log1p(FruitLength_cm),
     data = bird.data[bird.data$frugClass == "BIRDS",], 
     xlab="Log[ fruit length (cm) ]",
     ylab="Log[ median body mass (g) ]", 
     pch = 16, 
     col = "grey", 
     cex = cex, 
     xlim = c(0,4), 
     ylim = c(0,15),
     main = "Birds", 
     cex.lab = cex)
points(frugBodyMASS ~ FruitLength_cm, 
       data = matchTraitsAfr[matchTraitsAfr$frugClass == "BIRDS",],
       pch = 16, 
       col =  col.afro,
       xlim = c(0,4), 
       ylim = c(0,15),
       cex = 1.5)
ablineclip(LmBIRDAfr,
           col = col.afro, 
           x1 = 0.5,
           x2 = 2.2,
           lwd = 4,
           lty = 1)
# legend("bottomright", "p = 0.03", cex = cex)
# legend(-0.5,16, "e)", cex = cex, bty = "n")

# Plots mammal dataset 

mammal.data <- frugiv.af[frugiv.af$PALM %in% matchTraitsAfr$PALM[matchTraitsAfr$frugClass == "MAMMAL"],]

plot(log1p(frugBodyMASS) ~ log1p(FruitLength_cm),
     data = mammal.data[mammal.data$frugClass == "MAMMAL",], 
     xlab="Log[ fruit length (cm) ]", 
     ylab="Log[ median body mass (g) ]", 
     pch = 16,
     col = "grey", 
     cex = cex,
     xlim = c(0,4), 
     ylim = c(0,15),
     main = "Mammals",
     cex.lab = cex)
points(frugBodyMASS ~ FruitLength_cm, 
       data = matchTraitsAfr[matchTraitsAfr$frugClass == "MAMMAL",],
       pch = 16, 
       col =  col.afro,
       xlim = c(0,4),
       ylim = c(0,15),
       cex = 1.5)
ablineclip(LmMAMMALAfr,
           col = col.afro, 
           x1 = 0.2,
           x2 = 2.6, 
           lwd = 4, 
           lty = 1)
# legend("bottomright", "p=0.04", cex = cex)
# legend(-0.5,16, "f)", cex = cex, bty = "n")
dev.off()




#### Appendix figures


# 1)  Palm degree (i.e. frugivore interactions) in function of No.Studies

setEPS()
postscript("figs/InteractionsVsStudies.eps")
par(las = 1, mar = c(6,6,2,2))
plot(SCTableNeo$NoInteractions~SCTableNeo$NoStudies, 
     col = col.neo, 
     pch = 16, 
     cex = 1.5,
     xlab = "No. studies",
     ylab = "Palm degree",
     xlim = c(0,25),
     cex.lab = 1.5, 
     cex.axis = 1.5)
plotrix::ablineclip(LmSCTable$Neotropics,
                    x1 = 1,
                    x2 = max(SCTableNeo$NoStudies),
                    col = col.neo,
                    lwd = 3
)
points(SCTableAfr$NoInteractions~SCTableAfr$NoStudies, 
       col = col.afro, 
       pch = 16,
       cex = 1.5)
plotrix::ablineclip(LmSCTable$Afrotropics,
                    x1 = 1,
                    x2 = max(SCTableAfr$NoStudies),
                    col = col.afro,
                    lwd = 3)
legend("bottomright", 
       col = c(col.neo, col.afro), 
       pch = c(16,16),
       bty = "n",
       c("Neotropics" ,"Afrotropics"))
dev.off()
