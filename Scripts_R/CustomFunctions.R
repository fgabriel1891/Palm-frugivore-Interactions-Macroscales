## Functions 

source("Scripts_R/Beckett2016LPA_wb_plus.R")

### Helper function to calculate the Sampling completeness of a given dataset with a column named 
# PALM for palms and referenceKey for the studies, x represent the name of the palm to calculate the SC 


makeSC = function(dataset, x){ 
  SROm = dataset[dataset$PALM == x,]
  SROm = droplevels(SROm)
  SROm1 = table(SROm$referenceKey,SROm$FRUGIVORE)
  Acum = vegan::specpool(SROm1)[c("Species", "chao", "chao.se")]
  RetObj = c("SC" = round(Acum$Species/Acum$chao, 3) * 100, 
             "NoStudies" = length(unique(SROm$referenceKey)),
             "NoInteractions" = Acum$Species,
             "chao" = Acum$chao, 
             "chao.se" = Acum$chao.se)
  return(RetObj)
}

#### Helper function to calculate individual species richness accumulation tables 

makeSAC2 = function(dataset, x){ 
  SROm = droplevels(dataset[dataset$PALM == x,])
  SROm1 = table(SROm$referenceKey,SROm$FRUGIVORE)
  Acum = vegan::specaccum(SROm1, "random", 100)
  return(Acum)
}


## Function to make individual frugivore sampling completeness curve accumulation plots
makeFrugPlot = function(dataset, x){ 
  col.afro <- "#f1a340"
  col.neo  <- "#998ec3"
  
  SROm = dataset[dataset$PALM == x,]
  SROm = droplevels(SROm)
  SROm1 = table(SROm$referenceKey,SROm$FRUGIVORE)
  Acum = vegan::specpool(SROm1)[c("Species", "chao", "chao.se")]
  
  col = ifelse(unique(SROm$biogeographicRegion) == "Neotropics",
               col.neo, col.afro)
  ylim = Acum$chao + Acum$chao.se
  par(las = 1, mar = c(5,5,1,0.5))
  plot(vegan::specaccum(SROm1, "random", 100), 
       xlab = "Number of articles",
       ylab = "Frugivores", 
       ylim = c(0, ylim + ylim/3),
       xlim = c(1, length(rownames(SROm1))),
       col = col, 
       lwd = 2,
       xaxt = "n",
       main = x,
       cex.lab = 1.5,
       cex.main = 1.3)
  axis(1, 
       seq(1, length(rownames(SROm1)), 1),
       seq(1, length(rownames(SROm1)), 1))
  legend("topright",
         legend = unique(SROm$biogeographicRegion),
         pch = 16, 
         col = col,
         bty = "n", 
         cex = 1)
  abline(h=Acum$chao,
         col = "red", 
         lwd = 3)
  abline(h=c(Acum$chao-Acum$chao.se,
             Acum$chao+Acum$chao.se), 
         lty = 2, 
         col = "red") 
  legend("topleft" , 
         paste("SC = ", round(Acum$Species/Acum$chao, 3) * 100, "%"),
         bty = "n", cex = 1)}




## Function to calculate network indices and compare the two networks at the "meta-network" level 


MetaNetworkAnalisis =  function(net1, net2, repNull){ 
  # Make binary matrices
  bin1 = as.matrix(table(net1))
  bin1[bin1 >= 1] = 1 # transform matrix from semi-quantitative to binary
  
  bin2 = as.matrix(table(net2))
  bin2[bin2 >= 1] = 1 # transform matrix from semi-quantitative to binary
  
  # Number of species 
  
  SpeciesA = c(length(row.names(bin1)), 
               length(row.names(bin2)))
  
  SpeciesB = c(length(colnames(bin1)), 
               length(colnames(bin2)))
  
  # Species A ratio Species B 
  
  ratioSp = SpeciesA/SpeciesB
  
  # connectance
  
  con_b1 = length(which(bin1 > 0 ))/
    (nrow(bin1)*ncol(bin1))
  con_b2 = length(which(bin2 > 0 ))/
    (nrow(bin2)*ncol(bin2))
  
  # mean number of interactions per species 
  
  In_b1 = specieslevel(bin1, index = "degree")
  In_b2 = specieslevel(bin2, index = "degree")
  
  # Nestedness 
  
  nesA1 = bipartite::nested(web = bin1,
                            method = "NODF2")
  nesA2 = bipartite::nested(bin2, 
                            method = "NODF2")
  # Modularity 
  
  mod1 =  bipartite::LPA_wb_plus(bin1)
  mod2 =  bipartite::LPA_wb_plus(bin2)
  
  NMod1 = GetModularInformation(bin1,
                                mod1)$number_of_modules
  NMod2 = GetModularInformation(bin2,
                                mod2)$number_of_modules
  
  # Null Model 
  
  null1 = NullModSen(bin1, repNull)
  null2 = NullModSen(bin2, repNull)
  
  # Calculate  z-score Modularity
  zMod.net1 = (mod1$modularity  - null1[[1]]["NulMod.mean"]) / null1[[1]]["NulMod.sd"]
  zMod.net2 = (mod2$modularity - null2[[1]]["NulMod.mean"]) / null2[[1]]["NulMod.sd"]
  
  # Calculate z-score Nestedness 
  
  zNes.net1 =  (nesA1 - null1[[1]]["NulNest.mean"]) / null1[[1]]["NulNest.sd"]
  zNes.net2 =  (nesA2 - null2[[1]]["NulNest.mean"]) / null2[[1]]["NulNest.sd"]
  
  
  # return object 
  
  ret1 = data.frame("Variable" = c("Number of species A",
                                   "Number of species B",
                                   "Species Ratio A/B",
                                   "connectance", 
                                   "Mean Interactions Higher",
                                   "Mean Interactions Lower",
                                   "Nestedness",
                                   "Nestedness z-score",
                                   "Modularity", 
                                   "Modularity z-score",
                                   "Number of modules"),
                    
                    "Net1" = c(SpeciesA[1], 
                               SpeciesB[1],
                               ratioSp[1],
                               con_b1,
                               mean(In_b1[[1]]$degree),
                               mean(In_b1[[2]]$degree),
                               nesA1,
                               zNes.net1,
                               mod1$modularity,
                               zMod.net1,
                               NMod1),
                    "Net2" = c(SpeciesA[2], 
                               SpeciesB[2],
                               ratioSp[2], 
                               con_b2, 
                               mean(In_b2[[1]]$degree),
                               mean(In_b2[[2]]$degree),
                               nesA2,
                               zNes.net2,
                               mod2$modularity,
                               zMod.net2,
                               NMod2)
  )
  
  return(ret1)
}



####
### Functions to perform sensitivity analysis 

## a) Subsetting so the Neotropical dataset has the same number of interactions as the afrotropical dataset. 

#' subNetbyInt: Function to subsample one of two datasets of interactions so both have the same number of interactions 
#'
#' @param net1 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @param net2 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @return two networks with the same number of interactions as the smallest network between net1 and net2 
#' @examples sensAnNet(net1,net2)

subNetbyInt = function(net1,net2){
  # Check order of networks, and find the smaller one 
  toSample = min(dim(unique(net1))[1],
                 dim(unique(net2))[1])
  # subsample both networks 
  if (dim(unique(net1))[1] == toSample){
    eqNetworks = list("net1" = unique(net1),
                      "net2" = droplevels(unique(net2)[sample(toSample),]))
    print("Network 2 was subsampled")
  } else {
    eqNetworks = list("net1" = droplevels(unique(net1)[sample(toSample),]),
                      "net2" = unique(net2))
    print("Network 1 was subsampled")
  }
  print(paste("both networks have now",
              toSample, "interactions"))
  return(eqNetworks)
}


## b)  Subsetting so the Neotropical dataset has the same number of palm species as the afrotropical dataset.

#' subNetbySp: Function to subsample one of two datasets of interactions so both have the same number of species 
#'
#' @param net1 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @param net2 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @param colToSample Column where having the species of interest to standardize, default = 1 
#' @return two networks with the same number of species as the smallest network between net1 and net2
#' @examples sensAnNet(net1,net2)

subNetbySp = function(net1, net2, colToSample = colToSample){ 
  
  # Check richness of networks of species, and find the smaller one 
  toSample = min(length(unique(net1[,colToSample])),
                 length(unique(net2[,colToSample])))
  # subsample both networks 
  if (length(unique(net1[,colToSample])) == toSample){
    eqNetworks = list("net1" = droplevels(net1),
                      "net2" = droplevels(net2[net2[,colToSample] 
                                               %in% sample(unique(net2[,colToSample])
                                                           ,toSample),]))
    print("Network 2 was subsampled")
  } else {
    eqNetworks = list("net1" = droplevels(net1[net1[,colToSample] 
                                               %in% sample(unique(net1[,colToSample])
                                                           ,toSample),]),
                      "net2" = droplevels(net2))
    print("Network 1 was subsampled")
  }
  print(paste("both networks have now",toSample, 
              "species of", names(net1)[colToSample]))
  return(eqNetworks)
}



## c) Helper function to calculate null models for modularity and nestedness 

#' NullModSen: Helper Function to calculate mean and sd of null models of modularity and nestedness for a given matrix/network
#'
#' @param matrix A network organized as a matrix in with species A as columns  Species B as rows.
#' @param rep  integer specifying the number of replicates for the null model matrices, i.e. shuffles of the web, passes to bipartite::nullmodel
#' @return a list with the mean and sd of calculated modularity and nestedness for the null model.
#' @examples NullModSen(matrix,rep)

NullModSen =  function(matrix, rep){
  
  # Generate null models "r1 null model": non-sequential algorithm for binary matrices that preserves the site (row) frequencies, but uses column marginal frequencies as probabilities of selecting species.
  
  nnm = simulate(vegan::nullmodel(matrix, "r1"), rep)
  
  # Calculate nestedness 
  null.nest = sapply(1:rep, function(x) bipartite::nested(nnm[,,x], method = "NODF2"))
  
  # Calculate modularity 
  null.mod = sapply(1:rep, function(x) bipartite::LPA_wb_plus(nnm[,,x])$modularity)
  
  null = list(c("NulMod" = c("mean" = mean(null.mod), "sd" = sd(null.mod)),
                "NulNest" = c("mean" = mean(null.nest), "sd" = sd(null.nest))))
  return(null)
  
}



## d) Function to calculate network indexes

#' indexCal: Function to calculate the network indexes presented in the palm-frugivore interactions in revision for the JBI 
#'
#' @param sampleNet Object results from the subsampled networks either from subNetbySp or subNetbyInt functions 
#' @return a dataframe with the network index calculations "mostly coming from bipartite" 

indexCal = function(sampleNet, repNull){
  # Make binary matrices
  bin1 = as.matrix(table(sampleNet$net1))
  bin1[bin1 >= 1] = 1 # transform matrix from semi-quantitative to binary
  
  bin2 = as.matrix(table(sampleNet$net2))
  bin2[bin2 >= 1] = 1 # transform matrix from semi-quantitative to binary
  
  # Number of species 
  
  SpeciesA = c(length(row.names(bin1)), 
               length(row.names(bin2)))
  
  SpeciesB = c(length(colnames(bin1)), 
               length(colnames(bin2)))
  
  # Species A ratio Species B 
  
  ratioSp = SpeciesA/SpeciesB
  
  # connectance
  
  con_b1 = length(which(bin1 > 0 ))/
    (nrow(bin1)*ncol(bin1))
  con_b2 = length(which(bin2 > 0 ))/
    (nrow(bin2)*ncol(bin2))
  
  # mean number of interactions per species 
  
  In_b1 = specieslevel(bin1, index = "degree")
  In_b2 = specieslevel(bin2, index = "degree")
  
  # Nestedness 
  
  nesA1 = bipartite::nested(web = bin1,
                            method = "NODF2")
  nesA2 = bipartite::nested(bin2, 
                            method = "NODF2")
  # Modularity 
  
  mod1 =  bipartite::LPA_wb_plus(bin1)
  mod2 =  bipartite::LPA_wb_plus(bin2)
  
  NMod1 = GetModularInformation(bin1,
                                mod1)$number_of_modules
  NMod2 = GetModularInformation(bin2,
                                mod2)$number_of_modules
  
  # Null Model 
  
  null1 = NullModSen(bin1, repNull)
  null2 = NullModSen(bin2, repNull)
  
  # Calculate  z-score Modularity
  zMod.net1 = (mod1$modularity  - null1[[1]]["NulMod.mean"]) / null1[[1]]["NulMod.sd"]
  zMod.net2 = (mod2$modularity - null2[[1]]["NulMod.mean"]) / null2[[1]]["NulMod.sd"]
  
  # Calculate z-score Nestedness 
  
  zNes.net1 =  (nesA1 - null1[[1]]["NulNest.mean"]) / null1[[1]]["NulNest.sd"]
  zNes.net2 =  (nesA2 - null2[[1]]["NulNest.mean"]) / null2[[1]]["NulNest.sd"]
  
  
  
  # return object 
  
  ret1 = data.frame("Variable" = c("Number of species A",
                                   "Number of species B",
                                   "Species Ratio A/B",
                                   "connectance", 
                                   "Mean Interactions Higher",
                                   "Mean Interactions Lower",
                                   "Nestedness",
                                   "Nestedness z-score",
                                   "Modularity", 
                                   "Modularity z-score",
                                   "Number of modules"),
                    
                    "Net1" = c(SpeciesA[1], 
                               SpeciesB[1],
                               ratioSp[1],
                               con_b1,
                               mean(In_b1[[1]]$degree),
                               mean(In_b1[[2]]$degree),
                               nesA1,
                               zNes.net1,
                               mod1$modularity,
                               zMod.net1,
                               NMod1),
                    "Net2" = c(SpeciesA[2], 
                               SpeciesB[2],
                               ratioSp[2], 
                               con_b2, 
                               mean(In_b2[[1]]$degree),
                               mean(In_b2[[2]]$degree),
                               nesA2,
                               zNes.net2,
                               mod2$modularity,
                               zMod.net2,
                               NMod2)
  )
  
  return(ret1)
}


### e) Wrapper function to calculate the sensitivity analysis proposed by the reviewers 


#' sensAnNet: Function to subsample one of two datasets of interactions so both have the same number of species 
#'
#' @param net1 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @param net2 A network organized as an edgelist in a dataframe with species A and Species B separated in two columns.
#' @param NReps Integer specifying the number of subsampling replicates for the sensitivity analisis 
#' @param SpOrIn Specify wheter subsampling networks by having the same number of species A|B or having the same number of interactions
#' @param colToSample Column where having the species of interest to standardize, only relevant if SpOrIn = "Sp". default = 1. 
#' @param repNull Integer specifying the number of randomiced matrices to create the null model. Default = 10
#' @return a list for which $IndexTable =  Table with the mean, max and sd of the network indexes calculated for the sensitivity analisis; $ZscoresNullModel =  the z-scores of modularity and nested values; $SubsampledNetworks = an list contaning the subsampled networks 
#' @examples sensAnNet(net1,net2,NReps = 3, SpOrIn = "Sp", colToSample = 1)


sensAnNet = function(net1, net2, NReps = 2, SpOrIn = "Sp", colToSample = 1, repNull = 10) { 
  
  # Generate replicates  
  
  # If subset by species is selected 
  
  if(SpOrIn == "Sp"){
    reps = replicate(NReps,
                     indexCal(
                       subNetbySp(net1, net2,colToSample), repNull
                     ))
  }
  
  # If subset by interactions is selected 
  
  if(SpOrIn == "Int"){
    reps = replicate(NReps,
                     indexCal(
                       subNetbyInt(net1, net2), repNull
                     ))
    
  }
  
  
  aggData = data.frame("Var" = reshape2::melt(sapply(1:dim(reps)[2],
                                                     function(x) reshape2::melt(reps[,x]$Variable)))$value,
                       "Net1" = reshape2::melt(sapply(1:dim(reps)[2],
                                                      function(x) reshape2::melt(reps[,x]$Net1)))$value,
                       "Net2" = reshape2::melt(sapply(1:dim(reps)[2],
                                                      function(x) reshape2::melt(reps[,x]$Net2)))$value
  )
  # Calculate mean and SD from replicates
  # Mean: 
  aggData1 =  aggregate(aggData$Net1, by = list(aggData$Var), FUN = mean) 
  names(aggData1) = c("Variable", "Net1_mean")
  aggData1$Net2_mean =  aggregate(aggData$Net2, by = list(aggData$Var), FUN = mean)$x
  #   # SD: 
  aggData1$Net1_sd =  aggregate(aggData$Net1, by = list(aggData$Var), FUN = sd)$x
  aggData1$Net2_sd =  aggregate(aggData$Net2, by = list(aggData$Var), FUN = sd)$x
  #   # Max 
  aggData1$Net1_max =  aggregate(aggData$Net1, by = list(aggData$Var), FUN = max)$x
  aggData1$Net2_max =  aggregate(aggData$Net2, by = list(aggData$Var), FUN = max)$x
  
  return(aggData1)
  
}



# Custom function to find the mode of a vector in a dataframe (To extract most common type of disperser)
# From: http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



