#### Species Richness by polygon feature from moedeled species rasters
#### Developed by Felipe de Medeiros Costa

# load libraries
#
# # install.packages('SDMTools')
# library(SDMTools)
# # # install.packages("spatialEco")
# library(spatialEco)
# # # install.packages('raster')
# library(raster)
# # # install.packages('sp')
# library(sp)
# # # install.packages("GISTools")
# library(GISTools)
# # # install.packages('rgdal')
# library(rgdal)
# # # install.packages("rgeos")
# library(rgeos)
# # # install.packages('plyr')
# library(plyr)
# # # install.packages('psych')
# library(psych)
# # install.packages('sf')
# library(sf)


# find the folder that contains all the individual models; 
# creates folder for output saves;
#
getwd()
setwd("C:/Users/felipe/Documents/TCC_Felipe/script_R_felipe")
dir.create("specieShpfiles", showWarnings=T)

#### load rasters (.asc) in folder (wd) into a variable; 
#### load shape mask into a variable
#
myFiles <- list.files("C:/Users/felipe/Documents/Thresholded_Models", pattern=".asc$", full.names=T, recursive=F)
maskShp <- NULL
maskShp <- st_read(dsn="./mask_shp", layer="grid625kmSquare_BR_WGS84")
maskShp <- st_simplify(maskShp)
maskShp <- st_make_valid(maskShp)

# make a stack with all the individual files in it
#
myFiles.stack <- stack(myFiles)

#### run a loop to save individual point shapefiles for all species from the Thresholded rasters intersected by a polygon mask;
#### create a frequency table(1) of how many time a specie appears in each mask polygon's feature, creating unique pairs of species x polygon;
#### create a frequency table(2), based on the previously created frequency table(1), of how many times a polygon ID appears, wich also means the Species Richness for each polygon;
# 

#### The next lines are to set the script to run only for high threat categories' species;
#### Coment those if you want to run to all Thresholded species;
#### Notice that the following loop has to be modified too, commenting the first of it's 'if' clauses;
#
speciesRedBook <- read.csv(file = './input_tables/especiesavaliadaslv2013.csv', header = T, sep = ';')
VUENCR.speciesRedBook <- c()
VUENCR.speciesRedBook <- speciesRedBook[which(speciesRedBook[,'categoria'] %in% c('VU','EN','CR')),'nome.científico']
write.csv(VUENCR.speciesRedBook, file = './speciesRB_VUENCR.csv', row.names = F)
####

loopCounter <- nlayers(myFiles.stack)
newRst <- NULL
newShp <- NULL
specieGrid.count <- NULL
uniquePairs_model <- NULL
modelSpeciesGrid.rich <- NULL
tmp.rich <- NULL
modeled.species <- c()

#### Make sure to have large amounts of RAM avaible. Sometimes it may take over 6GB of memory depending on raster size
#
memory.size(max = T)
dir.create("specieShpfiles/species_grid", showWarnings=T)

for(i in 1:loopCounter){
  print(myFiles.stack[[i]]@data@names)
  if(gsub(' ', '_', myFiles.stack[[i]]@data@names) %in% gsub(' ', '_', gsub(' ', '_', VUENCR.speciesRedBook), fixed = T)){
    newRst <- raster(nrows=myFiles.stack[[i]]@nrows, ncols=myFiles.stack[[i]]@ncols, xmn = myFiles.stack[[i]]@extent@xmin, xmx = myFiles.stack[[i]]@extent@xmax, ymn = myFiles.stack[[i]]@extent@ymin, ymx = myFiles.stack[[i]]@extent@ymax, crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"), vals=F)
    newRst <- newRst + myFiles.stack[[i]]
    newShp <- NULL
    newShp <- rasterToPoints(newRst, fun = function(x){x==1}, spatial = T)
    
    if(NROW(newShp@data)>0) {
      #### This next 'if' is to store the names of the species considered by the modelling process
      #
      if(!gsub(' ','_',myFiles.stack[[i]]@data@names) %in% gsub(' ', '_', modeled.species)){
        modeled.species <- append(modeled.species, gsub(' ', '_', myFiles.stack[[i]]@data@names), after = length(modeled.species))
      }
      #
      newShp <- spTransform(newShp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      print(paste("CRS transform sucess ",i))
      newShp$specie <- myFiles.stack[[i]]@data@names
      print(paste("Specie name added as attribute ",i))
      sf_newShp <- NULL
      sf_newShp <- st_as_sf(newShp)
      sf_newShp <- st_simplify(sf_newShp)
      sf_newShp <- st_make_valid(sf_newShp)
      intersect_newShp <- st_intersection(sf_newShp, maskShp)
      print(paste("Intersect OK ", i))
      st_write(obj = intersect_newShp, dsn = gsub(' ', '', paste('./specieShpfiles/species_grid/', gsub(' ', '_', myFiles.stack[[i]]@data@names), '.shp')), layer = myFiles.stack[[i]]@data@names, driver = 'ESRI Shapefile', delete_dsn = T)
      print(paste("saved intersect SHP ", i))
      
      #### On the next line, vars considered in c() must be set to the columns desired to be unique combinations
      #intersect_newShp <- as.SpatialPointsDataFrame.ppp(intersect_newShp)
      specieGrid.count <- count(intersect_newShp, c("GRID_ID","specie"))
      
      print(paste('Accounted for: ', i))
      uniquePairs_model <- rbind.data.frame(uniquePairs_model, specieGrid.count)
      tmp.rich <- as.vector(specieGrid.count$GRID_ID, mode="any")
      modelSpeciesGrid.rich <- append(modelSpeciesGrid.rich, tmp.rich, after=length(modelSpeciesGrid.rich))
      print(paste('Vector append sucessfully ',i))
    }  else {
      print(paste("No points for specie ",i," - ",myFiles.stack[[i]]@data@names))}
  }
}

#### Save files;
#### The table() will show how many times each ID appeas in the vector, wich by logical inference is the species richness by polygon
#
write.csv(modeled.species, file = './modeledSpecies_Grid.csv', row.names = F)
write.csv(uniquePairs_model, file='./modelUniquePairs_species_Grid.csv', row.names = T)
modelSpeciesGrid.rich <- table(modelSpeciesGrid.rich)
modelSpeciesGrid.rich <- as.data.frame(modelSpeciesGrid.rich)
colnames(modelSpeciesGrid.rich) = c('GRID_ID','richness')
write.csv(modelSpeciesGrid.rich, file='./modelRichness_Grid.csv', row.names = T)

#### Subtract unique pairs of polygonID-species, retrieved from previous loop, from_
#### table that has unique pairs in columns with the same name as generated by previous loop;
#### Notice there's a line to consider only high threat degree species. If you want all species comment that line;
#
model_Specie.Grid <- read.csv(file = './modelUniquePairs_species_Grid.csv')
model_Specie.Grid <- subset.data.frame(model_Specie.Grid, select = -freq)
model_Specie.Grid <- as.matrix(model_Specie.Grid)
ptsRB_Specie.Grid <- read.csv(file = "./input_tables/uniquePairs_ptsRB_VUENCR_grid625.csv", header = T, sep = ";" )
ptsRB_Specie.Grid <- subset.data.frame(ptsRB_Specie.Grid, select = -FREQUENCY)
ptsRB_Specie.Grid <- as.matrix(ptsRB_Specie.Grid)
species_to_remove <- c()
observed.species <- c()

for(k in 1:NROW(ptsRB_Specie.Grid)){
  #### The next 'if' is to register considered oberved species,_
  #### even though they will only be fully added to final shape later in the code;
  #
  if(!gsub(' ','_',ptsRB_Specie.Grid[k,'specie']) %in% observed.species){
    observed.species <- append(observed.species, gsub(' ','_',ptsRB_Specie.Grid[k,'specie']), after = length(observed.species))
  }
  for(l in 1:NROW(model_Specie.Grid)){
    if (ptsRB_Specie.Grid[k,'GRID_ID'] == model_Specie.Grid[l,'GRID_ID'] && 
    gsub(" ","_",ptsRB_Specie.Grid[k,'specie'], fixed = T) == gsub(" ","_",model_Specie.Grid[l,'specie'], fixed = T)) {
      print(paste('Match: ',k,' ',l))
      species_to_remove <- append(species_to_remove, ptsRB_Specie.Grid[k,'OBJECTID'])
      print(paste('Vector with species to remove sucessfulyy appended ',k,' ',l))
    }
  }
}

write.csv(observed.species, file = './observedSpecies_Grid.csv', row.names = F)

#### The next loop creates a vector, and writes a .csv, for all the species considered in the script (unique names)
#
allConsidered.species <- c()
allConsidered.species <- append(allConsidered.species, modeled.species, after = length(allConsidered.species))
for(k in 1:NROW(observed.species)){
    if(!observed.species[k] %in% modeled.species){
      allConsidered.species<- append(allConsidered.species, gsub(' ','_', observed.species[k]), after = length(allConsidered.species))
    }
}
write.csv(allConsidered.species, file = './allConsideredSpecies_Grid.csv', row.names = F)
####

#### Next lines are the subtraction of unique pairs that appears on both Modelled species and in situ samples;
#
uniquePairs_rmvSpecies <- NULL
uniquePairs_rmvSpecies <- ptsRB_Specie.Grid[!ptsRB_Specie.Grid[,'OBJECTID'] %in% c(species_to_remove),]

tmp.rich_rmvSpecies <- NULL
tmp.rich_rmvSpecies <- uniquePairs_rmvSpecies[,'GRID_ID']
count.richness_rmvSpecies <- NULL
count.richness_rmvSpecies <- table(tmp.rich_rmvSpecies)
count.richness_rmvSpecies <- as.data.frame(count.richness_rmvSpecies, ncol=2)
colnames(count.richness_rmvSpecies) <- c('GRID_ID','richness')
write.csv(count.richness_rmvSpecies, file='./ObservedRichness_minus_PotentialRichness_Grid.csv', row.names = T)

#### Summ up Modeled Species Richness table and Observed Species Richness table without the already subtracted matching unique pairs;
#### Note: 1 matrix to potential (modeled) richness - "potRich"; 1 matrix to observed (samples) richness - "obsRich";
#### Note: 1 matrix to store IDs from "potRich" in match cases only - "summed.richness";
#### Note: 1 vector to store match cases' rowIDs from observed richness table - "tmp.Summed.Observed.uniquePairs";
#### Note: 1 variable to sum up richness for each row in "potRich" - "richness"; It gets 0 every round, and is only changed if there's a match;
#### Note: The 'richness' column value in "summed.richness" receives the sum of all "obsRich" with matching IDs;
#
potRich <- read.csv(file='./modelRichness_Grid.csv', stringsAsFactors = F)
obsRich <- read.csv(file='./ObservedRichness_minus_PotentialRichness_Grid.csv', stringsAsFactors = F)

tmp.Summed.Observed.uniquePairs <- NULL
tmp.Summed.Observed.uniquePairs <- c()

summed.richness <- NULL
summed.richness <- matrix(nrow = 0,ncol = 3)
colnames(summed.richness) <- colnames(potRich)

for(k in 1:NROW(potRich)){
  richness <- 0
  # print(paste('First loop: ', k))
  for(l in 1:NROW(obsRich)){
    # print(paste('Second loop: ', l))  
    if(potRich[k,'GRID_ID'] == obsRich[l,'GRID_ID']){
          print(paste('Match: ',k,' ',l))
          richness <- richness + as.numeric(obsRich[l,'richness'])
          # tmp.sum <- c(potRich[k,])
          summed.richness <- rbind(summed.richness, potRich[k,])
          summed.richness[summed.richness[,'GRID_ID'] %in% potRich[k,'GRID_ID'],'richness'] <- richness
          tmp.Summed.Observed.uniquePairs <- append(tmp.Summed.Observed.uniquePairs, obsRich[l,'X'], after = length(tmp.Summed.Observed.uniquePairs))
          next
      }
  }
}

dfOrder(summed.richness, columns = 'GRID_ID')
# write.csv(summed.richness, file = './testSummedRich_grid.csv', row.names = F)

#### Now sum up potential richness ("potRich") and all matching cases' richness (summed.richness")
#
pot_obs_Rich <- NULL
pot_obs_Rich <- potRich
pot_obs_Rich[which(pot_obs_Rich[,'GRID_ID'] %in% summed.richness[,'GRID_ID']),'richness'] <- (as.numeric(pot_obs_Rich[which(pot_obs_Rich[,'GRID_ID'] %in% summed.richness[,'GRID_ID']),'richness']) 
                                                                                        + as.numeric(summed.richness[which(summed.richness[,'GRID_ID'] %in% pot_obs_Rich[,'GRID_ID']),'richness']))
#### And what about no match cases? Here they are!
#
notSummed.ObservedRichness <- matrix(ncol=3)
colnames(notSummed.ObservedRichness) <- colnames(obsRich)
notSummed.ObservedRichness <- obsRich[!obsRich[,'X'] %in% c(tmp.Summed.Observed.uniquePairs),]
dfOrder(notSummed.ObservedRichness, columns = 'GRID_ID')
# write.csv(notSummed.ObservedRichness, file = './test_notSummedObserved.csv', row.names = F)

#### Bind it all together!
#
allRichness <- NULL
allRichness <- rbind(allRichness, pot_obs_Rich)
allRichness <- rbind(allRichness,notSummed.ObservedRichness)

#### Here, the previously constructed maskShp SpatialPolygonsDataFrame, with all Grids, receives the calculated data
#
maskShp <- merge(maskShp, allRichness, by= 'GRID_ID')
maskShp <- subset.data.frame(maskShp, select = -X)
dfOrder(maskShp, columns = 'GRID_ID')
dir.create('./shp_totalRichness_by_Grid', showWarnings = T)

#### Here's your final shapefile!
#
st_write(obj = maskShp, dsn = './shp_totalRichness_by_Grid/shp_totalRichness_grid.shp', layer = 'totalRichness_by_Grid', driver='ESRI Shapefile', delete_dsn = T)

#### End of script
