# Rapid-structural-network-changes-in-bird-communities

## R

The R scripts have been implemented on R version 4.2.2.

### Loading R packages

```{r setup, include=FALSE}

# Load packages

source("packages.R")

```

## Data

### Loading bird and geographic data

```{r setup, include=FALSE}

# Load data from the French Breeding Bird Survey (STOC-EPS)

dataprp <- readRDS("raw_data/dataprp.rds")

# Load geographical information on sites

square_centroid <- readRDS("raw_data/square_centroid.rds")

```

## Spatial species associations

Details on the calculation of species associations from spatial data are available in Rigal et al. (2022). The difference here is that spatial associations are not calculated for each year but for the entire period (2001-2017).

### Loading functions

```{r setup, include=FALSE}

# Load functions for spatial associations

source("function_spatial_associations.R")

```


### Calculating species spatial associations

Warning: it might take a while.

```{r}

# Parallelise computation as this might take a while

registerDoParallel(cores = c(detectCores()-1))

# Prepare data

dataperpoint <- dcast(dataprp, code_point + habit + zonebio ~ code_sp, fun.aggregate = sum,value.var="abond")

# Run the function

spatial_associations <- ddply(dataperpoint, .(habit, zonebio), .fun=main_fun_space, 
                           N.null=1000, occmin=0, 
                           .parallel = TRUE)
                           
# Keep all spatial associations                           

spatial_associations$ses_init <- spatial_associations$spatial_asso

# Set non-significant spatial associations to 0

spatial_associations$spatial_asso[spatial_associations$pval>0.05] <- NA

```

## Temporal species associations

### Loading functions

```{r setup, include=FALSE}

# Load functions for spatial associations

source("function_temporal_associations.R")

```

### Transforming data into time-series

```{r}
# Add 0 when site monitored but species not recorded

data_ts_0 <- as.data.frame(dataprp %>% group_by(code_square, zonebio, habit, code_sp, year) %>% summarize(abond=sum(abond), count=n()))

data_ts_0b <- ddply(data_ts_0, .(code_square, zonebio, habit), function(x){
  large_data <- dcast(droplevels(x), code_square + zonebio + habit + code_sp ~ year, fun.aggregate=sum, value.var="abond")
  long_data <- melt(large_data, id.vars=c("code_square","zonebio","habit","code_sp"))
  return(long_data)
})
names(data_ts_0b)[c(5,6)] <- c("year","abond")

data_ts_0b$year <- as.numeric(as.character(data_ts_0b$year))

data_ts_0b <- data_ts_0b[order(data_ts_0b[1],data_ts_0b[2],data_ts_0b[3],data_ts_0b[4],data_ts_0b[5]),]

# Calculate time-series length

data_ts_0b$seq_num <- cumsum(c(1, diff(data_ts_0b$year) != 1))

seq_length <- as.data.frame(data_ts_0b %>% group_by(code_square, zonebio, habit, code_sp, seq_num) %>% summarize(seq=n()))

data_ts_0b <- merge(data_ts_0b, seq_length, by=c("code_square", "zonebio", "habit", "code_sp", "seq_num"), all=T)

# Remove small time-series and species recorded less than 2 times (for using multispatial CCM, Clark et al. 2015)

data_ts_0c <- data_ts_0b[data_ts_0b$seq>5, c("code_square", "zonebio", "habit", "code_sp", "year", "abond","seq")]

to_remove <- data_ts_0c %>% group_by(code_square, zonebio, habit, code_sp) %>% summarize(nb_record = sum(sign(abond)))

data_ts_0c <- merge(data_ts_0c, to_remove, by=c("code_square", "zonebio", "habit", "code_sp"), all.x=T)

data_ts_0c <- data_ts_0c[data_ts_0c$nb_record>1,]

# Add coordinates of sites

data_ts_ready <- merge(data_ts_0c, square_centroid[,c("code_square", "long", "lat", "lon2", "lat2")], by="code_square", all.x=T)

# Detrend time-series if needed

data_ts_ready <- ddply(data_ts_ready, .(code_square, zonebio, habit, code_sp), .fun=detrend_data)
```

### Assess significance of temporal associations

```{r}

temporal_associations <- ddply(data_ts_ready, .(zonebio, habit), .fun=multisp_CCM, .progress="text")

# E dim

data_ts_ready <- readRDS("output/data_ts_ready.rds")
data_ts_ready$code_square <- sprintf("%06d", data_ts_ready$code_square)

Edim_temporal_associations <- ddply(data_ts_ready, .(zonebio, habit), .fun=Edim_multisp_CCM, .progress="text")
Edim_temporal_associations$species <- substr(Edim_temporal_associations$species,3,21)

temporal_associations_Edim <- merge(temporal_associations, Edim_temporal_associations, by=c("zonebio","habit","species"), all.x=T)

temporal_associations_signif <- droplevels(na.omit(temporal_associations_Edim[temporal_associations_Edim$rho<0.05,]))
temporal_associations_signif$spA <- substr(temporal_associations_signif$species,14,19)
temporal_associations_signif$spB <- substr(temporal_associations_signif$species,1,6)

```

### Merge spatial and temporal associations

```{r}

sp_tp_association <- merge(temporal_associations_signif, spatial_associations, by.x=c("habit","zonebio","spA","spB"), by.y=c("habit","zonebio","spi","spj"), all.x=T)

# also available 

sp_tp_association <- readRDS("output/sp_tp_association.rds")

```




## Species distances

### Functional distance

```{r}
# loading the trait matrix from Storchová & Hořák (2018) (https://doi.org/10.1111/geb.12709)

trait <- read.csv("raw_data/life_history_bird2.csv", header = TRUE)

species_name_data <- read.csv("raw_data/species_name_data.csv", header=T)
species_name_data$scientific_name2 <- sub(" ","_",species_name_data$scientific_name)

# removing duplicated traits and traits with multiple NA

table.na <- apply(trait,2,function(x){sum(is.na(x))})

row.names(trait) <- trait$Species

trait <- na.omit(trait[,c("LengthU_MEAN","WingU_MEAN","TailU_MEAN","BillU_MEAN","TarsusU_MEAN","WeightU_MEAN","Sexual.dimorphism","Clutch_MEAN","Broods.per.year","EggL_MEAN","EggW_MEAN","Egg_MASS","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.period","Incubation.sex","Association.outside.the.breeding.season","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")])


# transforming into a binomial format the relevant variables

trait[,c("Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")] <- apply(trait[,c("Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")], 2, function(x) replace(x,x>0,1))

trait[,c("Sexual.dimorphism","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.sex","Association.outside.the.breeding.season")] <- lapply(trait[,c("Sexual.dimorphism","Territoriality","Sedentary","Facultative.migrant","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_Y","Frugivore_Y","Granivore_Y","Arthropods_Y","Other.invertebrates_Y","Fish_Y","Other.vertebrates_Y","Carrion_Y","Omnivore_Y","Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B","Young","Association.during.nesting","Nest.type","Nest.building","Mating.system","Incubation.sex","Association.outside.the.breeding.season")], factor)

trait[,c("Incubation.period")] <- as.numeric(trait[,c("Incubation.period")])

# Construct a distance matrix using gower metric

mat.trait <- trait
mat.trait$code_sp <- species_name_data$code_sp[match(row.names(trait),species_name_data$scientific_name2)]
mat.trait <- subset(mat.trait, code_sp %in% levels(droplevels(dataprp$code_sp)))
mat.dist <- daisy(mat.trait[,1:58], metric="gower", type=list(asymm=c(22:58)))
attr(mat.dist, "Labels") <- mat.trait$code_sp

# selecting clustering algorithm

clust1 <- hclust(mat.dist, method = "ward.D")
clust2 <- hclust(mat.dist, method = "ward.D2")
clust3 <- hclust(mat.dist, method = "single")
clust4 <- hclust(mat.dist, method = "complete")
clust5 <- hclust(mat.dist, method = "average")
clust6 <- hclust(mat.dist, method = "mcquitty")
#clust7 <- hclust(mat.dist, method = "median") # this two methods are reported to crash
#clust8 <- hclust(mat.dist, method = "centroid")

list.clust <- list(clust1,clust2,clust3,clust4,clust5,clust6)

dissim <- data.frame(ward.D=rep(NA,(2^length(list.clust)-1)),ward.D2=rep(NA,(2^length(list.clust)-1)),single=rep(NA,(2^length(list.clust)-1)),
                   complete=rep(NA,(2^length(list.clust)-1)),average=rep(NA,(2^length(list.clust)-1)),mcquitty=rep(NA,(2^length(list.clust)-1)),
                   median=rep(NA,(2^length(list.clust)-1)),centroid=rep(NA,(2^length(list.clust)-1)),dissimilarity=rep(NA,(2^length(list.clust)-1)))

nb <- 0
for(i in 1:length(list.clust)){
  combinaison <- combn(length(list.clust),i)
  for(j in 1:ncol(combinaison)){
    nb <- nb+1
    print(nb)
    ensemble<-cl_ensemble(list=list.clust[combinaison[,j]])
    if(i==1){consensus <- list.clust[combinaison[,j]]
    }else{consensus <- cl_consensus(ensemble, method = "euclidean")}
    dissim[nb,c(combinaison[,j])] <- combinaison[,j]
    dissim[nb,9] <- cl_dissimilarity(mat.dist, consensus, method = "spectral")
  }
}

dissim[which.min(dissim$dissimilarity),]
dissim[order(dissim$dissimilarity),]
consensus <- cl_consensus(cl_ensemble(list=list.clust[c(5)]))
plot(consensus)

sigma <- var(mat.dist) + var(cophenetic(clust5))
threshold <- 2*sqrt(nrow(as.matrix(mat.dist))*sigma)

# as sigma < threshold we can continue and extract the cophenetic distances

coph.dist <- cophenetic(consensus)
coph.dist2 <- data.frame(t(combn(attr(coph.dist, "Labels"),2)), as.numeric(coph.dist))
names(coph.dist2) <- c("sp1","sp2","distance")

```

### Phylogenetic distance

```{r}

# Data from Ericson (https://birdtree.org/, info on https://doi.org/10.1016/j.cub.2014.03.011) sequenced species, 10 000 trees with 6670OTUS each

tree <- read.nexus("raw_data/output_10000.nex")

tree_cons <- ape::consensus(tree[c(1:2607,2609:8485,8487:10000)], p = 1, check.labels = TRUE)
tree_dist <- cophenetic.phylo(tree_cons)
tree_dist_long <- tabularize(tree_dist, name="distance")
tree_dist_long$distance2 <- tree_dist_long$distance*1e308*1e3 # sometime length issue as some species are very close

tree_dist_long$spi2 <- species_name_data$code_sp[match(tree_dist_long$spi, species_name_data$scientific_name2)]
ch <- levels(tree_dist_long$spi)[which(!(levels(tree_dist_long$spi) %in% levels(as.factor(species_name_data$scientific_name2))))]
for(i in 1:length(ch)){
  tree_dist_long$spi2[which(tree_dist_long$spi %in% ch[i])] <- species_name_data$code_sp[match(str_replace(ch[i], "_", " "), species_name_data$scientific_name)]
}

tree_dist_long$spj2 <- species_name_data$code_sp[match(tree_dist_long$spj, species_name_data$scientific_name2)]
ch <- levels(tree_dist_long$spj)[which(!(levels(tree_dist_long$spj) %in% levels(as.factor(species_name_data$scientific_name2))))]
for(i in 1:length(ch)){
  tree_dist_long$spj2[which(tree_dist_long$spj %in% ch[i])]<-species_name_data$code_sp[match(str_replace(ch[i], "_", " "), species_name_data$scientific_name)]
}

phylo_distance <- tree_dist_long
phylo_distance <- data.frame(phylo_distance, pair=paste(phylo_distance$spi2,phylo_distance$spj2,sep="_"))
phylo_distance <- data.frame(phylo_distance, pair2=paste(phylo_distance$spj2,phylo_distance$spi2,sep="_"))

```

### Niche overlap

```{r}

habitat_sp <- as.data.frame(dataprp %>% group_by(code_sp, habit) %>% summarize(ab=sum(abond)))
habitat_sp <- droplevels(subset(habitat_sp, !(habit=="NA_NA")))
habitat_sp_wide <- dcast(habitat_sp, code_sp~habit)
habitat_sp_wide[is.na(habitat_sp_wide)] <- 0
habitat_sp <- melt(habitat_sp_wide)
names(habitat_sp)[c(2:3)] <- c("habit","ab")
habitat_sp_count <- as.data.frame(habitat_sp %>% group_by(code_sp) %>% summarize(count=sum(ab)))
habitat_sp <- merge(habitat_sp,habitat_sp_count, by="code_sp",all=T)
habitat_sp$freq <- habitat_sp$ab/habitat_sp$count

habitat_sp_wide <- dcast(habitat_sp[,c("code_sp","habit","ab")], code_sp~habit)
row.names(habitat_sp_wide) <- habitat_sp_wide$code_sp 
habitat_sp_wide$code_sp <- NULL

# chi2 distance as we want to use frequency and comparing the distribution profile of each species in the different habitats, independently from the abundance.

habitat_dist_chi <- as.data.frame(habitat_sp_wide) 
habitat_dist_chi <- dudi.coa(habitat_dist_chi,scan=F)
habitat_dist_chi <- dist.dudi(habitat_dist_chi)
habitat_dist_chi <- tabularize(as.matrix(habitat_dist_chi))
habitat_dist_chi$assoc <- 1-habitat_dist_chi$assoc
names(habitat_dist_chi)[3] <- "dist"

```

### Specialisation distance

```{r}

# from Godet et al. (2015) (https://doi.org/10.1111/geb.12266)

sxi <- read.csv("raw_data/SXI_publi.csv")
sxi$code_sp <- as.factor(sxi$code_sp)

specialisation.dist <- function(x){
  sp <- levels(droplevels(x$code_sp))
  d_SGI <- data.frame(spj=sxi$code_sp[!(sxi$code_sp==sp)], dSGI=abs((sxi$SGIo[!(sxi$code_sp==sp)]-x$SGIo)))
  return(d_SGI)
}
specialisation_dist <- ddply(sxi, .(code_sp), .fun=specialisation.dist)

specialisation_dist$pair <- as.factor(paste0(specialisation_dist$code_sp,sep="_",specialisation_dist$spj))
names(specialisation_dist)[1] <- "spi"

```


### Combining species associations with species distances

```{r}

association_distance <- merge(sp_tp_association,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp1","sp2"), all.x=T)
association_distance <- merge(association_distance,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp2","sp1"), all.x=T)
association_distance$functional_distance <- rep(NA, nrow(association_distance))
association_distance$functional_distance[!is.na(association_distance$distance.x)] <- association_distance$distance.x[!is.na(association_distance$distance.x)]
association_distance$functional_distance[!is.na(association_distance$distance.y)] <- association_distance$distance.y[!is.na(association_distance$distance.y)]
association_distance$distance.x <- association_distance$distance.y <- NULL

association_distance <- merge(association_distance,phylo_distance[,c("distance","spi2","spj2")], by.x=c("spA","spB"), by.y=c("spi2","spj2"), all.x=T)
names(association_distance)[which(names(association_distance)=="distance")] <- "phylogenetic_distance"
association_distance$phylogenetic_distance[which(is.na(association_distance$phylogenetic_distance))] <- 0

association_distance <- merge(association_distance,specialisation_dist[,c("dSGI","spi","spj")], by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
names(association_distance)[which(names(association_distance)=="dSGI")] <- "specialisation_distance"

association_distance <- merge(association_distance,habitat_dist_chi, by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
association_distance$niche_overlap_distance <- -association_distance$dist
association_distance$dist <- NULL

# also available

association_distance <- readRDS("output/association_distance.rds")

```

### Combining species association by species pair

```{r}

# plot by biogeography and habitat

asso_to_plot <- association_distance
asso_to_plot$spatial_asso_scaled <- scale(asso_to_plot$spatial_asso,center=F)
asso_to_plot$temp_asso_scaled <- scale(asso_to_plot$temp_int2,center=F)

ggplot(droplevels(asso_to_plot), aes(x=spatial_asso_scaled, y=temp_asso_scaled))+
  geom_point(data=asso_to_plot,alpha=0.1,col="black")+
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D")+
  scale_x_continuous(name="Spatial associations (SES)", limits=c(-5,5))+
  scale_y_continuous(name="Temporal associations", limits=c(-5,5))+
  theme_classic() + theme(text=element_text(size=20),legend.position="none")

# plot by unique species pairs
  
asso_to_plot_grouped <- data.frame(association_distance %>% group_by(spA, spB) %>% summarize(temp_asso_grouped = mean(temp_int2), spatial_asso_grouped = mean(spatial_asso, na.rm=TRUE)))
asso_to_plot_grouped$spatial_asso_scaled <- scale(asso_to_plot_grouped$spatial_asso_grouped,center=F)
asso_to_plot_grouped$temp_asso_scaled <- scale(asso_to_plot_grouped$temp_asso_grouped,center=F)

ggplot(droplevels(asso_to_plot_grouped), aes(x=spatial_asso_scaled, y=temp_asso_scaled))+
  geom_point(data=asso_to_plot_grouped,alpha=0.1,col="black")+
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D")+
  scale_x_continuous(name="Spatial associations (SES)", limits=c(-5,5))+
  scale_y_continuous(name="Temporal associations", limits=c(-5,5))+
  theme_classic() + theme(text=element_text(size=20),legend.position="none")
  
# plot with non significant species pairs

pair_sp <- combn(levels(droplevels(dataprp$code_sp)),2)
asso_to_plot_na <- data.frame(spA=pair_sp[1,],spB=pair_sp[2,])
asso_to_plot_na <- merge(asso_to_plot_na, asso_to_plot_grouped[,c("spA","spB","temp_asso_grouped","spatial_asso_grouped")], by=c("spA","spB"), all.x=TRUE)
asso_to_plot_na <- merge(asso_to_plot_na, asso_to_plot_grouped[,c("spA","spB","temp_asso_grouped","spatial_asso_grouped")], by.x=c("spA","spB"),by.y=c("spB","spA"), all.x=TRUE)
asso_to_plot_na <- data.frame(asso_to_plot_na %>% replace(is.na(.), 0))

asso_to_plot_na$temp_asso_grouped <- (asso_to_plot_na$temp_asso_grouped.x + asso_to_plot_na$temp_asso_grouped.y)/2
asso_to_plot_na$spatial_asso_grouped <- (asso_to_plot_na$spatial_asso_grouped.x + asso_to_plot_na$spatial_asso_grouped.y)/2

asso_to_plot_na$spatial_asso_grouped.x <- asso_to_plot_na$spatial_asso_grouped.y <- asso_to_plot_na$temp_asso_grouped.x <- asso_to_plot_na$temp_asso_grouped.y <- NULL

asso_to_plot_na$spatial_asso_scaled <- scale(asso_to_plot_na$spatial_asso_grouped,center=F)
asso_to_plot_na$temp_asso_scaled <- scale(asso_to_plot_na$temp_asso_grouped,center=F)

ggplot(droplevels(asso_to_plot_na), aes(x=spatial_asso_scaled, y=temp_asso_scaled))+
  geom_point(data=asso_to_plot_na,alpha=0.1,col="black")+
  scale_color_viridis(alpha = 1, direction = 1,discrete = TRUE, option = "D")+
  scale_x_continuous(name="Spatial associations (SES)", limits=c(-5,5))+
  scale_y_continuous(name="Temporal associations", limits=c(-5,5))+
  theme_classic() + theme(text=element_text(size=20),legend.position="none")

```



## Testing the relationship between species associations and niche distances

### Preparing data for distances

```{r}

data_asso_dist <- droplevels(asso_to_plot_na[,c("spA","spB","temp_asso_grouped","spatial_asso_grouped")])

data_asso_dist <- merge(data_asso_dist,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp1","sp2"), all.x=T)
data_asso_dist <- merge(data_asso_dist,coph.dist2[,c("distance","sp1","sp2")], by.x=c("spA","spB"), by.y=c("sp2","sp1"), all.x=T)
data_asso_dist$functional_distance <- rep(NA, nrow(data_asso_dist))
data_asso_dist$functional_distance[!is.na(data_asso_dist$distance.x)] <- data_asso_dist$distance.x[!is.na(data_asso_dist$distance.x)]
data_asso_dist$functional_distance[!is.na(data_asso_dist$distance.y)] <- data_asso_dist$distance.y[!is.na(data_asso_dist$distance.y)]
data_asso_dist$distance.x <- data_asso_dist$distance.y <- NULL

data_asso_dist <- merge(data_asso_dist,phylo_distance[,c("distance2","spi2","spj2")], by.x=c("spA","spB"), by.y=c("spi2","spj2"), all.x=T)
names(data_asso_dist)[which(names(data_asso_dist)=="distance2")] <- "phylogenetic_distance"
data_asso_dist$phylogenetic_distance[which(is.na(data_asso_dist$phylogenetic_distance))] <- 0

data_asso_dist <- merge(data_asso_dist,specialisation_dist[,c("dSGI","spi","spj")], by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
names(data_asso_dist)[which(names(data_asso_dist)=="dSGI")] <- "specialisation_distance"

data_asso_dist <- merge(data_asso_dist,habitat_dist_chi, by.x=c("spA","spB"), by.y=c("spi","spj"), all.x=T)
data_asso_dist$niche_overlap_distance <- -data_asso_dist$dist
data_asso_dist$dist <- NULL

data_asso_dist_all <- data_asso_dist

# also available

data_asso_dist_all <- readRDS("output/data_asso_dist_all.rds")

```

### Preparing data for traits

```{r}

data_interaction_possible <- data_asso_dist_all

# add main habitat of species to control

max_habitat_sp <- data.frame(habitat_sp %>% group_by(code_sp) %>% slice(which.max(freq)))

data_interaction_possible <- merge(data_interaction_possible, max_habitat_sp[,c("code_sp","habit")], by.x="spA", by.y="code_sp", all.x=TRUE)
names(data_interaction_possible)[which(names(data_interaction_possible)=="habit")] <- "habitatA"
data_interaction_possible <- merge(data_interaction_possible, max_habitat_sp[,c("code_sp","habit")], by.x="spB", by.y="code_sp", all.x=TRUE)
names(data_interaction_possible)[which(names(data_interaction_possible)=="habit")] <- "habitatB"

# add nest type

data_interaction_possible <- merge(data_interaction_possible, mat.trait[,c("code_sp","Nest.type")], by.x="spA", by.y="code_sp", all.x=TRUE)
names(data_interaction_possible)[which(names(data_interaction_possible)=="Nest.type")] <- "NestA"
data_interaction_possible <- merge(data_interaction_possible, mat.trait[,c("code_sp","Nest.type")], by.x="spB", by.y="code_sp", all.x=TRUE)
names(data_interaction_possible)[which(names(data_interaction_possible)=="Nest.type")] <- "NestB"

# add similar diet

mat_trait_food <- mat.trait[,c("Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")]

mat_trait_food <- apply(mat_trait_food,2,as.numeric)

mat_dist_food <- daisy(mat_trait_food, metric="euclidean")
attr(mat_dist_food, "Labels") <- mat.trait$code_sp

mat_dist_food2 <- data.frame(t(combn(attr(mat_dist_food, "Labels"),2)), as.numeric(mat_dist_food))
names(mat_dist_food2) <- c("spA","spB","food_dist")
mat_dist_food2$food_eq <- ifelse(mat_dist_food2$food_dist==0,1,0)

data_interaction_possible <- merge(data_interaction_possible, mat_dist_food2, by.x=c("spA","spB"), by.y=c("spA","spB"), all.x=TRUE)
data_interaction_possible <- merge(data_interaction_possible, mat_dist_food2, by.x=c("spA","spB"), by.y=c("spB","spA"), all.x=TRUE)
data_interaction_possible$food_dist.x[which(is.na(data_interaction_possible$food_dist.x))] <- data_interaction_possible$food_dist.y[which(is.na(data_interaction_possible$food_dist.x))]
data_interaction_possible$food_eq.x[which(is.na(data_interaction_possible$food_eq.x))] <- data_interaction_possible$food_eq.y[which(is.na(data_interaction_possible$food_eq.x))]

data_interaction_possible$food_eq.y <- data_interaction_possible$food_dist.y <- NULL

names(data_interaction_possible)[which(names(data_interaction_possible)=="food_dist.x")] <- "food_dist"
names(data_interaction_possible)[which(names(data_interaction_possible)=="food_eq.x")] <- "same_food"

# finalise data preparation 

data_interaction_possible$obs_temp_asso <- ifelse(data_interaction_possible$temp_asso_grouped != 0, 1, 0)
data_interaction_possible$obs_asso <- ifelse(data_interaction_possible$temp_asso_grouped != 0 & data_interaction_possible$spatial_asso_grouped != 0, 1, 0)

data_interaction_possible$same_habitat <- ifelse(data_interaction_possible$habitatA == data_interaction_possible$habitatB,1,0)

data_interaction_possible$same_nest <- ifelse(data_interaction_possible$NestA == data_interaction_possible$NestB,1,0)

# also available

data_interaction_possible <- readRDS("output/data_interaction_possible.rds")

```

### Model

```{r}

# Test relationship between associations and ecological distances

mod1 <- glm(obs_temp_asso ~ same_habitat + same_nest + same_food, data=data_interaction_possible, family = binomial)
summary(mod1)
glm.diag.plots(mod1)

boxLabels <- c("Similar habitat","Similar nest type","Similar diet")
df <- data.frame(yAxis = length(boxLabels):1, 
                 boxOdds = exp(coef(mod1)[-1]), 
                 boxCILow = exp(coef(mod1)[-1]-1.96*summary(mod1)$coefficients[-1,2]), 
                 boxCIHigh = exp(coef(mod1)[-1]+1.96*summary(mod1)$coefficients[-1,2])
)

ggplot(df, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5, shape=1) +
  theme_modern()+
  ylab("") +
  xlab("Odds ratio") +
  annotate(geom = "text", y =0.5, x = 1, 
           label = paste0("McFadden R² = ", round(with(summary(mod1), 1 - deviance/null.deviance),2)), size = 3.5, hjust = -0.5)

sub_data <- data_interaction_possible[data_interaction_possible$phylogenetic_distance>0,]
sub_data[,c("functional_distance","phylogenetic_distance","niche_overlap_distance")] <- scale(sub_data[,c("functional_distance","phylogenetic_distance","niche_overlap_distance")])

# Test relationship between associations and specific traits

mod2 <- glm(obs_temp_asso ~ functional_distance + phylogenetic_distance +  niche_overlap_distance, data=sub_data, family = binomial)
summary(mod2)
glm.diag.plots(mod2)

boxLabels <- c("Functional distance","Phylogenetic distance","Niche overlap distance")
df <- data.frame(yAxis = length(boxLabels):1, 
                 boxOdds = exp(coef(mod2)[-1]), 
                 boxCILow = exp(coef(mod2)[-1]-1.96*summary(mod2)$coefficients[-1,2]), 
                 boxCIHigh = exp(coef(mod2)[-1]+1.96*summary(mod2)$coefficients[-1,2])
)

ggplot(df, aes(x = boxOdds, y = boxLabels)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5) +
  theme_modern()+
  ylab("") +
  xlab("Odds ratio") +
  annotate(geom = "text", y =0.5, x = 1, 
           label = paste0("McFadden R² = ", round(with(summary(mod2), 1 - deviance/null.deviance),2)), size = 3.5, hjust = 1.5)

```

### Process function and subset by habitat

```{r}

# Function to calculate metrics of species association networks

source("function_community_nb_link.R")

# Apply over all networks

network_structure_square <- ddply(droplevels(dataprp), .(code_square, year), .fun=community_nb_link, .parallel =F, .progress = "text")

network_structure_square <- merge(network_structure_square, square_centroid, by="code_square", all.x=T)

# Subset by main habitat type

habitat_square <- as.data.frame(dataprp %>% group_by(code_square, code_point, habit) %>% summarize(count1=n()))
habitat_square <- as.data.frame(habitat_square %>% group_by(code_square, habit) %>% summarize(count=n()))
max_habitat_square <- data.frame(habitat_square %>% group_by(code_square) %>% slice(which.max(count)))

max_habitat_square$habitat <- NA
max_habitat_square$habitat[which(max_habitat_square$habit %in% c("A_1","A_2","A_3"))] <- "Forest"
max_habitat_square$habitat[which(max_habitat_square$habit %in% c("B_1","B_2","B_3","C_1","C_2","C_4","F_1","G_1"))] <- "Natural openland"
max_habitat_square$habitat[which(max_habitat_square$habit %in% c("D_1","D_2","D_3","D_4","D_5"))] <- "Farmland"
max_habitat_square$habitat[which(max_habitat_square$habit %in% c("E_1","E_2","E_3"))] <- "Urban"

network_structure_square <- merge(network_structure_square, max_habitat_square[,c("code_square","habitat")], by="code_square", all.x=T)

```

## Model network dynamics

```{r}

mod_connectance_time <- gamm(connectance~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=network_structure_square, random = list(code_square=~1))
summary(mod_connectance_time$gam)

mod_asso_obs_time <- gamm(nb_association~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=network_structure_square, random = list(code_square=~1))
summary(mod_asso_obs_time$gam)

mod_modularity_time <- gamm(modularity~year+connectance+te(lon2,lat2,bs="tp",k=3), data=network_structure_square, random = list(code_square=~1))
summary(mod_modularity_time$gam)

mod_evenness_time <- gamm(evenness~year+connectance+te(lon2,lat2,bs="tp",k=3), data=network_structure_square, random = list(code_square=~1))
summary(mod_evenness_time$gam)

mod_species_time <- gamm(nb_sp~year+te(lon2,lat2,bs="tp",k=3), data=network_structure_square, random = list(code_square=~1))
summary(mod_species_time$gam)

```

# Differences in network structure 

```{r}

# Network structure on average in 2001

## nb_sp
nb_sp_2001 <- summary(mod_species_time$gam)$p.coef[2]*2001 + summary(mod_species_time$gam)$p.coef[1]
nb_sp_2001_up <- (summary(mod_species_time$gam)$p.coef[2]-1.96*summary(mod_species_time$gam)$se[2])*2001 + summary(mod_species_time$gam)$p.coef[1]+1.96*summary(mod_species_time$gam)$se[1]
nb_sp_2001_low <- (summary(mod_species_time$gam)$p.coef[2]+1.96*summary(mod_species_time$gam)$se[2])*2001 + summary(mod_species_time$gam)$p.coef[1]-1.96*summary(mod_species_time$gam)$se[1]

## nb_association
nb_asso_obs_2001 <- summary(mod_asso_obs_time$gam)$p.coef[2]*2001 + summary(mod_asso_obs_time$gam)$p.coef[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
nb_asso_obs_2001_up <- (summary(mod_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_asso_obs_time$gam)$se[2])*2001 + summary(mod_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_asso_obs_time$gam)$se[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
nb_asso_obs_2001_low <- (summary(mod_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_asso_obs_time$gam)$se[2])*2001 + summary(mod_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_asso_obs_time$gam)$se[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))

## connectance
connectance_2001 <- summary(mod_connectance_time$gam)$p.coef[2]*2001 + summary(mod_connectance_time$gam)$p.coef[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
connectance_2001_up <- (summary(mod_connectance_time$gam)$p.coef[2]-1.96*summary(mod_connectance_time$gam)$se[2])*2001 + summary(mod_connectance_time$gam)$p.coef[1]+1.96*summary(mod_connectance_time$gam)$se[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
connectance_2001_low <- (summary(mod_connectance_time$gam)$p.coef[2]+1.96*summary(mod_connectance_time$gam)$se[2])*2001 + summary(mod_connectance_time$gam)$p.coef[1]-1.96*summary(mod_connectance_time$gam)$se[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))

## modularity
modularity_2001 <- summary(mod_modularity_time$gam)$p.coef[2]*2001 + summary(mod_modularity_time$gam)$p.coef[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
modularity_2001_up <- (summary(mod_modularity_time$gam)$p.coef[2]-1.96*summary(mod_modularity_time$gam)$se[2])*2001 + summary(mod_modularity_time$gam)$p.coef[1]+1.96*summary(mod_modularity_time$gam)$se[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
modularity_2001_low <- (summary(mod_modularity_time$gam)$p.coef[2]+1.96*summary(mod_modularity_time$gam)$se[2])*2001 + summary(mod_modularity_time$gam)$p.coef[1]-1.96*summary(mod_modularity_time$gam)$se[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))

## evenness
evenness_2001 <- summary(mod_evenness_time$gam)$p.coef[2]*2001 + summary(mod_evenness_time$gam)$p.coef[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
evenness_2001_up <- (summary(mod_evenness_time$gam)$p.coef[2]-1.96*summary(mod_evenness_time$gam)$se[2])*2001 + summary(mod_evenness_time$gam)$p.coef[1]+1.96*summary(mod_evenness_time$gam)$se[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
evenness_2001_low <- (summary(mod_evenness_time$gam)$p.coef[2]+1.96*summary(mod_evenness_time$gam)$se[2])*2001 + summary(mod_evenness_time$gam)$p.coef[1]-1.96*summary(mod_evenness_time$gam)$se[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))

### example 
network_structure_square[which(network_structure_square$year <= 2002 & network_structure_square$nb_sp %in% c(31:34) & network_structure_square$connectance <= connectance_2001_up & network_structure_square$evenness >= evenness_2001_low ),]

# network structure on average in 2017

## nb_sp
nb_sp_2017 <- summary(mod_species_time$gam)$p.coef[2]*2017 + summary(mod_species_time$gam)$p.coef[1]
nb_sp_2017_low <- (summary(mod_species_time$gam)$p.coef[2]-1.96*summary(mod_species_time$gam)$se[2])*2017 + summary(mod_species_time$gam)$p.coef[1]+1.96*summary(mod_species_time$gam)$se[1]
nb_sp_2017_up <- (summary(mod_species_time$gam)$p.coef[2]+1.96*summary(mod_species_time$gam)$se[2])*2017 + summary(mod_species_time$gam)$p.coef[1]-1.96*summary(mod_species_time$gam)$se[1]

## nb_association
nb_asso_obs_2017 <- summary(mod_asso_obs_time$gam)$p.coef[2]*2017 + summary(mod_asso_obs_time$gam)$p.coef[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
nb_asso_obs_2017_low <- (summary(mod_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_asso_obs_time$gam)$se[2])*2017 + summary(mod_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_asso_obs_time$gam)$se[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
nb_asso_obs_2017_up <- (summary(mod_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_asso_obs_time$gam)$se[2])*2017 + summary(mod_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_asso_obs_time$gam)$se[1] + summary(mod_asso_obs_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))

## connectance
connectance_2017 <- summary(mod_connectance_time$gam)$p.coef[2]*2017 + summary(mod_connectance_time$gam)$p.coef[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
connectance_2017_low <- (summary(mod_connectance_time$gam)$p.coef[2]-1.96*summary(mod_connectance_time$gam)$se[2])*2017 + summary(mod_connectance_time$gam)$p.coef[1]+1.96*summary(mod_connectance_time$gam)$se[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))
connectance_2017_up <- (summary(mod_connectance_time$gam)$p.coef[2]+1.96*summary(mod_connectance_time$gam)$se[2])*2017 + summary(mod_connectance_time$gam)$p.coef[1]-1.96*summary(mod_connectance_time$gam)$se[1] + summary(mod_connectance_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$nb_sp))

## modularity
modularity_2017 <- summary(mod_modularity_time$gam)$p.coef[2]*2017 + summary(mod_modularity_time$gam)$p.coef[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
modularity_2017_low <- (summary(mod_modularity_time$gam)$p.coef[2]-1.96*summary(mod_modularity_time$gam)$se[2])*2017 + summary(mod_modularity_time$gam)$p.coef[1]+1.96*summary(mod_modularity_time$gam)$se[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
modularity_2017_up <- (summary(mod_modularity_time$gam)$p.coef[2]+1.96*summary(mod_modularity_time$gam)$se[2])*2017 + summary(mod_modularity_time$gam)$p.coef[1]-1.96*summary(mod_modularity_time$gam)$se[1] + summary(mod_modularity_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))

## evenness
evenness_2017 <- summary(mod_evenness_time$gam)$p.coef[2]*2017 + summary(mod_evenness_time$gam)$p.coef[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
evenness_2017_low <- (summary(mod_evenness_time$gam)$p.coef[2]-1.96*summary(mod_evenness_time$gam)$se[2])*2017 + summary(mod_evenness_time$gam)$p.coef[1]+1.96*summary(mod_evenness_time$gam)$se[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))
evenness_2017_up <- (summary(mod_evenness_time$gam)$p.coef[2]+1.96*summary(mod_evenness_time$gam)$se[2])*2017 + summary(mod_evenness_time$gam)$p.coef[1]-1.96*summary(mod_evenness_time$gam)$se[1] + summary(mod_evenness_time$gam)$p.coef[3]*mean(na.omit(network_structure_square$connectance))

### example 
network_structure_square[which(network_structure_square$year > 2015 & network_structure_square$nb_sp %in% c(27:30) & network_structure_square$connectance >= connectance_2017_low & network_structure_square$evenness <= evenness_2017_up ),]


```


## Differences in network structure by habitat

### Farmland

```{r}

# models for farmland

data_farmland <- network_structure_square[which(network_structure_square$habitat=="Farmland"),]

mod_farmland_connectance_time <- gamm(connectance~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_farmland, random = list(code_square=~1))
summary(mod_farmland_connectance_time$gam)

mod_farmland_asso_obs_time <- gamm(nb_association~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_farmland, random = list(code_square=~1))
summary(mod_farmland_asso_obs_time$gam)

mod_farmland_modularity_time <- gamm(modularity~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_farmland, random = list(code_square=~1))
summary(mod_farmland_modularity_time$gam)

mod_farmland_evenness_time <- gamm(evenness~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_farmland, random = list(code_square=~1))
summary(mod_farmland_evenness_time$gam)

mod_farmland_species_time <- gamm(nb_sp~year+te(lon2,lat2,bs="tp",k=3), data=data_farmland, random = list(code_square=~1))
summary(mod_farmland_species_time$gam)

# network structure on average in 2001 in farmland

## nb_sp
nb_sp_farmland_2001 <- summary(mod_farmland_species_time$gam)$p.coef[2]*2001 + summary(mod_farmland_species_time$gam)$p.coef[1]
nb_sp_farmland_2001_up <- (summary(mod_farmland_species_time$gam)$p.coef[2]-1.96*summary(mod_farmland_species_time$gam)$se[2])*2001 + summary(mod_farmland_species_time$gam)$p.coef[1]+1.96*summary(mod_farmland_species_time$gam)$se[1]
nb_sp_farmland_2001_low <- (summary(mod_farmland_species_time$gam)$p.coef[2]+1.96*summary(mod_farmland_species_time$gam)$se[2])*2001 + summary(mod_farmland_species_time$gam)$p.coef[1]-1.96*summary(mod_farmland_species_time$gam)$se[1]

## nb_association
nb_asso_obs_farmland_2001 <- summary(mod_farmland_asso_obs_time$gam)$p.coef[2]*2001 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
nb_asso_obs_farmland_2001_up <- (summary(mod_farmland_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_farmland_asso_obs_time$gam)$se[2])*2001 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_farmland_asso_obs_time$gam)$se[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
nb_asso_obs_farmland_2001_low <- (summary(mod_farmland_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_farmland_asso_obs_time$gam)$se[2])*2001 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_farmland_asso_obs_time$gam)$se[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))

## connectance
connectance_farmland_2001 <- summary(mod_farmland_connectance_time$gam)$p.coef[2]*2001 + summary(mod_farmland_connectance_time$gam)$p.coef[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
connectance_farmland_2001_up <- (summary(mod_farmland_connectance_time$gam)$p.coef[2]-1.96*summary(mod_farmland_connectance_time$gam)$se[2])*2001 + summary(mod_farmland_connectance_time$gam)$p.coef[1]+1.96*summary(mod_farmland_connectance_time$gam)$se[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
connectance_farmland_2001_low <- (summary(mod_farmland_connectance_time$gam)$p.coef[2]+1.96*summary(mod_farmland_connectance_time$gam)$se[2])*2001 + summary(mod_farmland_connectance_time$gam)$p.coef[1]-1.96*summary(mod_farmland_connectance_time$gam)$se[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))

## modularity
modularity_farmland_2001 <- summary(mod_farmland_modularity_time$gam)$p.coef[2]*2001 + summary(mod_farmland_modularity_time$gam)$p.coef[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
modularity_farmland_2001_up <- (summary(mod_farmland_modularity_time$gam)$p.coef[2]-1.96*summary(mod_farmland_modularity_time$gam)$se[2])*2001 + summary(mod_farmland_modularity_time$gam)$p.coef[1]+1.96*summary(mod_farmland_modularity_time$gam)$se[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
modularity_farmland_2001_low <- (summary(mod_farmland_modularity_time$gam)$p.coef[2]+1.96*summary(mod_farmland_modularity_time$gam)$se[2])*2001 + summary(mod_farmland_modularity_time$gam)$p.coef[1]-1.96*summary(mod_farmland_modularity_time$gam)$se[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))

## evenness
evenness_farmland_2001 <- summary(mod_farmland_evenness_time$gam)$p.coef[2]*2001 + summary(mod_farmland_evenness_time$gam)$p.coef[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
evenness_farmland_2001_up <- (summary(mod_farmland_evenness_time$gam)$p.coef[2]-1.96*summary(mod_farmland_evenness_time$gam)$se[2])*2001 + summary(mod_farmland_evenness_time$gam)$p.coef[1]+1.96*summary(mod_farmland_evenness_time$gam)$se[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
evenness_farmland_2001_low <- (summary(mod_farmland_evenness_time$gam)$p.coef[2]+1.96*summary(mod_farmland_evenness_time$gam)$se[2])*2001 + summary(mod_farmland_evenness_time$gam)$p.coef[1]-1.96*summary(mod_farmland_evenness_time$gam)$se[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))


# network structure on average in 2017 in farmland

## nb_sp
nb_sp_farmland_2017 <- summary(mod_farmland_species_time$gam)$p.coef[2]*2017 + summary(mod_farmland_species_time$gam)$p.coef[1]
nb_sp_farmland_2017_low <- (summary(mod_farmland_species_time$gam)$p.coef[2]-1.96*summary(mod_farmland_species_time$gam)$se[2])*2017 + summary(mod_farmland_species_time$gam)$p.coef[1]+1.96*summary(mod_farmland_species_time$gam)$se[1]
nb_sp_farmland_2017_up <- (summary(mod_farmland_species_time$gam)$p.coef[2]+1.96*summary(mod_farmland_species_time$gam)$se[2])*2017 + summary(mod_farmland_species_time$gam)$p.coef[1]-1.96*summary(mod_farmland_species_time$gam)$se[1]

## nb_association
nb_asso_obs_farmland_2017 <- summary(mod_farmland_asso_obs_time$gam)$p.coef[2]*2017 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
nb_asso_obs_farmland_2017_low <- (summary(mod_farmland_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_farmland_asso_obs_time$gam)$se[2])*2017 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_farmland_asso_obs_time$gam)$se[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
nb_asso_obs_farmland_2017_up <- (summary(mod_farmland_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_farmland_asso_obs_time$gam)$se[2])*2017 + summary(mod_farmland_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_farmland_asso_obs_time$gam)$se[1] + summary(mod_farmland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))

## connectance
connectance_farmland_2017 <- summary(mod_farmland_connectance_time$gam)$p.coef[2]*2017 + summary(mod_farmland_connectance_time$gam)$p.coef[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
connectance_farmland_2017_low <- (summary(mod_farmland_connectance_time$gam)$p.coef[2]-1.96*summary(mod_farmland_connectance_time$gam)$se[2])*2017 + summary(mod_farmland_connectance_time$gam)$p.coef[1]+1.96*summary(mod_farmland_connectance_time$gam)$se[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))
connectance_farmland_2017_up <- (summary(mod_farmland_connectance_time$gam)$p.coef[2]+1.96*summary(mod_farmland_connectance_time$gam)$se[2])*2017 + summary(mod_farmland_connectance_time$gam)$p.coef[1]-1.96*summary(mod_farmland_connectance_time$gam)$se[1] + summary(mod_farmland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_farmland$nb_sp))

## modularity
modularity_farmland_2017 <- summary(mod_farmland_modularity_time$gam)$p.coef[2]*2017 + summary(mod_farmland_modularity_time$gam)$p.coef[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
modularity_farmland_2017_low <- (summary(mod_farmland_modularity_time$gam)$p.coef[2]-1.96*summary(mod_farmland_modularity_time$gam)$se[2])*2017 + summary(mod_farmland_modularity_time$gam)$p.coef[1]+1.96*summary(mod_farmland_modularity_time$gam)$se[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
modularity_farmland_2017_up <- (summary(mod_farmland_modularity_time$gam)$p.coef[2]+1.96*summary(mod_farmland_modularity_time$gam)$se[2])*2017 + summary(mod_farmland_modularity_time$gam)$p.coef[1]-1.96*summary(mod_farmland_modularity_time$gam)$se[1] + summary(mod_farmland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))

## evenness
evenness_farmland_2017 <- summary(mod_farmland_evenness_time$gam)$p.coef[2]*2017 + summary(mod_farmland_evenness_time$gam)$p.coef[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
evenness_farmland_2017_low <- (summary(mod_farmland_evenness_time$gam)$p.coef[2]-1.96*summary(mod_farmland_evenness_time$gam)$se[2])*2017 + summary(mod_farmland_evenness_time$gam)$p.coef[1]+1.96*summary(mod_farmland_evenness_time$gam)$se[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))
evenness_farmland_2017_up <- (summary(mod_farmland_evenness_time$gam)$p.coef[2]+1.96*summary(mod_farmland_evenness_time$gam)$se[2])*2017 + summary(mod_farmland_evenness_time$gam)$p.coef[1]-1.96*summary(mod_farmland_evenness_time$gam)$se[1] + summary(mod_farmland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_farmland$connectance))

```

### Forest

```{r}

# models for forest

data_forest <- network_structure_square[which(network_structure_square$habitat=="Forest"),]

mod_forest_connectance_time <- gamm(connectance~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_forest, random = list(code_square=~1))
summary(mod_forest_connectance_time$gam)

mod_forest_asso_obs_time <- gamm(nb_association~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_forest, random = list(code_square=~1))
summary(mod_forest_asso_obs_time$gam)

mod_forest_modularity_time <- gamm(modularity~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_forest, random = list(code_square=~1))
summary(mod_forest_modularity_time$gam)

mod_forest_evenness_time <- gamm(evenness~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_forest, random = list(code_square=~1))
summary(mod_forest_evenness_time$gam)

mod_forest_species_time <- gamm(nb_sp~year+te(lon2,lat2,bs="tp",k=3), data=data_forest, random = list(code_square=~1))
summary(mod_forest_species_time$gam)

# network structure on average in 2001 in forest

## nb_sp
nb_sp_forest_2001 <- summary(mod_forest_species_time$gam)$p.coef[2]*2001 + summary(mod_forest_species_time$gam)$p.coef[1]
nb_sp_forest_2001_up <- (summary(mod_forest_species_time$gam)$p.coef[2]-1.96*summary(mod_forest_species_time$gam)$se[2])*2001 + summary(mod_forest_species_time$gam)$p.coef[1]+1.96*summary(mod_forest_species_time$gam)$se[1]
nb_sp_forest_2001_low <- (summary(mod_forest_species_time$gam)$p.coef[2]+1.96*summary(mod_forest_species_time$gam)$se[2])*2001 + summary(mod_forest_species_time$gam)$p.coef[1]-1.96*summary(mod_forest_species_time$gam)$se[1]

## nb_association
nb_asso_obs_forest_2001 <- summary(mod_forest_asso_obs_time$gam)$p.coef[2]*2001 + summary(mod_forest_asso_obs_time$gam)$p.coef[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
nb_asso_obs_forest_2001_up <- (summary(mod_forest_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_forest_asso_obs_time$gam)$se[2])*2001 + summary(mod_forest_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_forest_asso_obs_time$gam)$se[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
nb_asso_obs_forest_2001_low <- (summary(mod_forest_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_forest_asso_obs_time$gam)$se[2])*2001 + summary(mod_forest_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_forest_asso_obs_time$gam)$se[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))

## connectance
connectance_forest_2001 <- summary(mod_forest_connectance_time$gam)$p.coef[2]*2001 + summary(mod_forest_connectance_time$gam)$p.coef[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
connectance_forest_2001_up <- (summary(mod_forest_connectance_time$gam)$p.coef[2]-1.96*summary(mod_forest_connectance_time$gam)$se[2])*2001 + summary(mod_forest_connectance_time$gam)$p.coef[1]+1.96*summary(mod_forest_connectance_time$gam)$se[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
connectance_forest_2001_low <- (summary(mod_forest_connectance_time$gam)$p.coef[2]+1.96*summary(mod_forest_connectance_time$gam)$se[2])*2001 + summary(mod_forest_connectance_time$gam)$p.coef[1]-1.96*summary(mod_forest_connectance_time$gam)$se[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))

## modularity
modularity_forest_2001 <- summary(mod_forest_modularity_time$gam)$p.coef[2]*2001 + summary(mod_forest_modularity_time$gam)$p.coef[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
modularity_forest_2001_up <- (summary(mod_forest_modularity_time$gam)$p.coef[2]-1.96*summary(mod_forest_modularity_time$gam)$se[2])*2001 + summary(mod_forest_modularity_time$gam)$p.coef[1]+1.96*summary(mod_forest_modularity_time$gam)$se[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
modularity_forest_2001_low <- (summary(mod_forest_modularity_time$gam)$p.coef[2]+1.96*summary(mod_forest_modularity_time$gam)$se[2])*2001 + summary(mod_forest_modularity_time$gam)$p.coef[1]-1.96*summary(mod_forest_modularity_time$gam)$se[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))

## evenness
evenness_forest_2001 <- summary(mod_forest_evenness_time$gam)$p.coef[2]*2001 + summary(mod_forest_evenness_time$gam)$p.coef[1] + summary(mod_forest_evenness_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
evenness_forest_2001_up <- (summary(mod_forest_evenness_time$gam)$p.coef[2]-1.96*summary(mod_forest_evenness_time$gam)$se[2])*2001 + summary(mod_forest_evenness_time$gam)$p.coef[1]+1.96*summary(mod_forest_evenness_time$gam)$se[1] + summary(mod_forest_evenness_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
evenness_forest_2001_low <- (summary(mod_forest_evenness_time$gam)$p.coef[2]+1.96*summary(mod_forest_evenness_time$gam)$se[2])*2001 + summary(mod_forest_evenness_time$gam)$p.coef[1]-1.96*summary(mod_forest_evenness_time$gam)$se[1] + summary(mod_forest_evenness_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))


# network structure on average in 2017 in forest

## nb_sp
nb_sp_forest_2017 <- summary(mod_forest_species_time$gam)$p.coef[2]*2017 + summary(mod_forest_species_time$gam)$p.coef[1]
nb_sp_forest_2017_low <- (summary(mod_forest_species_time$gam)$p.coef[2]-1.96*summary(mod_forest_species_time$gam)$se[2])*2017 + summary(mod_forest_species_time$gam)$p.coef[1]+1.96*summary(mod_forest_species_time$gam)$se[1]
nb_sp_forest_2017_up <- (summary(mod_forest_species_time$gam)$p.coef[2]+1.96*summary(mod_forest_species_time$gam)$se[2])*2017 + summary(mod_forest_species_time$gam)$p.coef[1]-1.96*summary(mod_forest_species_time$gam)$se[1]

## nb_association
nb_asso_obs_forest_2017 <- summary(mod_forest_asso_obs_time$gam)$p.coef[2]*2017 + summary(mod_forest_asso_obs_time$gam)$p.coef[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
nb_asso_obs_forest_2017_low <- (summary(mod_forest_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_forest_asso_obs_time$gam)$se[2])*2017 + summary(mod_forest_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_forest_asso_obs_time$gam)$se[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
nb_asso_obs_forest_2017_up <- (summary(mod_forest_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_forest_asso_obs_time$gam)$se[2])*2017 + summary(mod_forest_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_forest_asso_obs_time$gam)$se[1] + summary(mod_forest_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))

## connectance
connectance_forest_2017 <- summary(mod_forest_connectance_time$gam)$p.coef[2]*2017 + summary(mod_forest_connectance_time$gam)$p.coef[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
connectance_forest_2017_low <- (summary(mod_forest_connectance_time$gam)$p.coef[2]-1.96*summary(mod_forest_connectance_time$gam)$se[2])*2017 + summary(mod_forest_connectance_time$gam)$p.coef[1]+1.96*summary(mod_forest_connectance_time$gam)$se[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))
connectance_forest_2017_up <- (summary(mod_forest_connectance_time$gam)$p.coef[2]+1.96*summary(mod_forest_connectance_time$gam)$se[2])*2017 + summary(mod_forest_connectance_time$gam)$p.coef[1]-1.96*summary(mod_forest_connectance_time$gam)$se[1] + summary(mod_forest_connectance_time$gam)$p.coef[3]*mean(na.omit(data_forest$nb_sp))

## modularity
modularity_forest_2017 <- summary(mod_forest_modularity_time$gam)$p.coef[2]*2017 + summary(mod_forest_modularity_time$gam)$p.coef[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
modularity_forest_2017_low <- (summary(mod_forest_modularity_time$gam)$p.coef[2]-1.96*summary(mod_forest_modularity_time$gam)$se[2])*2017 + summary(mod_forest_modularity_time$gam)$p.coef[1]+1.96*summary(mod_forest_modularity_time$gam)$se[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
modularity_forest_2017_up <- (summary(mod_forest_modularity_time$gam)$p.coef[2]+1.96*summary(mod_forest_modularity_time$gam)$se[2])*2017 + summary(mod_forest_modularity_time$gam)$p.coef[1]-1.96*summary(mod_forest_modularity_time$gam)$se[1] + summary(mod_forest_modularity_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))

## evenness
evenness_forest_2017 <- summary(mod_forest_evenness_time$gam)$p.coef[2]*2017 + summary(mod_forest_evenness_time$gam)$p.coef[1] + summary(mod_forest_evenness_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
evenness_forest_2017_low <- (summary(mod_forest_evenness_time$gam)$p.coef[2]-1.96*summary(mod_forest_evenness_time$gam)$se[2])*2017 + summary(mod_forest_evenness_time$gam)$p.coef[1]+1.96*summary(mod_forest_evenness_time$gam)$se[1] + summary(mod_forest_evenness_time$gam)$p.coef[3]*mean(na.omit(data_forest$connectance))
evenness_forest_2017_up <- (summary(mod_forest_evenness_time$gam)$p.coef[2]+1.96*summary(mod_forest_evenness_time$gam)$se[2])*2017 + summary(mod_forest_evenness_time$gam)$p.coef[1]-1.96*summary(mod_forest_evenness_time$gam)$se[1] + summary(mod_forest_evenness_time$gam)$p.coef[4]*mean(na.omit(data_forest$connectance))

```

### Urban

```{r}

# models for urban

data_urban <- network_structure_square[which(network_structure_square$habitat=="Urban"),]

mod_urban_connectance_time <- gamm(connectance~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_urban, random = list(code_square=~1))
summary(mod_urban_connectance_time$gam)

mod_urban_asso_obs_time <- gamm(nb_association~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_urban, random = list(code_square=~1))
summary(mod_urban_asso_obs_time$gam)

mod_urban_modularity_time <- gamm(modularity~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_urban, random = list(code_square=~1))
summary(mod_urban_modularity_time$gam)

mod_urban_evenness_time <- gamm(evenness~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_urban, random = list(code_square=~1))
summary(mod_urban_evenness_time$gam)

mod_urban_species_time <- gamm(nb_sp~year+te(lon2,lat2,bs="tp",k=3), data=data_urban, random = list(code_square=~1))
summary(mod_urban_species_time$gam)

# network structure on average in 2001 in urban

## nb_sp
nb_sp_urban_2001 <- summary(mod_urban_species_time$gam)$p.coef[2]*2001 + summary(mod_urban_species_time$gam)$p.coef[1]
nb_sp_urban_2001_up <- (summary(mod_urban_species_time$gam)$p.coef[2]-1.96*summary(mod_urban_species_time$gam)$se[2])*2001 + summary(mod_urban_species_time$gam)$p.coef[1]+1.96*summary(mod_urban_species_time$gam)$se[1]
nb_sp_urban_2001_low <- (summary(mod_urban_species_time$gam)$p.coef[2]+1.96*summary(mod_urban_species_time$gam)$se[2])*2001 + summary(mod_urban_species_time$gam)$p.coef[1]-1.96*summary(mod_urban_species_time$gam)$se[1]

## nb_association
nb_asso_obs_urban_2001 <- summary(mod_urban_asso_obs_time$gam)$p.coef[2]*2001 + summary(mod_urban_asso_obs_time$gam)$p.coef[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
nb_asso_obs_urban_2001_up <- (summary(mod_urban_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_urban_asso_obs_time$gam)$se[2])*2001 + summary(mod_urban_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_urban_asso_obs_time$gam)$se[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
nb_asso_obs_urban_2001_low <- (summary(mod_urban_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_urban_asso_obs_time$gam)$se[2])*2001 + summary(mod_urban_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_urban_asso_obs_time$gam)$se[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))

## connectance
connectance_urban_2001 <- summary(mod_urban_connectance_time$gam)$p.coef[2]*2001 + summary(mod_urban_connectance_time$gam)$p.coef[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
connectance_urban_2001_up <- (summary(mod_urban_connectance_time$gam)$p.coef[2]-1.96*summary(mod_urban_connectance_time$gam)$se[2])*2001 + summary(mod_urban_connectance_time$gam)$p.coef[1]+1.96*summary(mod_urban_connectance_time$gam)$se[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
connectance_urban_2001_low <- (summary(mod_urban_connectance_time$gam)$p.coef[2]+1.96*summary(mod_urban_connectance_time$gam)$se[2])*2001 + summary(mod_urban_connectance_time$gam)$p.coef[1]-1.96*summary(mod_urban_connectance_time$gam)$se[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))

## modularity
modularity_urban_2001 <- summary(mod_urban_modularity_time$gam)$p.coef[2]*2001 + summary(mod_urban_modularity_time$gam)$p.coef[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
modularity_urban_2001_up <- (summary(mod_urban_modularity_time$gam)$p.coef[2]-1.96*summary(mod_urban_modularity_time$gam)$se[2])*2001 + summary(mod_urban_modularity_time$gam)$p.coef[1]+1.96*summary(mod_urban_modularity_time$gam)$se[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
modularity_urban_2001_low <- (summary(mod_urban_modularity_time$gam)$p.coef[2]+1.96*summary(mod_urban_modularity_time$gam)$se[2])*2001 + summary(mod_urban_modularity_time$gam)$p.coef[1]-1.96*summary(mod_urban_modularity_time$gam)$se[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))

## evenness
evenness_urban_2001 <- summary(mod_urban_evenness_time$gam)$p.coef[2]*2001 + summary(mod_urban_evenness_time$gam)$p.coef[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
evenness_urban_2001_up <- (summary(mod_urban_evenness_time$gam)$p.coef[2]-1.96*summary(mod_urban_evenness_time$gam)$se[2])*2001 + summary(mod_urban_evenness_time$gam)$p.coef[1]+1.96*summary(mod_urban_evenness_time$gam)$se[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
evenness_urban_2001_low <- (summary(mod_urban_evenness_time$gam)$p.coef[2]+1.96*summary(mod_urban_evenness_time$gam)$se[2])*2001 + summary(mod_urban_evenness_time$gam)$p.coef[1]-1.96*summary(mod_urban_evenness_time$gam)$se[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))


# network structure on average in 2017 in urban

## nb_sp
nb_sp_urban_2017 <- summary(mod_urban_species_time$gam)$p.coef[2]*2017 + summary(mod_urban_species_time$gam)$p.coef[1]
nb_sp_urban_2017_low <- (summary(mod_urban_species_time$gam)$p.coef[2]-1.96*summary(mod_urban_species_time$gam)$se[2])*2017 + summary(mod_urban_species_time$gam)$p.coef[1]+1.96*summary(mod_urban_species_time$gam)$se[1]
nb_sp_urban_2017_up <- (summary(mod_urban_species_time$gam)$p.coef[2]+1.96*summary(mod_urban_species_time$gam)$se[2])*2017 + summary(mod_urban_species_time$gam)$p.coef[1]-1.96*summary(mod_urban_species_time$gam)$se[1]

## nb_association
nb_asso_obs_urban_2017 <- summary(mod_urban_asso_obs_time$gam)$p.coef[2]*2017 + summary(mod_urban_asso_obs_time$gam)$p.coef[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
nb_asso_obs_urban_2017_low <- (summary(mod_urban_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_urban_asso_obs_time$gam)$se[2])*2017 + summary(mod_urban_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_urban_asso_obs_time$gam)$se[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
nb_asso_obs_urban_2017_up <- (summary(mod_urban_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_urban_asso_obs_time$gam)$se[2])*2017 + summary(mod_urban_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_urban_asso_obs_time$gam)$se[1] + summary(mod_urban_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))

## connectance
connectance_urban_2017 <- summary(mod_urban_connectance_time$gam)$p.coef[2]*2017 + summary(mod_urban_connectance_time$gam)$p.coef[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
connectance_urban_2017_low <- (summary(mod_urban_connectance_time$gam)$p.coef[2]-1.96*summary(mod_urban_connectance_time$gam)$se[2])*2017 + summary(mod_urban_connectance_time$gam)$p.coef[1]+1.96*summary(mod_urban_connectance_time$gam)$se[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))
connectance_urban_2017_up <- (summary(mod_urban_connectance_time$gam)$p.coef[2]+1.96*summary(mod_urban_connectance_time$gam)$se[2])*2017 + summary(mod_urban_connectance_time$gam)$p.coef[1]-1.96*summary(mod_urban_connectance_time$gam)$se[1] + summary(mod_urban_connectance_time$gam)$p.coef[3]*mean(na.omit(data_urban$nb_sp))

## modularity
modularity_urban_2017 <- summary(mod_urban_modularity_time$gam)$p.coef[2]*2017 + summary(mod_urban_modularity_time$gam)$p.coef[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
modularity_urban_2017_low <- (summary(mod_urban_modularity_time$gam)$p.coef[2]-1.96*summary(mod_urban_modularity_time$gam)$se[2])*2017 + summary(mod_urban_modularity_time$gam)$p.coef[1]+1.96*summary(mod_urban_modularity_time$gam)$se[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
modularity_urban_2017_up <- (summary(mod_urban_modularity_time$gam)$p.coef[2]+1.96*summary(mod_urban_modularity_time$gam)$se[2])*2017 + summary(mod_urban_modularity_time$gam)$p.coef[1]-1.96*summary(mod_urban_modularity_time$gam)$se[1] + summary(mod_urban_modularity_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))

## evenness
evenness_urban_2017 <- summary(mod_urban_evenness_time$gam)$p.coef[2]*2017 + summary(mod_urban_evenness_time$gam)$p.coef[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
evenness_urban_2017_low <- (summary(mod_urban_evenness_time$gam)$p.coef[2]-1.96*summary(mod_urban_evenness_time$gam)$se[2])*2017 + summary(mod_urban_evenness_time$gam)$p.coef[1]+1.96*summary(mod_urban_evenness_time$gam)$se[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))
evenness_urban_2017_up <- (summary(mod_urban_evenness_time$gam)$p.coef[2]+1.96*summary(mod_urban_evenness_time$gam)$se[2])*2017 + summary(mod_urban_evenness_time$gam)$p.coef[1]-1.96*summary(mod_urban_evenness_time$gam)$se[1] + summary(mod_urban_evenness_time$gam)$p.coef[3]*mean(na.omit(data_urban$connectance))

```

### Natural openland

```{r}

# models for natural_openland

data_natural_openland <- network_structure_square[which(network_structure_square$habitat=="Natural openland"),]

mod_natural_openland_connectance_time <- gamm(connectance~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_natural_openland, random = list(code_square=~1))
summary(mod_natural_openland_connectance_time$gam)

mod_natural_openland_asso_obs_time <- gamm(nb_association~year+nb_sp+te(lon2,lat2,bs="tp",k=3), data=data_natural_openland, random = list(code_square=~1))
summary(mod_natural_openland_asso_obs_time$gam)

mod_natural_openland_modularity_time <- gamm(modularity~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_natural_openland, random = list(code_square=~1))
summary(mod_natural_openland_modularity_time$gam)

mod_natural_openland_evenness_time <- gamm(evenness~year+connectance+te(lon2,lat2,bs="tp",k=3), data=data_natural_openland, random = list(code_square=~1))
summary(mod_natural_openland_evenness_time$gam)

mod_natural_openland_species_time <- gamm(nb_sp~year+te(lon2,lat2,bs="tp",k=3), data=data_natural_openland, random = list(code_square=~1))
summary(mod_natural_openland_species_time$gam)

# network structure on average in 2001 in natural_openland

## nb_sp
nb_sp_natural_openland_2001 <- summary(mod_natural_openland_species_time$gam)$p.coef[2]*2001 + summary(mod_natural_openland_species_time$gam)$p.coef[1]
nb_sp_natural_openland_2001_up <- (summary(mod_natural_openland_species_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_species_time$gam)$se[2])*2001 + summary(mod_natural_openland_species_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_species_time$gam)$se[1]
nb_sp_natural_openland_2001_low <- (summary(mod_natural_openland_species_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_species_time$gam)$se[2])*2001 + summary(mod_natural_openland_species_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_species_time$gam)$se[1]

## nb_association
nb_asso_obs_natural_openland_2001 <- summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]*2001 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
nb_asso_obs_natural_openland_2001_up <- (summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[2])*2001 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
nb_asso_obs_natural_openland_2001_low <- (summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[2])*2001 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))

## connectance
connectance_natural_openland_2001 <- summary(mod_natural_openland_connectance_time$gam)$p.coef[2]*2001 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
connectance_natural_openland_2001_up <- (summary(mod_natural_openland_connectance_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_connectance_time$gam)$se[2])*2001 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_connectance_time$gam)$se[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
connectance_natural_openland_2001_low <- (summary(mod_natural_openland_connectance_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_connectance_time$gam)$se[2])*2001 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_connectance_time$gam)$se[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))

## modularity
modularity_natural_openland_2001 <- summary(mod_natural_openland_modularity_time$gam)$p.coef[2]*2001 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
modularity_natural_openland_2001_up <- (summary(mod_natural_openland_modularity_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_modularity_time$gam)$se[2])*2001 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_modularity_time$gam)$se[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
modularity_natural_openland_2001_low <- (summary(mod_natural_openland_modularity_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_modularity_time$gam)$se[2])*2001 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_modularity_time$gam)$se[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))

## evenness
evenness_natural_openland_2001 <- summary(mod_natural_openland_evenness_time$gam)$p.coef[2]*2001 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
evenness_natural_openland_2001_up <- (summary(mod_natural_openland_evenness_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_evenness_time$gam)$se[2])*2001 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_evenness_time$gam)$se[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
evenness_natural_openland_2001_low <- (summary(mod_natural_openland_evenness_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_evenness_time$gam)$se[2])*2001 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_evenness_time$gam)$se[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))


# network structure on average in 2017 in natural_openland

## nb_sp
nb_sp_natural_openland_2017 <- summary(mod_natural_openland_species_time$gam)$p.coef[2]*2017 + summary(mod_natural_openland_species_time$gam)$p.coef[1]
nb_sp_natural_openland_2017_low <- (summary(mod_natural_openland_species_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_species_time$gam)$se[2])*2017 + summary(mod_natural_openland_species_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_species_time$gam)$se[1]
nb_sp_natural_openland_2017_up <- (summary(mod_natural_openland_species_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_species_time$gam)$se[2])*2017 + summary(mod_natural_openland_species_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_species_time$gam)$se[1]

## nb_association
nb_asso_obs_natural_openland_2017 <- summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]*2017 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
nb_asso_obs_natural_openland_2017_low <- (summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[2])*2017 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
nb_asso_obs_natural_openland_2017_up <- (summary(mod_natural_openland_asso_obs_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[2])*2017 + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_asso_obs_time$gam)$se[1] + summary(mod_natural_openland_asso_obs_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))

## connectance
connectance_natural_openland_2017 <- summary(mod_natural_openland_connectance_time$gam)$p.coef[2]*2017 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
connectance_natural_openland_2017_low <- (summary(mod_natural_openland_connectance_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_connectance_time$gam)$se[2])*2017 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_connectance_time$gam)$se[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))
connectance_natural_openland_2017_up <- (summary(mod_natural_openland_connectance_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_connectance_time$gam)$se[2])*2017 + summary(mod_natural_openland_connectance_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_connectance_time$gam)$se[1] + summary(mod_natural_openland_connectance_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$nb_sp))

## modularity
modularity_natural_openland_2017 <- summary(mod_natural_openland_modularity_time$gam)$p.coef[2]*2017 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
modularity_natural_openland_2017_low <- (summary(mod_natural_openland_modularity_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_modularity_time$gam)$se[2])*2017 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_modularity_time$gam)$se[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
modularity_natural_openland_2017_up <- (summary(mod_natural_openland_modularity_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_modularity_time$gam)$se[2])*2017 + summary(mod_natural_openland_modularity_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_modularity_time$gam)$se[1] + summary(mod_natural_openland_modularity_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))

## evenness
evenness_natural_openland_2017 <- summary(mod_natural_openland_evenness_time$gam)$p.coef[2]*2017 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
evenness_natural_openland_2017_low <- (summary(mod_natural_openland_evenness_time$gam)$p.coef[2]-1.96*summary(mod_natural_openland_evenness_time$gam)$se[2])*2017 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1]+1.96*summary(mod_natural_openland_evenness_time$gam)$se[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))
evenness_natural_openland_2017_up <- (summary(mod_natural_openland_evenness_time$gam)$p.coef[2]+1.96*summary(mod_natural_openland_evenness_time$gam)$se[2])*2017 + summary(mod_natural_openland_evenness_time$gam)$p.coef[1]-1.96*summary(mod_natural_openland_evenness_time$gam)$se[1] + summary(mod_natural_openland_evenness_time$gam)$p.coef[3]*mean(na.omit(data_natural_openland$connectance))

```

## Plot results

### Number of species

```{r}

dd <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_2001,nb_sp_2017))
new_dd <- data.frame(year=c(2001:2017))
mod_nb_sp <- lm(nb_sp ~ year,dd)
new_dd$pred <- predict(mod_nb_sp, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_2001_low,nb_sp_2017_low))
mod_nb_sp <- lm(nb_sp ~ year,dd)
new_dd$lwr <- predict(mod_nb_sp, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_2001_up,nb_sp_2017_up))
mod_nb_sp <- lm(nb_sp ~ year,dd)
new_dd$upr <- predict(mod_nb_sp, newdata=new_dd)

dd_farm <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_farmland_2001,nb_sp_farmland_2017))
new_dd_farm <- data.frame(year=c(2001:2017))
mod_nb_sp <- lm(nb_sp ~ year,dd_farm)
new_dd_farm$pred <- predict(mod_nb_sp, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_farmland_2001_low,nb_sp_farmland_2017_low))
mod_nb_sp <- lm(nb_sp ~ year,dd_farm)
new_dd_farm$lwr <- predict(mod_nb_sp, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_farmland_2001_up,nb_sp_farmland_2017_up))
mod_nb_sp <- lm(nb_sp ~ year,dd_farm)
new_dd_farm$upr <- predict(mod_nb_sp, newdata=new_dd_farm)

dd_forest <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_forest_2001,nb_sp_forest_2017))
new_dd_forest <- data.frame(year=c(2001:2017))
mod_nb_sp <- lm(nb_sp ~ year,dd_forest)
new_dd_forest$pred <- predict(mod_nb_sp, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_forest_2001_low,nb_sp_forest_2017_low))
mod_nb_sp <- lm(nb_sp ~ year,dd_forest)
new_dd_forest$lwr <- predict(mod_nb_sp, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_forest_2001_up,nb_sp_forest_2017_up))
mod_nb_sp <- lm(nb_sp ~ year,dd_forest)
new_dd_forest$upr <- predict(mod_nb_sp, newdata=new_dd_forest)

dd_natural_open <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_natural_openland_2001,nb_sp_natural_openland_2017))
new_dd_natural_open <- data.frame(year=c(2001:2017))
mod_nb_sp <- lm(nb_sp ~ year,dd_natural_open)
new_dd_natural_open$pred <- predict(mod_nb_sp, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_natural_openland_2001_low,nb_sp_natural_openland_2017_low))
mod_nb_sp <- lm(nb_sp ~ year,dd_natural_open)
new_dd_natural_open$lwr <- predict(mod_nb_sp, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_natural_openland_2001_up,nb_sp_natural_openland_2017_up))
mod_nb_sp <- lm(nb_sp ~ year,dd_natural_open)
new_dd_natural_open$upr <- predict(mod_nb_sp, newdata=new_dd_natural_open)

dd_urban <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_urban_2001,nb_sp_urban_2017))
new_dd_urban <- data.frame(year=c(2001:2017))
mod_nb_sp <- lm(nb_sp ~ year,dd_urban)
new_dd_urban$pred <- predict(mod_nb_sp, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_urban_2001_low,nb_sp_urban_2017_low))
mod_nb_sp <- lm(nb_sp ~ year,dd_urban)
new_dd_urban$lwr <- predict(mod_nb_sp, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),nb_sp=c(nb_sp_urban_2001_up,nb_sp_urban_2017_up))
mod_nb_sp <- lm(nb_sp ~ year,dd_urban)
new_dd_urban$upr <- predict(mod_nb_sp, newdata=new_dd_urban)

ggplot(new_dd,aes(x=year,y=pred)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr),fill = "black", alpha=0.2) +
  geom_line() +
  geom_ribbon(data=new_dd_farm,aes(ymin=lwr, ymax=upr),fill = "#f5b041", alpha=0.2) +
  geom_line(data=new_dd_farm, col="#f5b041") +
  geom_ribbon(data=new_dd_forest,aes(ymin=lwr, ymax=upr),fill = "#52be80", alpha=0.2) +
  geom_line(data=new_dd_forest, col="#52be80") +
  geom_ribbon(data=new_dd_natural_open,aes(ymin=lwr, ymax=upr),fill = "#2C65FD", alpha=0.2) +
  geom_line(data=new_dd_natural_open, col="#2C65FD") +
  geom_ribbon(data=new_dd_urban,aes(ymin=lwr, ymax=upr),fill = "#FF4001", alpha=0.2) +
  geom_line(data=new_dd_urban, col="#FF4001") +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Number of species by network") + 
  theme_modern()
  

ggplot(new_dd[which(new_dd$year %in% c(2001,2017)),],aes(x=as.factor(year),y=pred,group=year)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(0.05), alpha=0.2) +
  geom_point(size=3,data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],col = "#f5b041") +
  geom_errorbar(data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#f5b041", alpha=0.2) +
  geom_point(size=3,data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],col = "#52be80") +
  geom_errorbar(data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#52be80", alpha=0.2) +
  geom_point(size=3,data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],col = "#2C65FD") +
  geom_errorbar(data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#2C65FD", alpha=0.2) +
  geom_point(size=3,data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],col = "#FF4001") +
  geom_errorbar(data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#FF4001", alpha=0.2) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Number of species by network") + 
  theme_modern()

```

### Number of association

```{r}

dd <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_2001,nb_asso_obs_2017))
new_dd <- data.frame(year=c(2001:2017))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd)
new_dd$pred <- predict(mod_nb_asso_obs, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_2001_low,nb_asso_obs_2017_low))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd)
new_dd$lwr <- predict(mod_nb_asso_obs, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_2001_up,nb_asso_obs_2017_up))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd)
new_dd$upr <- predict(mod_nb_asso_obs, newdata=new_dd)

dd_farm <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_farmland_2001,nb_asso_obs_farmland_2017))
new_dd_farm <- data.frame(year=c(2001:2017))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_farm)
new_dd_farm$pred <- predict(mod_nb_asso_obs, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_farmland_2001_low,nb_asso_obs_farmland_2017_low))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_farm)
new_dd_farm$lwr <- predict(mod_nb_asso_obs, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_farmland_2001_up,nb_asso_obs_farmland_2017_up))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_farm)
new_dd_farm$upr <- predict(mod_nb_asso_obs, newdata=new_dd_farm)

dd_forest <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_forest_2001,nb_asso_obs_forest_2017))
new_dd_forest <- data.frame(year=c(2001:2017))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_forest)
new_dd_forest$pred <- predict(mod_nb_asso_obs, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_forest_2001_low,nb_asso_obs_forest_2017_low))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_forest)
new_dd_forest$lwr <- predict(mod_nb_asso_obs, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_forest_2001_up,nb_asso_obs_forest_2017_up))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_forest)
new_dd_forest$upr <- predict(mod_nb_asso_obs, newdata=new_dd_forest)

dd_natural_open <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_natural_openland_2001,nb_asso_obs_natural_openland_2017))
new_dd_natural_open <- data.frame(year=c(2001:2017))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_natural_open)
new_dd_natural_open$pred <- predict(mod_nb_asso_obs, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_natural_openland_2001_low,nb_asso_obs_natural_openland_2017_low))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_natural_open)
new_dd_natural_open$lwr <- predict(mod_nb_asso_obs, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_natural_openland_2001_up,nb_asso_obs_natural_openland_2017_up))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_natural_open)
new_dd_natural_open$upr <- predict(mod_nb_asso_obs, newdata=new_dd_natural_open)

dd_urban <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_urban_2001,nb_asso_obs_urban_2017))
new_dd_urban <- data.frame(year=c(2001:2017))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_urban)
new_dd_urban$pred <- predict(mod_nb_asso_obs, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_urban_2001_low,nb_asso_obs_urban_2017_low))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_urban)
new_dd_urban$lwr <- predict(mod_nb_asso_obs, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),nb_asso_obs=c(nb_asso_obs_urban_2001_up,nb_asso_obs_urban_2017_up))
mod_nb_asso_obs <- lm(nb_asso_obs ~ year,dd_urban)
new_dd_urban$upr <- predict(mod_nb_asso_obs, newdata=new_dd_urban)

ggplot(new_dd,aes(x=year,y=pred)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr),fill = "black", alpha=0.2) +
  geom_line() +
  geom_ribbon(data=new_dd_farm,aes(ymin=lwr, ymax=upr),fill = "#f5b041", alpha=0.2) +
  geom_line(data=new_dd_farm, col="#f5b041") +
  geom_ribbon(data=new_dd_forest,aes(ymin=lwr, ymax=upr),fill = "#52be80", alpha=0.2) +
  geom_line(data=new_dd_forest, col="#52be80") +
  geom_ribbon(data=new_dd_natural_open,aes(ymin=lwr, ymax=upr),fill = "#2C65FD", alpha=0.2) +
  geom_line(data=new_dd_natural_open, col="#2C65FD") +
  geom_ribbon(data=new_dd_urban,aes(ymin=lwr, ymax=upr),fill = "#FF4001", alpha=0.2) +
  geom_line(data=new_dd_urban, col="#FF4001") +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Number of observed associations by network\ncontrolling for network size") + 
  theme_modern()
  
  
ggplot(new_dd[which(new_dd$year %in% c(2001,2017)),],aes(x=as.factor(year),y=pred,group=year)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(0.05), alpha=0.2) +
  geom_point(size=3,data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],col = "#f5b041") +
  geom_errorbar(data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#f5b041", alpha=0.2) +
  geom_point(size=3,data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],col = "#52be80") +
  geom_errorbar(data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#52be80", alpha=0.2) +
  geom_point(size=3,data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],col = "#2C65FD") +
  geom_errorbar(data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#2C65FD", alpha=0.2) +
  geom_point(size=3,data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],col = "#FF4001") +
  geom_errorbar(data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#FF4001", alpha=0.2) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") +  ylab("Number of observed associations by network\ncontrolling for network size") + 
  theme_modern()

```

### Connectance

```{r}

dd <- data.frame(year=c(2001,2017),connectance=c(connectance_2001,connectance_2017))
new_dd <- data.frame(year=c(2001:2017))
mod_connectance <- lm(connectance ~ year,dd)
new_dd$pred <- predict(mod_connectance, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),connectance=c(connectance_2001_low,connectance_2017_low))
mod_connectance <- lm(connectance ~ year,dd)
new_dd$lwr <- predict(mod_connectance, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),connectance=c(connectance_2001_up,connectance_2017_up))
mod_connectance <- lm(connectance ~ year,dd)
new_dd$upr <- predict(mod_connectance, newdata=new_dd)

dd_farm <- data.frame(year=c(2001,2017),connectance=c(connectance_farmland_2001,connectance_farmland_2017))
new_dd_farm <- data.frame(year=c(2001:2017))
mod_connectance <- lm(connectance ~ year,dd_farm)
new_dd_farm$pred <- predict(mod_connectance, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),connectance=c(connectance_farmland_2001_low,connectance_farmland_2017_low))
mod_connectance <- lm(connectance ~ year,dd_farm)
new_dd_farm$lwr <- predict(mod_connectance, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),connectance=c(connectance_farmland_2001_up,connectance_farmland_2017_up))
mod_connectance <- lm(connectance ~ year,dd_farm)
new_dd_farm$upr <- predict(mod_connectance, newdata=new_dd_farm)

dd_forest <- data.frame(year=c(2001,2017),connectance=c(connectance_forest_2001,connectance_forest_2017))
new_dd_forest <- data.frame(year=c(2001:2017))
mod_connectance <- lm(connectance ~ year,dd_forest)
new_dd_forest$pred <- predict(mod_connectance, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),connectance=c(connectance_forest_2001_low,connectance_forest_2017_low))
mod_connectance <- lm(connectance ~ year,dd_forest)
new_dd_forest$lwr <- predict(mod_connectance, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),connectance=c(connectance_forest_2001_up,connectance_forest_2017_up))
mod_connectance <- lm(connectance ~ year,dd_forest)
new_dd_forest$upr <- predict(mod_connectance, newdata=new_dd_forest)

dd_natural_open <- data.frame(year=c(2001,2017),connectance=c(connectance_natural_openland_2001,connectance_natural_openland_2017))
new_dd_natural_open <- data.frame(year=c(2001:2017))
mod_connectance <- lm(connectance ~ year,dd_natural_open)
new_dd_natural_open$pred <- predict(mod_connectance, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),connectance=c(connectance_natural_openland_2001_low,connectance_natural_openland_2017_low))
mod_connectance <- lm(connectance ~ year,dd_natural_open)
new_dd_natural_open$lwr <- predict(mod_connectance, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),connectance=c(connectance_natural_openland_2001_up,connectance_natural_openland_2017_up))
mod_connectance <- lm(connectance ~ year,dd_natural_open)
new_dd_natural_open$upr <- predict(mod_connectance, newdata=new_dd_natural_open)

dd_urban <- data.frame(year=c(2001,2017),connectance=c(connectance_urban_2001,connectance_urban_2017))
new_dd_urban <- data.frame(year=c(2001:2017))
mod_connectance <- lm(connectance ~ year,dd_urban)
new_dd_urban$pred <- predict(mod_connectance, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),connectance=c(connectance_urban_2001_low,connectance_urban_2017_low))
mod_connectance <- lm(connectance ~ year,dd_urban)
new_dd_urban$lwr <- predict(mod_connectance, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),connectance=c(connectance_urban_2001_up,connectance_urban_2017_up))
mod_connectance <- lm(connectance ~ year,dd_urban)
new_dd_urban$upr <- predict(mod_connectance, newdata=new_dd_urban)

ggplot(new_dd,aes(x=year,y=pred)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr),fill = "black", alpha=0.2) +
  geom_line() +
  geom_ribbon(data=new_dd_farm,aes(ymin=lwr, ymax=upr),fill = "#f5b041", alpha=0.2) +
  geom_line(data=new_dd_farm, col="#f5b041") +
  geom_ribbon(data=new_dd_forest,aes(ymin=lwr, ymax=upr),fill = "#52be80", alpha=0.2) +
  geom_line(data=new_dd_forest, col="#52be80") +
  geom_ribbon(data=new_dd_natural_open,aes(ymin=lwr, ymax=upr),fill = "#2C65FD", alpha=0.2) +
  geom_line(data=new_dd_natural_open, col="#2C65FD") +
  geom_ribbon(data=new_dd_urban,aes(ymin=lwr, ymax=upr),fill = "#FF4001", alpha=0.2) +
  geom_line(data=new_dd_urban, col="#FF4001") +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Network connectance\ncontrolling for network size") + 
  theme_modern()
  
ggplot(new_dd[which(new_dd$year %in% c(2001,2017)),],aes(x=as.factor(year),y=pred,group=year)) + 
  geom_point(size=3,data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],col = "#f5b041") +
  geom_errorbar(data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#f5b041", alpha=0.2) +
  geom_point(size=3,data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],col = "#52be80") +
  geom_errorbar(data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#52be80", alpha=0.2) +
  geom_point(size=3,data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],col = "#2C65FD") +
  geom_errorbar(data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#2C65FD", alpha=0.2) +
  geom_point(size=3,data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],col = "#FF4001") +
  geom_errorbar(data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#FF4001", alpha=0.2) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(0.05), alpha=0.2) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") +  ylab("Network connectance\ncontrolling for network size") + 
  theme_modern()

```

### Modularity

```{r}

dd <- data.frame(year=c(2001,2017),modularity=c(modularity_2001,modularity_2017))
new_dd <- data.frame(year=c(2001:2017))
mod_modularity <- lm(modularity ~ year,dd)
new_dd$pred <- predict(mod_modularity, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),modularity=c(modularity_2001_low,modularity_2017_low))
mod_modularity <- lm(modularity ~ year,dd)
new_dd$lwr <- predict(mod_modularity, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),modularity=c(modularity_2001_up,modularity_2017_up))
mod_modularity <- lm(modularity ~ year,dd)
new_dd$upr <- predict(mod_modularity, newdata=new_dd)

dd_farm <- data.frame(year=c(2001,2017),modularity=c(modularity_farmland_2001,modularity_farmland_2017))
new_dd_farm <- data.frame(year=c(2001:2017))
mod_modularity <- lm(modularity ~ year,dd_farm)
new_dd_farm$pred <- predict(mod_modularity, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),modularity=c(modularity_farmland_2001_low,modularity_farmland_2017_low))
mod_modularity <- lm(modularity ~ year,dd_farm)
new_dd_farm$lwr <- predict(mod_modularity, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),modularity=c(modularity_farmland_2001_up,modularity_farmland_2017_up))
mod_modularity <- lm(modularity ~ year,dd_farm)
new_dd_farm$upr <- predict(mod_modularity, newdata=new_dd_farm)

dd_forest <- data.frame(year=c(2001,2017),modularity=c(modularity_forest_2001,modularity_forest_2017))
new_dd_forest <- data.frame(year=c(2001:2017))
mod_modularity <- lm(modularity ~ year,dd_forest)
new_dd_forest$pred <- predict(mod_modularity, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),modularity=c(modularity_forest_2001_low,modularity_forest_2017_low))
mod_modularity <- lm(modularity ~ year,dd_forest)
new_dd_forest$lwr <- predict(mod_modularity, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),modularity=c(modularity_forest_2001_up,modularity_forest_2017_up))
mod_modularity <- lm(modularity ~ year,dd_forest)
new_dd_forest$upr <- predict(mod_modularity, newdata=new_dd_forest)

dd_natural_open <- data.frame(year=c(2001,2017),modularity=c(modularity_natural_openland_2001,modularity_natural_openland_2017))
new_dd_natural_open <- data.frame(year=c(2001:2017))
mod_modularity <- lm(modularity ~ year,dd_natural_open)
new_dd_natural_open$pred <- predict(mod_modularity, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),modularity=c(modularity_natural_openland_2001_low,modularity_natural_openland_2017_low))
mod_modularity <- lm(modularity ~ year,dd_natural_open)
new_dd_natural_open$lwr <- predict(mod_modularity, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),modularity=c(modularity_natural_openland_2001_up,modularity_natural_openland_2017_up))
mod_modularity <- lm(modularity ~ year,dd_natural_open)
new_dd_natural_open$upr <- predict(mod_modularity, newdata=new_dd_natural_open)

dd_urban <- data.frame(year=c(2001,2017),modularity=c(modularity_urban_2001,modularity_urban_2017))
new_dd_urban <- data.frame(year=c(2001:2017))
mod_modularity <- lm(modularity ~ year,dd_urban)
new_dd_urban$pred <- predict(mod_modularity, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),modularity=c(modularity_urban_2001_low,modularity_urban_2017_low))
mod_modularity <- lm(modularity ~ year,dd_urban)
new_dd_urban$lwr <- predict(mod_modularity, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),modularity=c(modularity_urban_2001_up,modularity_urban_2017_up))
mod_modularity <- lm(modularity ~ year,dd_urban)
new_dd_urban$upr <- predict(mod_modularity, newdata=new_dd_urban)

ggplot(new_dd,aes(x=year,y=pred)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr),fill = "black", alpha=0.2) +
  geom_line() +
  geom_ribbon(data=new_dd_farm,aes(ymin=lwr, ymax=upr),fill = "#f5b041", alpha=0.2) +
  geom_line(data=new_dd_farm, col="#f5b041") +
  geom_ribbon(data=new_dd_forest,aes(ymin=lwr, ymax=upr),fill = "#52be80", alpha=0.2) +
  geom_line(data=new_dd_forest, col="#52be80") +
  geom_ribbon(data=new_dd_natural_open,aes(ymin=lwr, ymax=upr),fill = "#2C65FD", alpha=0.2) +
  geom_line(data=new_dd_natural_open, col="#2C65FD") +
  geom_ribbon(data=new_dd_urban,aes(ymin=lwr, ymax=upr),fill = "#FF4001", alpha=0.2) +
  geom_line(data=new_dd_urban, col="#FF4001") +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Network modularity\ncontrolling for connectance")  +
  theme_modern()

ggplot(new_dd[which(new_dd$year %in% c(2001,2017)),],aes(x=as.factor(year),y=pred,group=year)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(0.05), alpha=0.2) +
  geom_point(size=3,data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],col = "#f5b041") +
  geom_errorbar(data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#f5b041", alpha=0.2) +
  geom_point(size=3,data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],col = "#52be80") +
  geom_errorbar(data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#52be80", alpha=0.2) +
  geom_point(size=3,data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],col = "#2C65FD") +
  geom_errorbar(data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#2C65FD", alpha=0.2) +
  geom_point(size=3,data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],col = "#FF4001") +
  geom_errorbar(data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#FF4001", alpha=0.2) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") +  ylab("Network modularity\ncontrolling for connectance") + 
  theme_modern()


```

### Evenness

```{r}

dd <- data.frame(year=c(2001,2017),evenness=c(evenness_2001,evenness_2017))
new_dd <- data.frame(year=c(2001:2017))
mod_evenness <- lm(evenness ~ year,dd)
new_dd$pred <- predict(mod_evenness, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),evenness=c(evenness_2001_low,evenness_2017_low))
mod_evenness <- lm(evenness ~ year,dd)
new_dd$lwr <- predict(mod_evenness, newdata=new_dd)
dd <- data.frame(year=c(2001,2017),evenness=c(evenness_2001_up,evenness_2017_up))
mod_evenness <- lm(evenness ~ year,dd)
new_dd$upr <- predict(mod_evenness, newdata=new_dd)

dd_farm <- data.frame(year=c(2001,2017),evenness=c(evenness_farmland_2001,evenness_farmland_2017))
new_dd_farm <- data.frame(year=c(2001:2017))
mod_evenness <- lm(evenness ~ year,dd_farm)
new_dd_farm$pred <- predict(mod_evenness, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),evenness=c(evenness_farmland_2001_low,evenness_farmland_2017_low))
mod_evenness <- lm(evenness ~ year,dd_farm)
new_dd_farm$lwr <- predict(mod_evenness, newdata=new_dd_farm)
dd_farm <- data.frame(year=c(2001,2017),evenness=c(evenness_farmland_2001_up,evenness_farmland_2017_up))
mod_evenness <- lm(evenness ~ year,dd_farm)
new_dd_farm$upr <- predict(mod_evenness, newdata=new_dd_farm)

dd_forest <- data.frame(year=c(2001,2017),evenness=c(evenness_forest_2001,evenness_forest_2017))
new_dd_forest <- data.frame(year=c(2001:2017))
mod_evenness <- lm(evenness ~ year,dd_forest)
new_dd_forest$pred <- predict(mod_evenness, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),evenness=c(evenness_forest_2001_low,evenness_forest_2017_low))
mod_evenness <- lm(evenness ~ year,dd_forest)
new_dd_forest$lwr <- predict(mod_evenness, newdata=new_dd_forest)
dd_forest <- data.frame(year=c(2001,2017),evenness=c(evenness_forest_2001_up,evenness_forest_2017_up))
mod_evenness <- lm(evenness ~ year,dd_forest)
new_dd_forest$upr <- predict(mod_evenness, newdata=new_dd_forest)

dd_natural_open <- data.frame(year=c(2001,2017),evenness=c(evenness_natural_openland_2001,evenness_natural_openland_2017))
new_dd_natural_open <- data.frame(year=c(2001:2017))
mod_evenness <- lm(evenness ~ year,dd_natural_open)
new_dd_natural_open$pred <- predict(mod_evenness, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),evenness=c(evenness_natural_openland_2001_low,evenness_natural_openland_2017_low))
mod_evenness <- lm(evenness ~ year,dd_natural_open)
new_dd_natural_open$lwr <- predict(mod_evenness, newdata=new_dd_natural_open)
dd_natural_open <- data.frame(year=c(2001,2017),evenness=c(evenness_natural_openland_2001_up,evenness_natural_openland_2017_up))
mod_evenness <- lm(evenness ~ year,dd_natural_open)
new_dd_natural_open$upr <- predict(mod_evenness, newdata=new_dd_natural_open)

dd_urban <- data.frame(year=c(2001,2017),evenness=c(evenness_urban_2001,evenness_urban_2017))
new_dd_urban <- data.frame(year=c(2001:2017))
mod_evenness <- lm(evenness ~ year,dd_urban)
new_dd_urban$pred <- predict(mod_evenness, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),evenness=c(evenness_urban_2001_low,evenness_urban_2017_low))
mod_evenness <- lm(evenness ~ year,dd_urban)
new_dd_urban$lwr <- predict(mod_evenness, newdata=new_dd_urban)
dd_urban <- data.frame(year=c(2001,2017),evenness=c(evenness_urban_2001_up,evenness_urban_2017_up))
mod_evenness <- lm(evenness ~ year,dd_urban)
new_dd_urban$upr <- predict(mod_evenness, newdata=new_dd_urban)

ggplot(new_dd,aes(x=year,y=pred)) + 
  geom_ribbon(aes(ymin=lwr, ymax=upr),fill = "black", alpha=0.2) +
  geom_line() +
  geom_ribbon(data=new_dd_farm,aes(ymin=lwr, ymax=upr),fill = "#f5b041", alpha=0.2) +
  geom_line(data=new_dd_farm, col="#f5b041") +
  geom_ribbon(data=new_dd_forest,aes(ymin=lwr, ymax=upr),fill = "#52be80", alpha=0.2) +
  geom_line(data=new_dd_forest, col="#52be80") +
  geom_ribbon(data=new_dd_natural_open,aes(ymin=lwr, ymax=upr),fill = "#2C65FD", alpha=0.2) +
  geom_line(data=new_dd_natural_open, col="#2C65FD") +
  geom_ribbon(data=new_dd_urban,aes(ymin=lwr, ymax=upr),fill = "#FF4001", alpha=0.2) +
  geom_line(data=new_dd_urban, col="#FF4001") +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") + ylab("Network eveness\ncontrolling for connectance") + 
  theme_modern()
  
  
ggplot(new_dd[which(new_dd$year %in% c(2001,2017)),],aes(x=as.factor(year),y=pred,group=year)) + 
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2, position=position_dodge(0.05), alpha=0.2) +
  geom_point(size=3,data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],col = "#f5b041") +
  geom_errorbar(data=new_dd_farm[which(new_dd_farm$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#f5b041", alpha=0.2) +
  geom_point(size=3,data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],col = "#52be80") +
  geom_errorbar(data=new_dd_forest[which(new_dd_forest$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#52be80", alpha=0.2) +
  geom_point(size=3,data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],col = "#2C65FD") +
  geom_errorbar(data=new_dd_natural_open[which(new_dd_natural_open$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#2C65FD", alpha=0.2) +
  geom_point(size=3,data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],col = "#FF4001") +
  geom_errorbar(data=new_dd_urban[which(new_dd_urban$year %in% c(2001,2017)),],aes(ymin=lwr, ymax=upr),width=.2, position=position_dodge(0.05),col = "#FF4001", alpha=0.2) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  xlab("Years") +  ylab("Network eveness\ncontrolling for connectance") + 
  theme_modern()

```

### Network example

```{r}

x <- droplevels(dataprp[which(dataprp$code_point==unique(dataprp$code_point)[10] & dataprp$year == 2008),])
# example for 2001
x <- droplevels(dataprp[which(dataprp$code_square==31782 & dataprp$year == 2002),])
#example for 2017
x <- droplevels(dataprp[which(dataprp$code_square==641279 & dataprp$year == 2016),])
# example no change in species but change in metrics
x <- droplevels(dataprp[which(dataprp$code_square==690715 & dataprp$year == 2016),])
# example no change in metric but change in species
x <- droplevels(dataprp[which(dataprp$code_square==520389 & dataprp$year == 2007),])



b <- levels(droplevels(x$code_sp))

datasso <- data_interaction_possible[which(data_interaction_possible$spA %in% b & data_interaction_possible$spB %in% b),c("spA","spB","obs_temp_asso")]

adj_mat <- datasso
for(i in 1:length(b)){
  adj_mat <- rbind(adj_mat, data.frame(spA=b[i],spB=b[i],obs_temp_asso=0))
}
adj_mat <- dcast(adj_mat, spA ~ spB, value.var="obs_temp_asso")
row.names(adj_mat) <- adj_mat$spA
adj_mat$spA <- NULL
adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]

library(network)
library(ggnet)

net <- network(as.matrix(adj_mat), directed = FALSE)
ggnet2(net, node.size = 6, node.color = "black", edge.size = 0.5, edge.color = "grey")

```

### Correlation between metrics

```{r}

cormat <- na.omit(network_structure_square[,c("nb_association","nb_sp","connectance","modularity","evenness")])

names(cormat) <- c("Number of associations", "Number of species","Connectance","Modularity","Evenness")

cormat <- round(cor(cormat),2)

# Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}
upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Pearson\nCorrelation") +
  theme_minimal()+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none") +
 coord_fixed()

```

## References

Clark, A. T., Ye, H., Isbell, F., Deyle, E. R., Cowles, J., Tilman, G. D., & Sugihara, G. (2015). Spatial convergent cross mapping to detect causal relationships from short time series. Ecology, 96(5), 1174-1181. https://doi.org/10.1890/14-1479.1

Rigal, S., Devictor, V., Gaüzère, P., Kéfi, S., Forsman, J. T., Kajanus, M. H., ... & Dakos, V. (2022). Biotic homogenisation in bird communities leads to large‐scale changes in species associations. Oikos, 2022(3), e08756. https://doi.org/10.1111/oik.08756
