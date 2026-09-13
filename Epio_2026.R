##
# load packages
##
pacman::p_load(ggcube,MASS,readxl,phangorn,fpc,ggnewscale,Morpho,ggtree,ape,
               patchwork,gridExtra,stringr,geomorph,scales,viridis,phylosignal,
               phylobase,mclust,factoextra,cluster,pvclust,tidyr,dplyr,gridExtra)


################################################################################
################################################################################
## --- Data organization --- ##
################################################################################
################################################################################

setwd('~/Box/Myprojects/Epioblasma/')

set.seed(117)

# load data
raw <- read_xlsx("data/raw/Procrustes-AllEpioRawData_Aug2026.xlsx",
                 sheet='Sheet1')

# Inspect raw data
raw <-as.data.frame(raw)
x<-c(raw[,9],raw[,11],raw[,13],raw[,15],
     raw[,17],raw[,19],raw[,21],raw[,23],raw[,25],
     raw[,27],raw[,29],raw[,31],raw[,33],raw[,35],
     raw[,37],raw[,39],raw[,41],raw[,43],raw[,45],
     raw[,47],raw[,49],raw[,51],raw[,53],raw[,55],
     raw[,57],raw[,59],raw[,61],raw[,63],raw[,65],
     raw[,67],raw[,69],raw[,71],raw[,73])
y<-c(raw[,10],raw[,12],raw[,14],
                 raw[,16],raw[,18],raw[,20],raw[,22],raw[,24],
                 raw[,26],raw[,28],raw[,30],raw[,32],raw[,34],
                 raw[,36],raw[,38],raw[,40],raw[,42],raw[,44],
                 raw[,46],raw[,48],raw[,50],raw[,52],raw[,54],
                 raw[,56],raw[,58],raw[,60],raw[,62],raw[,64],
                 raw[,66],raw[,68],raw[,70],raw[,72],raw[,74])

lot <- rep(raw$ID, times=33)
raw_check <- cbind.data.frame(lot,x,y)

#LMs are have been transformed by GPA
ggplot(raw_check,aes(x=x,y=y,color=as.factor(lot)))+
  geom_point()+
  theme(legend.position = 'off')

# Correct species names in raw file
raw$Species <- gsub("ahlstedi"  , "ahlstedti" , raw$Species , fixed = TRUE)

# load environmental data
stream_df <- read_xlsx("data/raw/StreamData_byLot.xlsx",
                 sheet='Sheet1')
stream_df <- stream_df %>% filter(!QEMA == '?')
stream_df$ln_flow <- log(as.numeric(stream_df$QEMA)*0.0283168)
# stream_df <- stream_df |> filter(!is.na(QEMA)) |> mutate(QEMA_ms = as.numeric(QEMA)*0.0283168)
# range(stream_df$QEMA_ms)
quantile(stream_df$QEMA_ms, prob=c(0.25, 0.5, 0.75))
IQR(stream_df$QEMA_ms,na.rm = TRUE)
# median(stream_df$QEMA_ms)
stream_df <- stream_df %>% dplyr::select(`Lot Unlettter`,ln_flow)
stream_df <- stream_df %>% dplyr::distinct(`Lot Unlettter`,ln_flow, .keep_all = TRUE)
stream_df <- stream_df %>% rename(Lot_no = `Lot Unlettter`)
# Join ln flow rate data to raw df by lot number
raw <- raw %>% rename(Lot_no = `Lot Number`)
raw <- raw %>% left_join(as.data.frame(stream_df), by = 'Lot_no')

## Collate each dataset
# all specimens with environmental data
raw_env <- raw %>% filter(!is.na(ln_flow))
# all specimens with size data (all rows also have ln Flow data)
raw_sz <- raw %>% filter(!is.na(area_mm2),Species %in% c('ahlstedti','brevidens','flexuosa'))
# remove species that do not have at least 5 specimens for each sex
species_list <- raw |> group_by(Species,Sex) |> count() |> 
          filter(n <= 5) |> ungroup() |> select(Species)
raw_all_lowsamps <- raw |> filter(!Species %in% species_list$Species)



# Organize landmarks
n_landmarks <- 33
n_dim <- 2
n_specimens_all <- nrow(raw) #1382
length(unique(raw$Lot_no)) # 242 lots
nrow(raw |> filter(Sex == 'M')) #813
nrow(raw |> filter(Sex == 'F')) #565
n_specimens_env <- nrow(raw_env) #1,168
n_specimens_sz <- nrow(raw_sz)
table(raw_sz$Species)
n_specimens_all_lowsamps <- nrow(raw_all_lowsamps)

# Extract landmark columns (X1-Y33)
coord_cols <- grep("^X|^Y", colnames(raw))
coords_mat_all <- as.matrix(raw[, coord_cols])
dim(coords_mat_all)
# Add rownames for all-specimens PCA
rownames(coords_mat_all) <-paste0(raw$Species,'.',raw$Sex,'.',raw$ID)
coords_mat_env <- as.matrix(raw_env[, coord_cols])
coords_mat_sz <- as.matrix(raw_sz[, coord_cols])
coords_mall_lowsamps <- as.matrix(raw_all_lowsamps[, coord_cols])

# Convert to 3D array: landmarks × dimensions × specimens
coords_all <- array(NA, dim = c(n_landmarks, n_dim, n_specimens_all))
for(i in 1:n_specimens_all){
  coords_all[,,i] <- matrix(coords_mat_all[i,], ncol = n_dim, byrow = TRUE)
}
coords_env <- array(NA, dim = c(n_landmarks, n_dim, n_specimens_env))
for(i in 1:n_specimens_env){
  coords_env[,,i] <- matrix(coords_mat_env[i,], ncol = n_dim, byrow = TRUE)
}
coords_sz <- array(NA, dim = c(n_landmarks, n_dim, n_specimens_sz))
for(i in 1:n_specimens_sz){
  coords_sz[,,i] <- matrix(coords_mat_sz[i,], ncol = n_dim, byrow = TRUE)
}
coords_all_lowsamps <- array(NA, dim = c(n_landmarks, n_dim, n_specimens_all_lowsamps))
for(i in 1:n_specimens_all_lowsamps){
  coords_all_lowsamps[,,i] <- matrix(coords_mall_lowsamps[i,], ncol = n_dim, byrow = TRUE)
}

## Create Geomorph dataframes for all subsets of all-specimens data

# Geomorph data frame - ALL SPECIMENS
Epio.df_all <- geomorph.data.frame(
  coords = coords_all,
  species = raw$Species,
  subgenus = raw$Subgenus,
  sex = raw$Sex)

# Geomorph data frame - ENVIRONMENT
Epio.df_env <- geomorph.data.frame(
  coords = coords_env,
  species = raw_env$Species,
  subgenus = raw_env$Subgenus,
  sex = raw_env$Sex,
  ln_flow = raw_env$ln_flow)

# Geomorph data frame - SIZE
Epio.df_sz <- geomorph.data.frame(
  coords = coords_sz,
  species = raw_sz$Species,
  subgenus = raw_sz$Subgenus,
  sex = raw_sz$Sex,
  ln_flow = raw_sz$ln_flow,
  ln_size = log(raw_sz$area_mm2))

# Geomorph data frame - ALL SPECIMENS - LOW - for trajectory analyses
Epio.df_all_s <- geomorph.data.frame(
  coords = coords_all_lowsamps,
  species = raw_all_lowsamps$Species,
  subgenus = raw_all_lowsamps$Subgenus,
  sex = raw_all_lowsamps$Sex)

## Compute Mean Shape for each sex of each species from all-specimens df
### Compute mean shapes ###
gp <- interaction(Epio.df_all$species, Epio.df_all$sex)
levels(gp) <- gsub(" ", "", levels(gp))
group_levels <- levels(gp)

mean_shapes <- array(NA, dim = c(n_landmarks, n_dim, length(group_levels)),
                     dimnames = list(NULL, NULL, group_levels))

for(i in seq_along(group_levels)){
  specimens <- which(gp == group_levels[i])
  mean_shapes[,,i] <- mshape(coords_all[,,specimens])
}
# Flatten mean shapes into a matrix
mean_shapes_mat <- t(apply(mean_shapes, 3, function(x) as.vector(x)))
species_sex <- rownames(mean_shapes_mat)

species_mean <- gsub("\\..*", "", species_sex)
sex_mean <- sub(".*\\.","",species_sex)

# Mean shape geomorph dataframe
Epio.df.mean <- geomorph.data.frame(
  coords = mean_shapes,
  species = species_mean,
  sex = sex_mean)

# Calculate mean Epioblasma shape from species*sex mean specimens
ref <- mshape(Epio.df.mean$coords)

# # Create links to define outline of specimens
# links <- define.links(ref, ptsize = 1, links = NULL)
links <- array(dim = c(30,2))
links[,1] <- c(1:29,1) 
links[,2] <- c(2:30,30)
Epio.df.mean$links <- links
# # Save mean data geomorph dataframe so links don't have to be defined each R session
# saveRDS(Epio.df.mean, file = 'data/raw/Epio.df.mean.Rdata')

# # Load mean Epio geomorph dataframe
# Epio.df.mean <- readRDS(file = 'data/raw/Epio.df.mean.Rdata')



# Check how many lots of Epioblasma exists at NMNH, OSUM, and UMMZ relative to other colls
# pfeiff_epio <- read_xlsx('~/Box/Myprojects/river_map/MusselMuseum/Pfeiffer2024.xlsx',
#                         sheet = '3c_all_records_converted')
# pfeiff_epio <- pfeiff_epio %>% filter(genus == 'Epioblasma')
# pfeiff_table <- pfeiff_epio %>% dplyr::count(institutionCode, sort = T)
# 
# sum(pfeiff_table %>% filter(institutionCode %in% c('OSUM','UMMZ', 'USNM')) %>% 
#       dplyr::select(n))/sum(pfeiff_table$n)
# 0.370666

# Check sample sizes of complete (raw) dataset for Table 1

no_males <- raw %>% filter(Sex == 'M') %>% dplyr::select(Species, Subgenus, Sex)
subgenus <- no_males %>% dplyr::select(Species, Subgenus)
subgenus <- distinct(subgenus, Species, Subgenus)
no_males <- no_males %>% dplyr::count(Species, name = 'males', sort = T)
no_females <- raw %>% filter(Sex == 'F') %>% dplyr::select(Species, Sex)
no_females <- no_females %>% dplyr::count(Species, name = 'females', sort = T)

table1 <- subgenus %>% left_join(no_females, by='Species') %>% 
  left_join(no_males,by='Species')
table1 <- table1 %>% mutate(print_species = paste0('E. ',Species)) |> 
  arrange(Subgenus)
# write.csv(table1,'figures/tables/table1.csv')



################################################################################
################################################################################
## --- AIM 1: Testing for Sexual Dimorphism --- ##
################################################################################
################################################################################


################################################################################
## --- AIM 1A: Morphological Disparity --- ##
################################################################################


### Procrustes ANOVA / linear modeling ###
# Does shell shape differ between species, sexes, and their interaction?

# Size/Allometry
# Fit_sz <- procD.lm(coords ~ species * sex * ln_size, data = Epio.df_sz, iter = 9999)
# saveRDS(Fit_sz, file = 'data/raw/MANOVA/Fit_sz.Rdata')
# Load MANOVA
Fit_sz <- readRDS(file = 'data/raw/MANOVA/Fit_sz.Rdata')
summary(Fit_sz)
# Df      SS       MS     Rsq        F       Z Pr(>F)    
# species               2 0.48854 0.244269 0.39727 116.2540 10.0135  1e-04 ***
#   sex                   1 0.21294 0.212942 0.17316 101.3448  6.3728  1e-04 ***
#   ln_size               1 0.04408 0.044081 0.03585  20.9793  5.7465  1e-04 ***
#   species:sex           2 0.13018 0.065092 0.10586  30.9792  8.6278  1e-04 ***
#   species:ln_size       2 0.01813 0.009065 0.01474   4.3141  4.3663  1e-04 ***
#   sex:ln_size           1 0.00894 0.008945 0.00727   4.2570  3.1427  4e-04 ***
#   species:sex:ln_size   2 0.03275 0.016373 0.02663   7.7925  5.9568  1e-04 ***
#   Residuals           140 0.29416 0.002101 0.23921                            
# Total               151 1.22973     


# Stream size effects
# Fit_env <- procD.lm(coords ~ species * sex + ln_flow, data = Epio.df_env, iter = 9999)
# saveRDS(Fit_env, file = 'data/raw/MANOVA/Fit_env.Rdata')
# Load MANOVA
Fit_env <- readRDS(file = 'data/raw/MANOVA/Fit_env.Rdata')
summary(Fit_env)
# Df     SS      MS     Rsq        F      Z Pr(>F)    
# species       26 3.5680 0.13723 0.40098  44.2608 33.968 0.0001 ***
#   sex            1 1.1378 1.13783 0.12787 366.9781 10.285 0.0001 ***
#   ln_flow        1 0.0041 0.00408 0.00046   1.3145  0.773 0.2213    
# species:sex   26 0.7375 0.02836 0.08288   9.1483 22.011 0.0001 ***
#   Residuals   1113 3.4509 0.00310 0.38781                           
# Total       1167 8.8983 


# Sex * Species effects - ALL SPECIMENS
# Fit_all <- procD.lm(coords ~ species * sex, data = Epio.df_all, iter = 9999)
# saveRDS(Fit_all, file = 'data/raw/MANOVA/Fit_all.Rdata')
# Load MANOVA
Fit_all <- readRDS(file = 'data/raw/MANOVA/Fit_all.Rdata')
summary(Fit_all)
# Df      SS      MS     Rsq        F      Z    Pr(>F)    
# species       27  4.2311 0.15671 0.39845  48.6772 35.475 9.999e-05 ***
#   sex            1  1.2529 1.25290 0.11799 389.1792  9.985 9.999e-05 ***
#   species:sex   27  0.8662 0.03208 0.08157   9.9651 18.654 9.999e-05 ***
#   Residuals   1326  4.2688 0.00322 0.40200                              
# Total       1381 10.6190  


# Sex * Species effects - SPP W/ LOW N OMITTED
# Fit_all_s <- procD.lm(coords ~ species * sex, data = Epio.df_all_s, iter = 9999)
# saveRDS(Fit_all_s, file = 'data/raw/MANOVA/Fit_all_s.Rdata')
# Load MANOVA
Fit_all_s <- readRDS(file = 'data/raw/MANOVA/Fit_all_s.Rdata')
summary(Fit_all_s)
# Df      SS      MS     Rsq       F      Z Pr(>F)    
# species       23  4.1577 0.18077 0.39961  56.823 32.103  1e-04 ***
#   sex            1  1.2301 1.23007 0.11823 386.660 10.392  1e-04 ***
#   species:sex   23  0.8459 0.03678 0.08130  11.561 23.370  1e-04 ***
#   Residuals   1311  4.1706 0.00318 0.40086                          
# Total       1358 10.4043 


####
### PCA on all specimens ###
####
# all_shapes_mat <- t(apply(Epio.df_all$coords, 3, function(x) as.vector(x)))
pca_all <- prcomp(coords_mat_all, scale. = FALSE)
pc_scores_all <- data.frame(pca_all$x)
pc_scores_all$Group <- rownames(pc_scores_all)
pc_scores_all <- pc_scores_all %>% mutate(Species = sapply(strsplit(as.character(pc_scores_all$Group), "\\."), `[`, 1))
pc_scores_all <- pc_scores_all %>% mutate(Sex = sapply(strsplit(as.character(pc_scores_all$Group), "\\."), `[`, 2))
var_explained_all <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)
var_all <- 100 * pca_all$sdev^2 / sum(pca_all$sdev^2)
var_all <- data.frame(var_all = var_all[1:8])
var_all$PC <- c(1:8)

p_allpca_12 <- ggplot() +
  geom_point(data = pc_scores_all, aes(x = PC1, y = PC2,
                                       shape = Sex, color = Species, fill = Species),
             size = 0.8, alpha = 0.3, stroke = 0.9) +
  scale_shape_manual(labels = c('Female','Male'), values=c(0,19)) +
  scale_color_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  scale_fill_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  theme_linedraw() +
  ylim(min(pc_scores_all$PC2,pc_scores_all$PC3),max(pc_scores_all$PC2,pc_scores_all$PC3)) +
  xlim(min(pc_scores_all$PC1,pc_scores_all$PC3),max(pc_scores_all$PC1,pc_scores_all$PC3)) +
  coord_fixed(ratio=1) +
  labs(x = paste0("PC1 (", var_explained_all[1], "%)"),
       y = paste0("PC2 (", var_explained_all[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10, face = 'italic')) +
  guides(color = guide_legend(ncol = 3),
         fill = guide_legend(ncol = 3))

p_allpca_23 <- ggplot() +
  geom_point(data = pc_scores_all, aes(x = PC3, y = PC2,
                                       shape = Sex, color = Species, fill = Species),
             size = 0.8, alpha = 0.3, stroke = 0.9) +
  scale_shape_manual(labels = c('Female','Male'), values=c(0,19)) +
  scale_color_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  scale_fill_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  theme_linedraw() +
  ylim(min(pc_scores_all$PC2,pc_scores_all$PC3),max(pc_scores_all$PC2,pc_scores_all$PC3)) +
  xlim(min(pc_scores_all$PC1,pc_scores_all$PC3),max(pc_scores_all$PC1,pc_scores_all$PC3)) +
  coord_fixed(ratio=1) +
  labs(x = paste0("PC3 (", var_explained_all[3], "%)"),
       y = paste0("PC2 (", var_explained_all[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10, face = 'italic')) +
  guides(color = guide_legend(ncol = 3),
         fill = guide_legend(ncol = 3))

p_allpca_13 <- ggplot() +
  geom_point(data = pc_scores_all, aes(x = PC1, y = PC3,
                                       shape = Sex, color = Species, fill = Species),
             size = 0.8, alpha = 0.3, stroke = 0.9) +
  scale_shape_manual(labels = c('Female','Male'), values=c(0,19)) +
  scale_color_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  scale_fill_discrete(labels = unique(paste('E.',pc_scores_all$Species))) +
  theme_linedraw() +
  ylim(min(pc_scores_all$PC2,pc_scores_all$PC3),max(pc_scores_all$PC2,pc_scores_all$PC3)) +
  xlim(min(pc_scores_all$PC1,pc_scores_all$PC3),max(pc_scores_all$PC1,pc_scores_all$PC3)) +
  coord_fixed(ratio=1) +
  labs(x = paste0("PC1 (", var_explained_all[1], "%)"),
       y = paste0("PC3 (", var_explained_all[3], "%)")) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10, face = 'italic')) +
  guides(color = guide_legend(ncol = 3),
         fill = guide_legend(ncol = 3))

# plot determine PC variation
bar_var <- ggplot(data=var_all) +
  geom_bar(aes(x=PC, y=var_all),stat = "identity") +
  labs(y='Variance explained (%)',x='Principal component', 
       title = 'Variance per PC') +
  theme_linedraw() +
  theme(legend.position = 'off',
        plot.title = element_text(hjust = 0.5,size=11),
        axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        axis.text = element_text(size=9)) 

p_allpca <- p_allpca_12 + p_allpca_23 + p_allpca_13 + guide_area() + 
  plot_layout(ncol=2, guides = "collect")

p_allpca <- p_allpca + inset_element(bar_var, left = 0.6,bottom = -0.25,right = 1,top=0.25) +
  plot_annotation(tag_levels = 'A') + theme(plot.tag.position = c(0, 1),
                                            plot.tag = element_text(size = 8, hjust = 0, vjust = 0))

# ggsave('figures/pca_all.png', p_allpca, width = 12, height = 10, dpi = 600, bg = "white")


###Plot disparity landscape for all specimens of each sex

####
### Topography of sex occupation within morphospace
####
x.range <- range(pc_scores_all$PC1)
y.range <- range(pc_scores_all$PC2)
grid.n <- 175
grid <- expand.grid(x = seq(x.range[1], x.range[2], length.out = grid.n),
                    y = seq(y.range[1], y.range[2], length.out = grid.n))
pc.diag <- sqrt(diff(x.range)^2 + diff(y.range)^2)

# radius size
# was 0.05
# 7 % of diagonal of morphospace
radius <- pc.diag * 0.07

# Surface function, includes quartic kernel for smoothing
local_occupancy <- function(points, query, radius){
  R2 <- radius^2
  z <- numeric(nrow(query))
  for(i in seq_len(nrow(query))){
    dx <- points$PC1 - query$PC1[i]
    dy <- points$PC2 - query$PC2[i]
    d2 <- dx^2 + dy^2
    inside <- d2 < R2
    if(any(inside)){
      u2 <- d2[inside] / R2
      weights <- (1 - u2)^2
      z[i] <- sum(weights)
    }
  }
  z
}

local_surface <- function(points, grid, radius){
  grid$z <- local_occupancy(points = points,query = data.frame(PC1 = grid$x,
      PC2 = grid$y),radius = radius)
  grid
}

# filter pc_scores by sex
pc_scores_all_f <- pc_scores_all |> filter(Sex == 'F')
pc_scores_all_m <- pc_scores_all |> filter(Sex == 'M')

# Calculate male and female occupancy surfaces
male.surface <- local_surface(pc_scores_all_m, grid, radius)
female.surface <- local_surface(pc_scores_all_f, grid, radius)
z.max <- max(male.surface$z,female.surface$z)


pc_scores_all_f$z <- local_occupancy(points = pc_scores_all_f,
                                     query = pc_scores_all_f,radius = radius)
pc_scores_all_m$z <- local_occupancy(points = pc_scores_all_m,
                                     query = pc_scores_all_m,radius = radius)

female_surface_pc_scores <- data.frame(x = pc_scores_all_f$PC1,
                                       y= pc_scores_all_f$PC2,
                                       z = pc_scores_all_f$z)
male_surface_pc_scores <- data.frame(x = pc_scores_all_m$PC1,
                                     y= pc_scores_all_m$PC2,
                                     z = pc_scores_all_m$z)

lim_x <-  max(pc_scores_all$PC1) - min(pc_scores_all$PC1)
lim_y <- max(pc_scores_all$PC2) - min(pc_scores_all$PC2)


#plot surfaces
f_disparity <- ggplot() +
  geom_surface_3d(data=female.surface, aes(x=x, y=y, z=z, fill = z, color = z)) +
  scale_fill_gradient(low = "white",high = "tomato", limits = c(0, z.max),
                      breaks = c(0, 40, 80)) + 
  scale_color_gradient(low = "white",high = "tomato", limits = c(0, z.max),
                       breaks = c(0, 40, 80)) +
  geom_point(data=female_surface_pc_scores, aes(x=x, y=y, z=z),
             shape=21,size=0.55,color='black', fill='white',alpha=0.7) +
  scale_z_continuous(limits = c(0, z.max)) +
  coord_3d(ratio = c(lim_x*5, lim_y*5, 1), expand = FALSE, panels = "zmin",
           light = light(method = "diffuse", mode = "hsv", 
                         direction = c(0, 1, 0), contrast = 0.5)) +
  guides(fill = guide_colorbar_3d()) +
  theme_light() +
  labs(x = paste0("PC1 (", var_explained_all[1], "%)"),
       y = paste0("PC2 (", var_explained_all[2], "%)"),
       fill = 'No. of\nspecimens',
       color = 'No. of\nspecimens') +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = c(0.90, 0.65))

m_disparity <- ggplot() +
  geom_surface_3d(data=male.surface, aes(x=x, y=y, z=z, fill = z, color = z)) +
  scale_fill_gradient(low = "white",high = "steelblue", limits = c(0, z.max),
                      breaks = c(0, 40, 80)) + 
  scale_color_gradient(low = "white",high = "steelblue", limits = c(0, z.max),
                       breaks = c(0, 40, 80)) +
  geom_point(data=male_surface_pc_scores, aes(x=x, y=y, z=z),
             shape=21,size=0.55,color='black', fill='white',alpha=0.7) +
  scale_z_continuous(limits = c(0, z.max)) +
  coord_3d(ratio = c(lim_x*5, lim_y*5, 1), expand = FALSE, panels = "zmin",
           light = light(method = "diffuse", mode = "hsv", 
                         direction = c(0, 1, 0), contrast = 0.5)) +
  guides(fill = guide_colorbar_3d()) +
  theme_light() +
  labs(x = paste0("PC1 (", var_explained_all[1], "%)"),
       y = paste0("PC2 (", var_explained_all[2], "%)"),
       fill = 'No. of\nspecimens',
       color = 'No. of\nspecimens') +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=12),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = c(0.90, 0.65))

p_disparity <- f_disparity/m_disparity

# ggsave('figures/disparity/disparity_surfaces.png', p_disparity, width = 6, height = 10, dpi = 600, bg = "white")


# Morphological disparity - WITHIN SPP


# Procrustes variance function
proc_var <- function(shape_mat) {
  n <- nrow(shape_mat)
  if (n < 2) return(NA)
  centered <- scale(shape_mat, center = TRUE, scale = FALSE) # deviations from group mean shape
  sum(centered^2) / n   # trace(SumofSqs)/n = trace(pop cov matrix)
}

# Bootstrap: resample specimens with replacement, recompute disparity
boot_disparity <- function(shape_mat, n_boot = 2000, conf = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(shape_mat)
  obs <- proc_var(shape_mat)
  boots <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    proc_var(shape_mat[idx, , drop = FALSE])
  })
  alpha <- (1 - conf) / 2
  ci <- quantile(boots, probs = c(alpha, 1 - alpha), na.rm = TRUE)
  list(obs = obs, boot_vals = boots, n = n,
       lower = unname(ci[1]), upper = unname(ci[2]))
}

# Compute Procrustes variance estimate + bootstrap 95% CI
species_list <- levels(factor(Epio.df_all_s$species))
sexes <- levels(factor(Epio.df_all_s$sex))
min_n <- 5
n_boot <- 2000

disparity_table <- data.frame()

for (sp in species_list) {
  for (sx in sexes) {
    keep <- which(Epio.df_all_s$species == sp & Epio.df_all_s$sex == sx)
    if (length(keep) < min_n) next
    shape_mat <- two.d.array(Epio.df_all_s$coords[,, keep, drop = FALSE])
    bd <- boot_disparity(shape_mat, n_boot = n_boot, seed = 1)
    disparity_table <- rbind(disparity_table,
                             data.frame(species = sp, sex = sx, n = bd$n,
                                        disparity = bd$obs, lower = bd$lower, upper = bd$upper))
  }
}

# Procrustes variance should be identical to morphol.disparity
# morphol.disparity(coords ~ species * sex, groups = ~ species * sex, data = Epio.df_all_s, partial = FALSE, iter = 9999)
# Procrustes variances for defined groups
# ahlstedti.F      ahlstedti.M    arcaeformis.F    arcaeformis.M   biemarginata.F   biemarginata.M      brevidens.F 
# 0.004588880      0.002065677      0.004022008      0.002888892      0.002845275      0.002517012      0.002588757 
# brevidens.M   capsaeformis.F   capsaeformis.M       flexuosa.F       flexuosa.M     florentina.F     florentina.M 
# 0.002440148      0.004072219      0.002682580      0.003275886      0.005311100      0.002807849      0.002738183 
# gubernaculum.F   gubernaculum.M       haysiana.F       haysiana.M         lenior.F         lenior.M    metastriata.F 
# 0.002602224      0.001810229      0.002442801      0.003102362      0.001717045      0.001499363      0.002689643 
# metastriata.M      obliquata.F      obliquata.M othcaloogensis.F othcaloogensis.M         penita.F         penita.M 
# 0.002851164      0.002078687      0.002346836      0.002112382      0.002057228      0.003743583      0.002900229 
# perobliqua.F     perobliqua.M      personata.F      personata.M      propinqua.F      propinqua.M       rangiana.F 
# 0.001684752      0.002879755      0.002527748      0.002844089      0.002374588      0.002802421      0.004174029 
# rangiana.M      sampsonii.F      sampsonii.M   stewardsonii.F   stewardsonii.M       torulosa.F       torulosa.M 
# 0.002429569      0.004312252      0.003094311      0.004997448      0.001795762      0.003890421      0.003028666 
# triquetra.F      triquetra.M      turgidula.F      turgidula.M        walkeri.F        walkeri.M 
# 0.004322706      0.003404852      0.003653104      0.002842571      0.002381096      0.002264065 
disparity_table |> select(species,sex,disparity) # looks good

# Reformat table
disparity_t_wide <- disparity_table |> pivot_wider(names_from = sex, values_from = c(n, disparity, lower, upper))
disparity_t_wide <- as.data.frame(disparity_t_wide)

# which species have non-overlapping CI
rng <-  cbind(pmin(disparity_t_wide[,6], disparity_t_wide[,8]), pmax(disparity_t_wide[,6], disparity_t_wide[,8]),
              pmin(disparity_t_wide[,7], disparity_t_wide[,9]), pmax(disparity_t_wide[,7], disparity_t_wide[,9]))
disparity_t_wide$overlap <-  (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])

# scale x axis
x_range <- range(c(disparity_table$lower, disparity_table$upper), na.rm = TRUE)
# add overlap column to disparity table df
disparity_table <- disparity_table |> left_join(disparity_t_wide |> select(species,overlap),by='species')

mean_male_disp <- disparity_table |> filter(sex=='M')
mean_male_disp <- mean(mean_male_disp$disparity)
mean_female_disp <- disparity_table |> filter(sex=='F')
mean_female_disp <- mean(mean_female_disp$disparity)

# plot
intra_disparity <- ggplot(disparity_table, aes(x = disparity, y = species)) +
  geom_pointrange(aes(xmin = lower, xmax = upper,color=overlap), size = 0.4) +
  scale_y_discrete(limits = rev, labels = paste0('E. ', rev(unique(disparity_table$species)))) +
  scale_color_manual(values = c('red','black'),guide = 'none') +
  facet_wrap(~ sex, nrow = 1,
             labeller = as_labeller(c('M' = "Male", 'F' = "Female"))) + 
  geom_vline(data=filter(disparity_table, sex=="M"),
    aes(xintercept = mean_male_disp), colour = "black", 
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(data=filter(disparity_table, sex=="F"),
             aes(xintercept = mean_female_disp), colour = "black", 
             linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(limits = x_range) +
  labs(title = "Intraspecific shell shape disparity",
    x = "Procrustes variance (disparity)", y = "Species") +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=16),
        strip.background = element_rect(fill = "gray99"),
        strip.text = element_text(size=14,color='black'),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=11,face = 'italic'),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=16))

# ggsave('figures/disparity/intra_disparity.png', intra_disparity, width = 10, height = 6, dpi = 600, bg = "white")



# Check to see if disparity is influenced by size
morphol.disparity(coords ~ species * sex * ln_size, groups = ~ species * sex, data = Epio.df_sz, partial = FALSE,
                  iter = 9999)
# Procrustes variances for defined groups
# ahlstedti.F ahlstedti.M brevidens.F brevidens.M  flexuosa.F  flexuosa.M 
# 0.002468188 0.001969745 0.001977984 0.001652691 0.001283358 0.002044183  
morphol.disparity(coords ~ species * sex, groups = ~ species * sex, data = Epio.df_sz, partial = FALSE,
                  iter = 9999)
# Procrustes variances for defined groups
# ahlstedti.F ahlstedti.M brevidens.F brevidens.M  flexuosa.F  flexuosa.M 
# 0.003796457 0.002052038 0.002071779 0.002360958 0.003094152 0.003029872 

# size adjusted residuals from species * sex * ln_size MANOVA
resid_sz <- Fit_sz$LM$residuals

species_list <- levels(factor(raw_sz$Species))
sexes        <- levels(factor(raw_sz$Sex))
min_n  <- 10
n_boot <- 2000

intra_disparity_sz <- data.frame()

for (sp in species_list) {
  for (sx in sexes) {
    keep <- which(raw_sz$Species == sp & raw_sz$Sex == sx)
    if (length(keep) < min_n) next
    
    shape_raw <- two.d.array(Epio.df_sz$coords[,, keep, drop = FALSE])
    shape_adj <- resid_sz[keep, , drop = FALSE]
    
    bd_raw <- boot_disparity(shape_raw, n_boot = n_boot, seed = 1)
    bd_adj <- boot_disparity(shape_adj, n_boot = n_boot, seed = 1)
    
    intra_disparity_sz <- rbind(intra_disparity_sz,
                     data.frame(species = sp, sex = sx, n = length(keep), type = "Raw",
                                disparity = bd_raw$obs, lower = bd_raw$lower, upper = bd_raw$upper),
                     data.frame(species = sp, sex = sx, n = length(keep), type = "Corrected",
                                disparity = bd_adj$obs, lower = bd_adj$lower, upper = bd_adj$upper))
  }
}


intra_disparity_sz_p <- ggplot(intra_disparity_sz, aes(x = disparity, y = species)) +
  geom_pointrange(aes(xmin = lower, xmax = upper,color=type), size = 0.6,
                  position = position_jitter(height = 0.1)) +
  scale_y_discrete(limits = rev, labels = paste0('E. ', rev(unique(intra_disparity_sz$species)))) +
  scale_color_manual(values = c('red','black'), labels = c('Size corrected', 'Raw')) +
  facet_wrap(~ sex, nrow = 1,
             labeller = as_labeller(c('M' = "Male", 'F' = "Female"))) + 
  labs(title = "Intraspecific shell shape disparity",
       x = "Procrustes variance (disparity)", y = "Species",
       color = 'Disparity \nmeasurement') +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size=16),
        strip.background = element_rect(fill = "gray99"),
        strip.text = element_text(size=14,color='black'),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=11,face = 'italic'),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=16))

# ggsave('figures/disparity/intra_disparity_sz_p.png', intra_disparity_sz_p, width = 10, height = 6, dpi = 600, bg = "white")




# Morphological disparity - ACROSS SPP
genus_disp <- morphol.disparity(coords ~ sex, groups = ~ sex, data = Epio.df.mean, partial = FALSE,
                                iter = 9999)
genus_disp$Procrustes.var
# F           M 
# 0.004974355 0.003180057 
genus_disp$PV.dist.Pval
# F      M
# F 1.0000 0.0257
# M 0.0257 1.0000

################################################################################
## --- AIM 1B: Magnitude & Direction of Sexual Dimorphism --- ##
################################################################################

## Procrustes distances - Magnitude of sexual dimorph / Hypothesis testing

spsex_interaction <- interaction(Epio.df_all_s$species, Epio.df_all_s$sex)
PW_dist <- pairwise(Fit_all_s, groups = spsex_interaction, covariate = NULL)

# Pairwise distances between means
options(max.print = 8000)
pw_summ_out <- capture.output(summary(PW_dist, confidence = 0.95, test.type = "dist"),
                         file=NULL,append=T)

pw_summ_mat <- do.call(rbind, strsplit(trimws(pw_summ_out[13:1140]), "\\s+"))
pw_summ <- data.frame(
  species.a = pw_summ_mat[,1],
  d = as.numeric(pw_summ_mat[,2]),
  UCL_95 = as.numeric(pw_summ_mat[,3]),
  Z = as.numeric(pw_summ_mat[,4]),
  pvalue = as.numeric(pw_summ_mat[,5]),
  stringsAsFactors = FALSE)

# Get Pairwise data better formatted to only get sexual dimorph between same sp
pw_summ <- separate(pw_summ, species.a, c('species.a', "species.b"),sep = ':')
pw_summ <- separate_wider_delim(pw_summ, cols = species.a, delim=".", names = c('species.a', "sex.a"))
pw_summ <- separate_wider_delim(pw_summ, cols = species.b, delim=".", names = c('species.b', "sex.b"))

### Procrustes distances between sexes for each species
procrust_d <- pw_summ %>% filter(species.a == species.b)
procrust_d <- procrust_d %>% mutate(proc_d_sig = if_else(pvalue > 0.05, "notsig", "sig"))
# Add omitted species to get them to be "nonsig"
omitted_spp <- table1 |> anti_join(procrust_d, by = c('Species' = 'species.a'))
omitted_spp <- data.frame(
  species.a = omitted_spp$Species,
  sex.a = rep('F',nrow(omitted_spp)),
  species.b = omitted_spp$Species,
  sex.b = rep('M',nrow(omitted_spp)),
  d = rep(NA,nrow(omitted_spp)),
  UCL_95 = rep(NA,nrow(omitted_spp)),
  Z = rep(NA,nrow(omitted_spp)),
  pvalue = rep(NA,nrow(omitted_spp)),
  proc_d_sig = rep("notsig",nrow(omitted_spp)),
  stringsAsFactors = FALSE)
procrust_d <- rbind.data.frame(procrust_d,omitted_spp)
# View(procrust_d)
# write.csv(procrust_d,'figures/procrustes_dist/proc_dist_results.csv')

##
# Get Pairwise Slope Vectors -> Direction of sexual dimorphism vectors
##
traj_ana <- trajectory.analysis(Fit_all_s, groups = Epio.df_all_s$species, 
                          traj.pts = Epio.df_all_s$sex, print.progress = FALSE)
traj_ana$TC
summary(traj_ana, attribute = "MD") # magnitude 
# Check to see above if Procrustes distance match between traj_ana & procrust_d
# ahlstedti    arcaeformis   biemarginata      brevidens   capsaeformis       flexuosa     florentina   gubernaculum 
# 0.13322239     0.05389930     0.08513379     0.05149672     0.08886272     0.08739962     0.09940638     0.13046359 
# haysiana         lenior    metastriata      obliquata othcaloogensis         penita     perobliqua      personata 
# 0.03582515     0.06496220     0.04565258     0.09055491     0.05738516     0.05702073     0.07083081     0.05047891 
# propinqua       rangiana      sampsonii   stewardsonii       torulosa      triquetra      turgidula        walkeri 
# 0.05268189     0.08970722     0.04895832     0.05839677     0.09900996     0.05761347     0.09362319     0.08309580 
procrust_d |> select(species.a,d) # looks good

#####
# VECTOR CORRELATIONS -> DIRECTION OF PHENOTYPIC CHANGE
#####

# summary(traj_ana, attribute = "TC") 
options(max.print = 8000)
vectors_out <- capture.output(summary(traj_ana, attribute = "TC"),
                              file=NULL,append=T)

vectors_out_mat <- do.call(rbind, strsplit(trimws(vectors_out[13:288]), "\\s+"))
vectors_out_df <- data.frame(
  species.comp = vectors_out_mat[,1],
  vector_cor = as.numeric(vectors_out_mat[,2]),
  angle = as.numeric(vectors_out_mat[,3]),
  UCL_95 = as.numeric(vectors_out_mat[,4]),
  Z = as.numeric(vectors_out_mat[,5]),
  pvalue = as.numeric(vectors_out_mat[,6]),
  stringsAsFactors = FALSE)


p_vectors <- ggplot(data=vectors_out_df,aes(vector_cor)) +
  # geom_histogram(binwidth = 0.05,fill = 'lightblue') +
  # geom_freqpoly(fill='gray') +
  geom_density(colour = "black",fill = "lightblue", linewidth = 0.5) +
  geom_vline(xintercept = 0, colour = "red", 
             linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = mean(vectors_out_df$vector_cor), colour = "black", 
             linetype = "solid", linewidth = 0.8) +
  scale_x_continuous(limits = c(-1,1)) +
  # lims(x=c(-1,1)) +
  theme_linedraw() +
  labs(title = "Density of sexual dimorphism vector correlations",
    x = 'Vector correlation (r)',
    y = 'Density of pairwise species comparisons') +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=16))
# ggsave('figures/direction/pheno_direction.png', p_vectors, width = 8, height = 6, dpi = 600, bg = "white")

vectors_out_df |> filter(vector_cor <= 0)
# 14 pairwise comparisons negative, 10 include E. arcaeformis
nrow(vectors_out_df)

###
# Combine Proc distances & Proc var into a single table
###
procrustes_table <- procrust_d |> mutate(species = species.a,
                     Distance = d,
                     Zscore = Z) |> 
  select(species, Distance, Zscore, pvalue) |> 
  left_join(
    disparity_t_wide |> select(species,disparity_F,lower_F,upper_F,
                               disparity_M,lower_M,upper_M,overlap),
    by='species') |> left_join(
      table1,by=c('species'='Species')
    ) |> mutate (Species = print_species,
                 Females = females,
                 Males = males) |> 
  select(Species,Subgenus,Females,disparity_F,lower_F,upper_F,
         Males,disparity_M,lower_M,upper_M,overlap,Distance,Zscore,pvalue) |> 
  arrange(Subgenus,Species)
# write.csv(procrustes_table,'figures/tables/procrustes_table.csv')


### PCA on mean shapes ###
# Perform PCA
pca <- prcomp(mean_shapes_mat, scale. = FALSE)

# Put PCA scores into a data frame
pc_scores <- data.frame(pca$x)

# Add the group names as a column
pc_scores$Group <- rownames(pc_scores)

# Extract species label from group names
# Assumes format: "SpeciesName_Sex", e.g., "Epio1.M" or "Epio2.F"
pc_scores$SpeciesLabel <- sapply(pc_scores$Group, function(x) {
  # Split by "_" if it exists
  parts <- strsplit(x, "_")[[1]]
  # Take second part if "_" exists, otherwise first part
  label <- if(length(parts) > 1) parts[2] else parts[1]
  # Remove the ".M" or ".F" at the end
  gsub("\\..$", "", label)
})

# Keep species for joining lines and sex
pc_scores$Species <- sapply(strsplit(as.character(pc_scores$Group), "\\."), `[`, 1)
pc_scores$Sex <- sapply(strsplit(as.character(pc_scores$Group), "\\."), `[`, 2)

# Separate male and female points
male_scores <- pc_scores %>% filter(Sex == "M")
female_scores <- pc_scores %>% filter(Sex == "F")

# Create dataframe for connecting lines
lines_df <- left_join(female_scores, male_scores, by = "Species", suffix = c(".F", ".M"))
lines_df <- lines_df %>% left_join(procrust_d, by = c("SpeciesLabel.M" = "species.a"))
# variance explained by pc1 and 2
var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
xy.limits <- range(pc_scores$PC1)
unique(lines_df$proc_d_sig)
unique(pc_scores$Sex)
pc_scores_all$Sex <- as.factor(pc_scores_all$Sex)
# Below - only needed if plotting density in background
# pc_scores_all$Sex =  with(pc_scores_all, factor(Sex, levels = rev(levels(Sex))))
# pc_scores$Sex <- as.factor(pc_scores$Sex)
# pc_scores$Sex =  with(pc_scores, factor(Sex, levels = rev(levels(Sex))))
# Plot with PC axis labels including variance
pca_FM <- ggplot() +
  # Density background of ALL PCA plot
  # stat_density_2d(data = pc_scores_all, aes(x = PC1, y = PC2, fill = Sex),
  #                 geom = "polygon", color="white", alpha=0.05,
  #                 linewidth=0.05) +
  # Lines connecting male-female pairs
  geom_segment(data = lines_df,
               aes(x = PC1.F, y = PC2.F, xend = PC1.M, yend = PC2.M,
                   linetype = proc_d_sig), 
               color = 'gray10', linewidth=0.15) +
  # scale_alpha_manual(values=c(0.6,0.6),guide = 'none') +
  # scale_linetype_manual(values=c('dashed','solid'), guide = 'none') +
  scale_linetype_manual(values=c('79','solid'), guide = 'none') +
  geom_point(data = pc_scores, aes(x = PC1, y = PC2,
                                   shape = Sex, fill = Sex), color='black', size = 4.7) +
  scale_shape_manual(labels = c('Female','Male'), values = c(22,21), guide = 'legend') +
  scale_fill_manual(labels = c('Female','Male'), values = c('tomato','steelblue')) +
  # # Female points
  # geom_point(data = female_scores, aes(x = PC1, y = PC2), shape = 21, fill = "tomato", color='black', size = 4.7) +
  # # Male points
  # geom_point(data = male_scores, aes(x = PC1, y = PC2), shape = 22, fill = "steelblue",color='black', size = 4.7) +
  # # Species names at female dot
  # # geom_text(data = female_scores,
  # #           aes(x = PC1, y = PC2, label = SpeciesLabel),
  # #           size = 3, vjust = -1, hjust = 0.5) +
  # scale_x_continuous(limits=xy.limits) + 
  # scale_y_continuous(limits=xy.limits) + 
  theme_linedraw() +
  coord_fixed(ratio=1) +
  labs(
    # title = "PCA of Mean Shapes with Sexual Dimorphism Vectors",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)")) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = 'left')

# ggsave('figures/sex_dimorph_PCA/pca_FM.png', pca_FM, width = 10, height = 7, dpi = 300, bg = "white")

#######
## Create custom TPS heatmaps of Jacobian determinants for extremes on PCA
#######


# Convert GDF to 2D to find species*sexes at min and max of PC1 & 2
epio.df.mean_2d <- data.frame(two.d.array(Epio.df.mean$coords, sep = ".")) %>% 
  tibble::rownames_to_column("specimen") 

# PC1 min = 'propinqua.M'
PC1_min <- rownames(subset(pc_scores, PC1 == min(PC1)))
target <- matrix(unlist(epio.df.mean_2d %>% 
  filter(specimen == PC1_min) %>% 
  dplyr::select(!specimen)), ncol = 2, byrow = TRUE)
# target <- Epio.df.mean$coords[,,14]

# PC1 max = 'ahlstedti.F
PC1_max <- rownames(subset(pc_scores, PC1 == max(PC1)))
target <- matrix(unlist(epio.df.mean_2d %>% 
                   filter(specimen == PC1_max) %>% 
                   dplyr::select(!specimen)), ncol = 2, byrow = TRUE)
# PC2 min = 'personata.F'
PC2_min <- rownames(subset(pc_scores, PC2 == min(PC2)))
target <- matrix(unlist(epio.df.mean_2d %>% 
                   filter(specimen == PC2_min) %>% 
                   dplyr::select(!specimen)), ncol = 2, byrow = TRUE)
# PC2 max = 'lenior.M'
PC2_max <- rownames(subset(pc_scores, PC2 == max(PC2)))
target <- matrix(unlist(epio.df.mean_2d %>% 
                   filter(specimen == PC2_max) %>% 
                   dplyr::select(!specimen)), ncol = 2, byrow = TRUE)

ref <- mshape(Epio.df.mean$coords)

tps_fit <- function(ref, target){
  n <- nrow(ref)
  ## radial basis matrix
  r <- as.matrix(dist(ref))
  K <- r^2
  K[K==0] <- 1
  K <- K*log(K)
  K[is.na(K)] <- 0
  ## affine part
  P <- cbind(1, ref)
  O <- matrix(0,3,3)
  L <- rbind(cbind(K,P),cbind(t(P),O))
  Y <- rbind(target,matrix(0,3,2))
  coef <- solve(L,Y)
  list(ref=ref,coef=coef)}

#step 2
tps_predict <- function(model,newpts){
  newpts <- as.matrix(newpts)
  ref <- as.matrix(model$ref)
  coef <- model$coef
  n <- nrow(ref)
  r <- as.matrix(dist(rbind(newpts,ref)))
  r <- r[1:nrow(newpts),
         (nrow(newpts)+1):(nrow(newpts)+n)]
  U <- r^2
  U[U==0] <- 1
  U <- U*log(U)
  U[is.na(U)] <- 0
  P <- cbind(1,newpts)
  cbind(U,P)%*%coef}

#step 3
pad <- 0.02
x <- seq(min(ref[,1])-pad,
         max(ref[,1])+pad,
         length=50) # this determines resolution of grid
y <- seq(min(ref[,2])-pad,
         max(ref[,2])+pad,
         length=50) # this determines resolution of grid
grid <- expand.grid(x=x,y=y)

#step 4
fit <- tps_fit(ref,target)
warp <- tps_predict(fit,grid)

# step 5 - magnitude
grid$dx <- warp[,1]-grid$x
grid$dy <- warp[,2]-grid$y
grid$disp <- sqrt(grid$dx^2+grid$dy^2)

#step 7 - warped TPS grid
grid$wx <- warp[,1]
grid$wy <- warp[,2]
gridlines <- grid |>
  mutate(ix=rep(1:50,each=50),
         iy=rep(1:50,50))

# target specimen landmarks
target.df <- data.frame(
  x = target[,1],
  y = target[,2],
  landmark = 1:nrow(target))

# links
links <- Epio.df.mean$links
link.df <- data.frame(
  x    = target[links[,1],1],
  y    = target[links[,1],2],
  xend = target[links[,2],1],
  yend = target[links[,2],2])


# loop to build warped polygons
poly.list <- list()
id <- 1

for(i in 1:(max(gridlines$ix)-1)){
  
  for(j in 1:(max(gridlines$iy)-1)){
    
    p1 <- filter(gridlines, ix==i,   iy==j)
    p2 <- filter(gridlines, ix==i+1, iy==j)
    p3 <- filter(gridlines, ix==i+1, iy==j+1)
    p4 <- filter(gridlines, ix==i,   iy==j+1)
    # needed for Jacobian Determinant
    v1 <- c(p2$x-p1$x,
            p2$y-p1$y)
    v2 <- c(p4$x-p1$x,
            p4$y-p1$y)
    V1 <- c(p2$wx-p1$wx,
            p2$wy-p1$wy)
    V2 <- c(p4$wx-p1$wx,
            p4$wy-p1$wy)
    A <- cbind(v1,v2)
    B <- cbind(V1,V2)
    C <- B %*% solve(A)
    jac <- det(C)
    
    poly.list[[id]] <-
      data.frame(
        x = c(p1$wx,p2$wx,p3$wx,p4$wx),
        y = c(p1$wy,p2$wy,p3$wy,p4$wy),
        group = id,
        disp = mean(c(p1$disp,p2$disp,p3$disp,p4$disp)),
        jac = log(jac))
    id <- id+1
  }
}
poly.df <- bind_rows(poly.list)


p_jac <- ggplot() +
  ## Heat map
  # geom_polygon(data=poly.df, aes(x,y,
  #                  group=group, fill=disp),
  #              colour='gray40', linewidth = 0.4, alpha = 0.5) +
  geom_polygon(data=poly.df, aes(x,y,
                                 group=group, fill=jac),
               colour='gray30', linewidth = 0.4, alpha = 0.95) +
  ## Links
  geom_segment(data = link.df,
               aes(x, y, xend = xend, yend = yend),
               linewidth = 2.3, colour = "black") +
  ## Landmarks
  geom_point(data = target.df, aes(x, y),
             shape = 21, fill = "white",
             colour = "black",stroke = 1.8, size = 7) +
  coord_equal() +
  # scale_fill_viridis_c(option="turbo") +
  # scale_fill_gradient2(low = "#053061",
  #         mid = "#F7F7F7", high = "#67001F", midpoint = 0) +
  scale_fill_gradient2(low = "dodgerblue4",
                       mid = "white", high = "firebrick3", midpoint = 0) +
  theme_void() +
  theme(legend.position = 'off')

# ggsave('figures/pca_FM_PC2_max.png', p_jac, width = 10, height = 7, dpi = 300, bg = "white")




################################################################################
################################################################################
## --- AIM 3: Patterns of female shape --- ##
################################################################################
################################################################################



## FEMALE ONLY
# Perform PCA
mean_shapes_2d <- data.frame(two.d.array(mean_shapes, sep = "."))
female <- mean_shapes_2d %>% filter(str_detect(row.names(mean_shapes_2d), ".F"))

pca_F <- prcomp(female, scale. = FALSE)

# Put PCA scores into a data frame
pc_scores_F <- data.frame(pca_F$x)

# Variation explained by PCs
var_explained_F <- round(100 * pca_F$sdev^2 / sum(pca_F$sdev^2), 1)
var_unrounded <- pca_F$sdev^2 / sum(pca_F$sdev^2)

sum(var_unrounded[1:6])
# 0.9604824
# Retain only high variation PCs 95%
pc_scores_F_95 <- pc_scores_F[,1:6]


#
### Hierarchical clustering ----
#

## Determine number of clusters
# Gap statistic
gap_stat <- clusGap(pc_scores_F_95, 
               FUN = function(x, k){hcut(x, k = k, hc_method = "ward.D2")},
  K.max = 9, B = 1000)
fviz_gap_stat(gap_stat)
maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
# 4 clusters

# Silhouette width
fviz_nbclust(pc_scores_F_95, FUNcluster = hcut, 
             method = "silhouette", 
             k.max=9)
# 4 clusters

## Assess stability of clusters
#clusterboot
cb_4 <- clusterboot(pc_scores_F_95, bootmethod = "boot", B = 1000,
  distances = FALSE, clustermethod = hclustCBI, method = "ward.D2",
  scaling = FALSE, k = 4)
cb_4
cb_4$result

# Dendrogram creation
# pvclust
female_dend <- pvclust(t(pc_scores_F_95), method.hclust = "ward.D2", 
                       method.dist = "euclidean", nboot = 1000)
plot(female_dend, main = "Female only", xlab = "", sub = "", print.num=FALSE)
# pvrect(female_dend)
# pvpick(female_dend)

# hclust
# Compute Dissimilarity Matrix
female_dist <- dist(pc_scores_F_95, method = "euclidean")
female_hclust <- hclust(female_dist, method = "ward.D2")
female_hclust$labels <- gsub("\\.F$", "", female_hclust$labels)
fviz_dend(female_hclust, k = 4, rect = TRUE, show_labels = TRUE, horiz = TRUE, type = "rectangle", lower_rect = -.1, main = NULL)
# pvclust and hclust produce same topology (euclidean distance with ward.D2 linkage)
plot(female_dend)
fviz_dend(female_hclust, k = 4, rect = TRUE, show_labels = TRUE, horiz = TRUE, type = "rectangle", lower_rect = -.1, main = NULL)

####
## Combine dendrogram with pvclust AU bootstraps for downstream ggtree
####

# create df with necessary metadata for plotting tree
species.subgenus <- unique(raw %>% mutate(species.subgenus = paste(raw$Species,raw$Subgenus)) %>% 
  dplyr::select(species.subgenus)) 

# Clean up data for plotting
species.subgenus <- separate(species.subgenus, species.subgenus, c('species', "subgenus"),sep = ' ')
species.subgenus <- species.subgenus %>% mutate(full.name = paste0('E. (',species.subgenus$subgenus,') ',species.subgenus$species))

# dendrogram object
hc <- female_dend$hclust
# convert dendrogram to phylo
phy <- as.phylo(hc)
# create function to get descendants / taxon sets
get_hclust_members <- function(hc) {
  n <- length(hc$order)
  merges <- hc$merge
  members <- vector("list", nrow(merges))
  for(i in seq_len(nrow(merges))) {
    left <- merges[i,1]
    right <- merges[i,2]
    left_members <- if(left < 0) {
      hc$labels[-left]
    } else {
      members[[left]]
    }
    right_members <- if(right < 0) {
      hc$labels[-right]
    } else {
      members[[right]]
    }
    members[[i]] <- sort(c(left_members, right_members))
  }
  members
}
# taxa sets of pvclust dendrogram
hc_members <- get_hclust_members(hc)
# get nodes of converted dendrogram
phy_nodes <- (Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
# taxa sets of converted phylo object
phy_members <- lapply(phy_nodes, function(node) {
  sort(extract.clade(phy, node)$tip.label)
})
# node matches
matches <- sapply(phy_members, function(x) {
  which(sapply(hc_members, identical, x))
})
# AU bootstraps
au <- female_dend$edges$au
# df to add to ggtree
node_df <- data.frame(
  node = phy_nodes,
  AU = au[matches])

# Assemble additional metadata for PCA & ggtree plotting
# Final k to be used for plotting
k <- 4
# add clusters
female_clusters <- cutree(female_hclust, k = k)  # tip name -> cluster number
species.subgenus$cluster <- factor(female_clusters[species.subgenus$species]) 
# Add additional F traits (if needed)
# Ftraits <- Ftraits %>% dplyr::select(taxa,denticles,recurve,print_names)
# Add tip.labels to df
species.subgenus$label <- paste0(species.subgenus$species,'.F')
# label column needs to be first in df
species.subgenus <- species.subgenus %>% dplyr::select(label,full.name, cluster)

###
## Configure ggtree
###
base_ggtree <- ggtree(phy, ladderize = TRUE) %<+% node_df %<+% species.subgenus
base_ggtree@data
# Check nodes for flipping rotation
base_ggtree + geom_nodelab(aes(label=node), hjust=-.3) +
  geom_tiplab(aes(label=full.name,color=cluster),
              fontface=4, size = 3,offset = 0.005)
# Swivel nodes for better labeling
#tree_plot <- flip(tree_plot, 30, 31)
# tree_plot <- flip(base_ggtree, 46, 45)
# tree_plot <- flip(tree_plot, 30, 31)
# tree_plot <- flip(tree_plot, 34, 35)
tree_plot <- rotate(base_ggtree, 44)

# clade colors
cluster_colors <- c(
  "1" = "#0274BD", # 33 = cluster 1
  "2" = "#25A18E", # 32 = cluster 2
  "3" = "#AC855E", # 35 = cluster 3
  "4" = "#E37C78") # 34 = cluster 4

# saveRDS(tree_plot, file = 'data/raw/clust_dendro_ggtree.Rdata')

# Load mean Epio geomorph dataframe
tree_plot <- readRDS(file = 'data/raw/clust_dendro_ggtree.Rdata')


tree_plot + geom_nodelab(aes(label=node), hjust=-.3)
# need to change these groups below manually to match clusters
tree_fig <- groupClade(tree_plot, c(32,33,35,34), group_name = "group") + 
  geom_tree(aes(color = group), linewidth = 4) +
  scale_color_manual(values = c("black","#25A18E", "#0274BD","#E37C78","#AC855E")) +
  new_scale_color() +
  geom_tiplab(aes(label=full.name,
                  color=cluster),
              fontface=4, size = 8,
              offset = 0.0025,
              align= FALSE) +
  geom_text2(aes(subset = !isTip, label = round(AU,2)),
             hjust = 1.2,
             vjust = -0.7,
             size = 5.5) +
  scale_color_manual(values = cluster_colors, na.value = "black") +
  theme_tree2() +
  theme(legend.position = 'off',
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  coord_cartesian(clip = "off") +
  xlim_tree(0.27)

ggsave(filename = "figures/tree_fig_test.png", plot = tree_fig, 
       width = 10, height = 20, dpi = 300)


#####
## MEAN SHAPES FOR TPS GRIDS OF FIGURE 3
#####

mean_shapes_2d <- data.frame(two.d.array(mean_shapes, sep = "."))
#female mean shapes
female <- mean_shapes_2d %>% filter(str_detect(row.names(mean_shapes_2d), ".F"))
female$label <- rownames(female)
# not sure if this line works
female <- female |> left_join(tree_plot@data |> dplyr::select(label, cluster), by = 'label')
# male mean shapes
male <- mean_shapes_2d %>% filter(str_detect(row.names(mean_shapes_2d), ".M"))
male$label <- paste0(sapply(rownames(male), gsub, pattern='.M', replacement='',USE.NAMES = F),'.F')
male <- male |> left_join(female |> dplyr::select(label, cluster), by='label')


### cluster 1
clust1 <- as.matrix(female |> filter(cluster == 1) |> dplyr::select(X1:X66))
clust1_n <- nrow(clust1)
clust1_array <- array(NA, dim = c(33, 2, clust1_n))
for(i in 1:clust1_n){
  clust1_array[,,i] <- matrix(clust1[i,], ncol = 2, byrow = TRUE)}
# female mean shape
mshape_1_tor_f <- mshape(clust1_array)
# male mean shape
clust1_m <- as.matrix(male |> filter(cluster == 1) |> dplyr::select(X1:X66))
clust1_m_array <- array(NA, dim = c(33, 2, clust1_n))
for(i in 1:clust1_n){
  clust1_m_array[,,i] <- matrix(clust1_m[i,], ncol = 2, byrow = TRUE)}
mshape_1_tor_m <- mshape(clust1_m_array)


### cluster 2
clust2 <- as.matrix(female |> filter(cluster == 2) |> dplyr::select(X1:X66))
clust2_n <- nrow(clust2)
clust2_array <- array(NA, dim = c(33, 2, clust2_n))
for(i in 1:clust2_n){
  clust2_array[,,i] <- matrix(clust2[i,], ncol = 2, byrow = TRUE)}
# female mean shape
mshape_2_plag_f <- mshape(clust2_array)
# male mean shape
clust2_m <- as.matrix(male |> filter(cluster == 2) |> dplyr::select(X1:X66))
clust2_m_array <- array(NA, dim = c(33, 2, clust2_n))
for(i in 1:clust2_n){
  clust2_m_array[,,i] <- matrix(clust2_m[i,], ncol = 2, byrow = TRUE)}
mshape_2_plag_m <- mshape(clust2_m_array)


### cluster 3
clust3 <- as.matrix(female |> filter(cluster == 3) |> dplyr::select(X1:X66))
clust3_n <- nrow(clust3)
clust3_array <- array(NA, dim = c(33, 2, clust3_n))
for(i in 1:clust3_n){
  clust3_array[,,i] <- matrix(clust3[i,], ncol = 2, byrow = TRUE)}
# female mean shape
mshape_3_epi_f <- mshape(clust3_array)
# male mean shape
clust3_m <- as.matrix(male |> filter(cluster == 3) |> dplyr::select(X1:X66))
clust3_m_array <- array(NA, dim = c(33, 2, clust3_n))
for(i in 1:clust3_n){
  clust3_m_array[,,i] <- matrix(clust3_m[i,], ncol = 2, byrow = TRUE)}
mshape_3_epi_m <- mshape(clust3_m_array)


### cluster 4
clust4 <- as.matrix(female |> filter(cluster == 4) |> dplyr::select(X1:X66))
clust4_n <- nrow(clust4)
clust4_array <- array(NA, dim = c(33, 2, clust4_n))
for(i in 1:clust4_n){
  clust4_array[,,i] <- matrix(clust4[i,], ncol = 2, byrow = TRUE)}
# female mean shape
mshape_4_pil_f <- mshape(clust4_array)
# male mean shape
clust4_m <- as.matrix(male |> filter(cluster == 4) |> dplyr::select(X1:X66))
clust4_m_array <- array(NA, dim = c(33, 2, clust4_n))
for(i in 1:clust4_n){
  clust4_m_array[,,i] <- matrix(clust4_m[i,], ncol = 2, byrow = TRUE)}
mshape_4_pil_m <- mshape(clust4_m_array)


#### Mean shape of grouped clusters

## clusters 2,3,4
clust234 <- as.matrix(female |> filter(cluster %in% c(2,3,4)) |> dplyr::select(X1:X66))
clust234_n <- nrow(clust234)
clust234_array <- array(NA, dim = c(33, 2, clust234_n))
for(i in 1:clust234_n){
  clust234_array[,,i] <- matrix(clust234[i,], ncol = 2, byrow = TRUE)}
mshape_clust234 <- mshape(clust234_array)

## clusters 1,3,4
clust134 <- as.matrix(female |> filter(cluster %in% c(1,3,4)) |> dplyr::select(X1:X66))
clust134_n <- nrow(clust134)
clust134_array <- array(NA, dim = c(33, 2, clust134_n))
for(i in 1:clust134_n){
  clust134_array[,,i] <- matrix(clust134[i,], ncol = 2, byrow = TRUE)}
mshape_clust134 <- mshape(clust134_array)

## clusters 1,2,4
clust124 <- as.matrix(female |> filter(cluster %in% c(1,2,4)) |> dplyr::select(X1:X66))
clust124_n <- nrow(clust124)
clust124_array <- array(NA, dim = c(33, 2, clust124_n))
for(i in 1:clust124_n){
  clust124_array[,,i] <- matrix(clust124[i,], ncol = 2, byrow = TRUE)}
mshape_clust124 <- mshape(clust124_array)

## clusters 1,2,3
clust123 <- as.matrix(female |> filter(cluster %in% c(1,2,3)) |> dplyr::select(X1:X66))
clust123_n <- nrow(clust123)
clust123_array <- array(NA, dim = c(33, 2, clust123_n))
for(i in 1:clust123_n){
  clust123_array[,,i] <- matrix(clust123[i,], ncol = 2, byrow = TRUE)}
mshape_clust123 <- mshape(clust123_array)


# Create TPS grids for -----
# mshape_1_tor_f / mshape_1_tor_m
# mshape_1_tor_f / clust234_f
# mshape_2_plag_f / mshape_2_plag_m
# mshape_2_plag_f / clust134_f
# mshape_3_epi_f / mshape_3_epi_m
# mshape_3_epi_f / clust124_f
# mshape_4_pil_f / mshape_4_pil_m
# mshape_4_pil_f / clust123_f


target.list <- array(NA, dim = c(33, 2, 8))

target.list[,,1] <- mshape_1_tor_f
target.list[,,2] <- mshape_1_tor_f
target.list[,,3] <- mshape_2_plag_f
target.list[,,4] <- mshape_2_plag_f
target.list[,,5] <- mshape_3_epi_f
target.list[,,6] <- mshape_3_epi_f
target.list[,,7] <- mshape_4_pil_f
target.list[,,8] <- mshape_4_pil_f

ref.list <- array(NA, dim = c(33, 2, 8))

ref.list[,,1] <- mshape_1_tor_m
ref.list[,,2] <- mshape_clust234
ref.list[,,3] <- mshape_2_plag_m
ref.list[,,4] <- mshape_clust134
ref.list[,,5] <- mshape_3_epi_m
ref.list[,,6] <- mshape_clust124
ref.list[,,7] <- mshape_4_pil_m
ref.list[,,8] <- mshape_clust123

for(z in 1:8){
  
ref <- ref.list[,,z]
target <- target.list[,,z]

tps_fit <- function(ref, target){
  n <- nrow(ref)
  ## radial basis matrix
  r <- as.matrix(dist(ref))
  K <- r^2
  K[K==0] <- 1
  K <- K*log(K)
  K[is.na(K)] <- 0
  ## affine part
  P <- cbind(1, ref)
  O <- matrix(0,3,3)
  L <- rbind(cbind(K,P),cbind(t(P),O))
  Y <- rbind(target,matrix(0,3,2))
  coef <- solve(L,Y)
  list(ref=ref,coef=coef)}

#step 2
tps_predict <- function(model,newpts){
  newpts <- as.matrix(newpts)
  ref <- as.matrix(model$ref)
  coef <- model$coef
  n <- nrow(ref)
  r <- as.matrix(dist(rbind(newpts,ref)))
  r <- r[1:nrow(newpts),
         (nrow(newpts)+1):(nrow(newpts)+n)]
  U <- r^2
  U[U==0] <- 1
  U <- U*log(U)
  U[is.na(U)] <- 0
  P <- cbind(1,newpts)
  cbind(U,P)%*%coef}

#step 3
pad <- 0.02
x <- seq(min(ref[,1])-pad,
         max(ref[,1])+pad,
         length=50) # this determines resolution of grid
y <- seq(min(ref[,2])-pad,
         max(ref[,2])+pad,
         length=50) # this determines resolution of grid
grid <- expand.grid(x=x,y=y)

#step 4
fit <- tps_fit(ref,target)
warp <- tps_predict(fit,grid)

# step 5 - magnitude
grid$dx <- warp[,1]-grid$x
grid$dy <- warp[,2]-grid$y
grid$disp <- sqrt(grid$dx^2+grid$dy^2)

#step 7 - warped TPS grid
grid$wx <- warp[,1]
grid$wy <- warp[,2]
gridlines <- grid |>
  mutate(ix=rep(1:50,each=50),
         iy=rep(1:50,50))

# target specimen landmarks
target.df <- data.frame(
  x = target[,1],
  y = target[,2],
  landmark = 1:nrow(target))

# links
links <- Epio.df.mean$links
link.df <- data.frame(
  x    = target[links[,1],1],
  y    = target[links[,1],2],
  xend = target[links[,2],1],
  yend = target[links[,2],2])


# loop to build warped polygons
poly.list <- list()
id <- 1

for(i in 1:(max(gridlines$ix)-1)){
  
  for(j in 1:(max(gridlines$iy)-1)){
    
    p1 <- filter(gridlines, ix==i,   iy==j)
    p2 <- filter(gridlines, ix==i+1, iy==j)
    p3 <- filter(gridlines, ix==i+1, iy==j+1)
    p4 <- filter(gridlines, ix==i,   iy==j+1)
    # needed for Jacobian Determinant
    v1 <- c(p2$x-p1$x,
            p2$y-p1$y)
    v2 <- c(p4$x-p1$x,
            p4$y-p1$y)
    V1 <- c(p2$wx-p1$wx,
            p2$wy-p1$wy)
    V2 <- c(p4$wx-p1$wx,
            p4$wy-p1$wy)
    A <- cbind(v1,v2)
    B <- cbind(V1,V2)
    C <- B %*% solve(A)
    jac <- det(C)
    
    poly.list[[id]] <-
      data.frame(
        x = c(p1$wx,p2$wx,p3$wx,p4$wx),
        y = c(p1$wy,p2$wy,p3$wy,p4$wy),
        group = id,
        disp = mean(c(p1$disp,p2$disp,p3$disp,p4$disp)),
        jac = log(jac))
    id <- id+1
  }
}
poly.df <- bind_rows(poly.list)


p_jac <- ggplot() +
  ## Heat map
  # geom_polygon(data=poly.df, aes(x,y,
  #                  group=group, fill=disp),
  #              colour='gray40', linewidth = 0.4, alpha = 0.5) +
  geom_polygon(data=poly.df, aes(x,y,
                                 group=group, fill=jac),
               colour='gray30', linewidth = 0.4, alpha = 0.95) +
  ## Links
  geom_segment(data = link.df,
               aes(x, y, xend = xend, yend = yend),
               linewidth = 2.3, colour = "black") +
  ## Landmarks
  geom_point(data = target.df, aes(x, y),
             shape = 21, fill = "white",
             colour = "black",stroke = 1.8, size = 7) +
  coord_equal() +
  # scale_fill_viridis_c(option="turbo") +
  # scale_fill_gradient2(low = "#053061",
  #         mid = "#F7F7F7", high = "#67001F", midpoint = 0) +
  scale_fill_gradient2(low = "dodgerblue4",
                       mid = "white", high = "firebrick3", midpoint = 0) +
  theme_void() +
  theme(legend.position = 'off')

# ggsave(paste0('figures/TPS/',z,'.png'), p_jac, width = 10, height = 7, dpi = 300,
#        bg = "white")
}


####
# FEMALE PCA - Plots
####
# Variation explained by PCs
var_unrounded <- pca_F$sdev^2 / sum(pca_F$sdev^2)
var_unrounded <- data.frame(var_unrounded = var_unrounded[1:20])
var_unrounded$cum_pro <- c(sum(var_unrounded$var_unrounded[1]),
                           sum(var_unrounded$var_unrounded[1:2]),
                           sum(var_unrounded$var_unrounded[1:3]),
                           sum(var_unrounded$var_unrounded[1:4]),
                           sum(var_unrounded$var_unrounded[1:5]),
                           sum(var_unrounded$var_unrounded[1:6]),
                           sum(var_unrounded$var_unrounded[1:7]),
                           sum(var_unrounded$var_unrounded[1:8]),
                           sum(var_unrounded$var_unrounded[1:9]),
                           sum(var_unrounded$var_unrounded[1:10]),
                           sum(var_unrounded$var_unrounded[1:11]),
                           sum(var_unrounded$var_unrounded[1:12]),
                           sum(var_unrounded$var_unrounded[1:13]),
                           sum(var_unrounded$var_unrounded[1:14]),
                           sum(var_unrounded$var_unrounded[1:15]),
                           sum(var_unrounded$var_unrounded[1:16]),
                           sum(var_unrounded$var_unrounded[1:17]),
                           sum(var_unrounded$var_unrounded[1:18]),
                           sum(var_unrounded$var_unrounded[1:19]),
                           sum(var_unrounded$var_unrounded[1:20]))
var_unrounded$PC <- c(1:20)
# Scree plot to determine how many PCs
scree_p <- ggplot(data=var_unrounded) +
  geom_line(aes(x=PC, y=var_unrounded, group = 1)) +
  geom_point(aes(x=PC, y=var_unrounded),fill='black',  color='black',shape=21,
             size=3) +
  geom_line(aes(x=PC, y=cum_pro, group = 1),linetype = 'dashed') +
  geom_point(aes(x=PC, y=cum_pro),fill='white',  color='black',shape=21,
             size=3) +
  geom_vline(xintercept = 6,color='red') +
  labs(y='Variance explained',x='Principal component', 
       title = 'Variance (proportion & cumulative) per PC') +
  theme_linedraw() +
  theme(legend.position = 'off',
        plot.title = element_text(hjust = 0.5,size=16),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=16))

# Create convex hulls
pc_scores_F_95$label <- rownames(pc_scores_F_95)
pc_scores_F_95 <- pc_scores_F_95 |> left_join(tree_plot@data, by = 'label')

hull_data <- pc_scores_F_95 %>% group_by(cluster) %>% 
  slice(chull(PC1, PC2))
hull_data_2 <- pc_scores_F_95 %>% group_by(cluster) %>% 
  slice(chull(PC1, PC3))
hull_data_3 <- pc_scores_F_95 %>% group_by(cluster) %>% 
  slice(chull(PC3, PC2))

# PCA 1 & 2
PCA_12_F <- ggplot(data = pc_scores_F_95, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill=cluster), 
             size = 4.5, shape=21,
             color='black', alpha = 0.9) +
  scale_color_manual(values = cluster_colors, na.value = "black") +
  scale_fill_manual(values = cluster_colors, na.value = "black") +
  geom_polygon(data = hull_data, aes(fill = cluster, color=cluster), alpha = 0.2) +
  labs(y=paste0("PC2 (", var_explained_F[2], "%)"),
       x=paste0("PC1 (", var_explained_F[1], "%)")) +
  theme_linedraw()+
  ylim(min(pc_scores_F_95$PC2,pc_scores_F_95$PC3),max(pc_scores_F_95$PC2,pc_scores_F_95$PC3)) +
  xlim(min(pc_scores_F_95$PC1,pc_scores_F_95$PC3),max(pc_scores_F_95$PC1,pc_scores_F_95$PC3)) +
  coord_fixed(ratio=1) +
  theme(legend.position = 'off',
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=16))

# PCA 1 & 3 (for supplementary)
PCA_13_F <- ggplot(data = pc_scores_F_95, aes(x = PC1, y = PC3)) +
  geom_point(aes(fill=cluster), 
             size = 4.5, shape=21,
             color='black', alpha = 0.9) +
  scale_color_manual(values = cluster_colors, na.value = "black") +
  scale_fill_manual(values = cluster_colors, na.value = "black") +
  geom_polygon(data = hull_data_2, aes(fill = cluster, color=cluster), alpha = 0.2) +
  labs(y=paste0("PC3 (", var_explained_F[3], "%)"),
       x=paste0("PC1 (", var_explained_F[1], "%)")) +
  theme_linedraw()+
  ylim(min(pc_scores_F_95$PC2,pc_scores_F_95$PC3),max(pc_scores_F_95$PC2,pc_scores_F_95$PC3)) +
  xlim(min(pc_scores_F_95$PC1,pc_scores_F_95$PC3),max(pc_scores_F_95$PC1,pc_scores_F_95$PC3)) +
  coord_fixed(ratio=1) +
  theme(legend.position = 'off',
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=16))

# PCA 2 & 3 (for supplementary)
PCA_23_F <- ggplot(data = pc_scores_F_95, aes(x = PC3, y = PC2)) +
  geom_point(aes(fill=cluster), 
             size = 4.5, shape=21,
             color='black', alpha = 0.9) +
  scale_color_manual(values = cluster_colors, na.value = "black") +
  scale_fill_manual(values = cluster_colors, na.value = "black") +
  geom_polygon(data = hull_data_3, aes(fill = cluster, color=cluster), alpha = 0.2) +
  labs(y=paste0("PC2 (", var_explained_F[2], "%)"),
       x=paste0("PC3 (", var_explained_F[3], "%)")) +
  theme_linedraw()+
  ylim(min(pc_scores_F_95$PC2,pc_scores_F_95$PC3),max(pc_scores_F_95$PC2,pc_scores_F_95$PC3)) +
  xlim(min(pc_scores_F_95$PC1,pc_scores_F_95$PC3),max(pc_scores_F_95$PC1,pc_scores_F_95$PC3)) +
  coord_fixed(ratio=1) +
  theme(legend.position = 'off',
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=16))
unique(Epio.df_all$species)

female_PCA_patch <- PCA_12_F + PCA_23_F + PCA_13_F + plot_spacer() 
female_PCA_patch <- female_PCA_patch + 
  inset_element(scree_p, left = 0.0, bottom = 0,right = 1,top=1) +
  plot_annotation(tag_levels = 'A')
# ggsave('figures/female_PCA_patch.png', female_PCA_patch, width = 12, height = 10, dpi = 300, bg = "white")