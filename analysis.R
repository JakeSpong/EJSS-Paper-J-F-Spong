usethis::use_git()

####load libraries----

#for generating maps of sample locations
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggmap)
#for data manipulation
library(dplyr)
library(tidyverse) 
#for statistics
library(stats)
library(car)
#for fitting decay curves
library(broom)
#for diversity analysis
library(vegan)
#for Moran's I to test pseudoreplication
library(spdep)
#to decompose drivers of betadiversity
library(betapart)
#for graphing
library(ggplot2)
library(ggpubr)
#for significance letters
library(multcompView)
#for file path referencing
library(here)
#to display numbers as desired when graphing
library(scales)
#to save figures as .svg
library(svglite)

#### generate maps of sample locations (Figure 1) ----
#to get google satellite layer, Go to the Google Cloud Console.  Create a new project (or use an existing one). Enable the Maps Static API and Geocoding API.  Get your API key.

d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 
# Register your Google API key
register_google(key = "AIzaSyA42ZmX_RGv9dTA7PFy4UR0YIjxW3x0rUE")

# Load the UK map from the rnaturalearth package
uk_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "United Kingdom")

#extract he sample coordinates from the dataframe
sample_coordinates <- data.frame(
  Longitude = d$LongitudeE,
  Latitude = d$LatitudeN,
  Habitat = d$Habitat,
  Vegetation = d$Vegetation
)
# Calculate the bounding box surrounding the 30 sample coordinates
longitude_range_sample <- range(sample_coordinates$Longitude)
latitude_range_sample <- range(sample_coordinates$Latitude)

# Plot the UK map with only the bounding box
UK_Hawes <- ggplot(data = uk_map) +
  geom_sf(fill = "lightblue") +  # UK map with light blue fill
  geom_rect(
    aes(xmin = longitude_range_sample[1], xmax = longitude_range_sample[2], ymin = latitude_range_sample[1], ymax = latitude_range_sample[2]), 
    color = "red", fill = NA, size = 2) +  # Add bounding box around the points
  theme_minimal()
#show the map
#show(UK_Hawes)
#save the UK map with Haweswater highlighted
#ggsave("Figures/UK_Haweswater_boundingbox.svg", width = 8, height = 6)


# Get the Google Satellite map for the zoomed-in area based on the 30 sample coordinates
satellite_map <- get_map(location = c(lon = mean(sample_coordinates$Longitude), 
                                      lat = mean(sample_coordinates$Latitude)),
                         zoom = 15,  # Adjust zoom level for the desired zoom
                         maptype = "satellite",  # Use the "satellite" map type
                         source = "google")
#plot the map
sample_coords_map <- ggmap(satellite_map) +
  geom_point(data = sample_coordinates, aes(x = Longitude, y = Latitude, 
                                            color = Habitat, shape = Vegetation), size = 2) +
  scale_color_manual(values = c("Grassland" = "green", "Heathland" = "purple", "Woodland" = "brown")) +  # Color by habitat
  scale_shape_manual(values = c("Bracken Present" = 16, "Bracken Absent" = 17)) +  # Shape by vegetation
  theme_minimal() + 
  theme(legend.position = "right") 
#show the map
#show(sample_coords_map)
# Save the final map to a file
#ggsave("Figures/sample_locations.svg", plot = last_plot(), width = 7.5, height = 6, dpi = 300)


#combine the two maps into a single figure
maps_figure <- ggarrange(UK_Hawes, sample_coords_map, 
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1,
                         #the width of each panel of the multifigure plot
                         widths = c(2,5))
#show the plot in the Plots window
#show(maps_figure)
#save the figure
ggsave("Figures/Maps_panel_figure.svg", plot = last_plot(), width = 7.5, height = 6, dpi = 300)


#### vegetation alpha diversity (Figure 2, Supp. Figure 1) ----
d <- readr::read_csv(here::here("Data", "Vegetation-Survey-Data.csv")) 

#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)

#Simpson, Shannon and richness of plant communities

#rename Nonbracken so that it is plotted first, and bracken second
d$Vegetation[d$Vegetation == 'Non Bracken'] <- 'Bracken Absent'
d$Vegetation[d$Vegetation == 'Bracken'] <- 'Bracken Present'
#create a subset containing only species abundance values
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0

#look to see which species are growing next to bracken

veg_d <- d[c(1,2,3,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)]
veg_d <- veg_d[veg_d$Vegetation != "Bracken Absent", ]
#keep only columns (species) with at least one data value
veg_d <- veg_d[, colSums(veg_d != 0, na.rm = TRUE) > 0]
#reorder species from most to least abundant
veg_d <- cbind(
  veg_d[, 1:3],
  veg_d[, 4:16][, order(colSums(veg_d[, 4:16], na.rm = TRUE), decreasing = TRUE)],
  veg_d[, -(1:16), drop = FALSE]
)
print(colSums(veg_d[, 4:16], na.rm = TRUE))



#just the species counts
spe <- d[,-(1:9)]
spe <- spe[,-(23:24)]

#species richness
d$richness <- apply(spe[,]>0,1,sum)
#calculate diversity
d$shannon <- diversity(spe[,], "shannon")
d$simpson <- diversity(spe[,], "simpson")
d$evenness <- d$shannon / log(d$richness)
#replace NaN with 0 for evenness, where we have only 1 species
d$evenness[is.nan(d$evenness)] <- 0

#boxplot of grassland vs heathland, species richness
richness_bxp <- ggboxplot(d, x = "Habitat", y = "richness", color = "Vegetation", ylab = "Species Richness", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 


#graphing boxplots with bracken split from nonbracken
shannon_bxp <-ggboxplot(d, x = "Habitat", y = "shannon", color = "Vegetation", ylab = "Shannon Diversity", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#graphing boxplots with bracken split from nonbracken
simpson_bxp <- ggboxplot(d, x = "Habitat", y = "simpson", color = "Vegetation", ylab = "Simpson Diversity", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#boxplot of grassland vs heathland, species richness
evenness_bxp <- ggboxplot(d, x = "Habitat", y = "evenness", color = "Vegetation", ylab = "Pielou's Evenness", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#create panel figures
all_bxp1 <- ggarrange(richness_bxp, shannon_bxp, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1,
                      common.legend = TRUE)

all_bxp2 <- ggarrange(evenness_bxp, simpson_bxp, 
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1,
                      common.legend = TRUE)
show(all_bxp1)

ggsave(path = "Figures", paste0(Sys.Date(), "_Figure-2.svg"), all_bxp1)
ggsave(path = "Figures", paste0(Sys.Date(), "_Supp-1-Figure.svg"), all_bxp2)


#### vegetation beta diversity (Figure 3) ----
d <- readr::read_csv(here::here("Data", "Vegetation-Survey-Data.csv")) 

#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#replace row index with sample names
rownames(d) <- c("GB1", "GB2", "GB3", "GB4", "GB5", "GN10", "GN3", "GN4", "GN5","GN9", "UB1","UB10","UB4","UB6", "UB9", "UN1", "UN10", "UN3", "UN7", "UN8", "WB3","WB4","WB6","WB7","WB8","WN2","WN4","WN5","WN6","WN8")
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0
#make sure our variables are coded as factors
d$Habitat <- factor(d$Habitat, levels = c("Grassland", "Heathland", "Woodland"), labels = c("Grassland", "Heathland", "Woodland"))

#just the species counts
spe <- d[,-(1:9)]
spe <- spe[,-(23:24)]

#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#add extra space to the right of the plot so the legend fits
#par(mar=c(5, 4, 4, 15), xpd=TRUE)
#plot the NMDS
plot(example_NMDS, col = "white")


#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#point shapes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5)) 
text(-1.8,1.3, paste("Stress = ", round(example_NMDS$stress, 3)))

for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }
#specify legend manually for text sample IDs
#legend(1.7,0.9, legend = c("Grassland Bracken Present", "Grassland Bracken Asbent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"))

#save the file using Export -> Save As Image -> Width = 655, Height = 500 
#save the legend as a separate file and combine in Inkscape
plot(example_NMDS, col = "white")
#legend for point codes
legend(-1.7,0.9, legend=c("Grassland Bracken Present", "Grassland Bracken Absent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), pch = c(15, 0,16,1,17,2))


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- d[,(2:3)]
#run the permanova
veg_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
veg_permanova
#veg_permanova indicates that Habitat and Vegetation have significant effects, with habitat explainng 26.7% of the variation and Vegetation explaining 36.0 %

#run an anosim - when grouping by habitat
ano = anosim(as.matrix(spe), grouping = d$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(spe), grouping = d$Vegetation, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)



#### soil pH, moisture, total and water-extracable organic carbon, WEN and C:N (Figure 4)----
#order samples by ID alphabetically
d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)

#plot pH
pH_bxp <- ggboxplot(d, x = "Habitat", y = 'pH', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = "Soil pH") + theme(
    #remove x axis label
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#plot soil moisture
gsmc_bxp <- ggboxplot(d, x = "Habitat", y = "Water content as percentage wet soil mass", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = "Water Content (%)") + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    #remove legend
    legend.position = "none",
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#plot carbon
tsc_bxp <- ggboxplot(d, x = "Habitat", y = 'TotalCarbon', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = expression("Total Soil Carbon (g kg"^-1*")")) + theme(
    
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 
#plot water-extractable organic carbon
npoc_bxp <- ggboxplot(d, x = "Habitat", y = "`DOC Concentration (mgCperg)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = expression("WEOC Concentration (mg C g"^-1*")")) + theme(
    #remove x axis label
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 
#graphing boxplots. Units are ppm aka mg per kg
tnb_bxp <- ggboxplot(d, x = "Habitat", y = "`TNb Concentration (mgNperg)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = expression("TNb Concentration (mg N g"^-1*")")) + theme(
    #remove x axis label
    axis.title.x=element_blank(),
    # Remove panel border
    
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.75),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 
#C:N ratio
cnr_bxp <- ggboxplot(d, x = "Habitat", y = 'CNratio', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(x = "Habitat x Vegetation",
       y = "C:N ratio") + theme(
         #remove x axis label, tickes, labels
         axis.title.x=element_blank(),
         # Remove panel border
         panel.border = element_blank(),  
         # Remove panel grid lines
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         # Remove panel background
         panel.background = element_blank(),
         # Add axis line
         axis.line = element_line(colour = "black", linewidth = 0.5),
         #change colour and thickness of axis ticks
         axis.ticks = element_line(colour = "black", linewidth = 0.5),
         #change axis labels colour
         axis.title.y = element_text(colour = "black"),
         #change tick labels colour
         axis.text.y = element_text(colour = "black"),
       ) 

#save 770 wide, 960 high
panel_bxp <- ggarrange(pH_bxp, gsmc_bxp, tsc_bxp, npoc_bxp, tnb_bxp, cnr_bxp, 
                       labels = c("A", "B", "C", "D", "E", "F"),
                       ncol = 2, nrow = 3,
                       common.legend = TRUE, legend="top")

show(panel_bxp)

#### SUVA, alpha, total nitrogen(Supp. Figure 2) ----
#quality needs to be analysed; easiest thing to do is: check data for nice curve, fit 2 component exponential decay curve through data.  Get slope from this (spectral slope), and the beta paramter from this slope
#read in the processed, STANDARDISED absorbance data
d <- readr::read_csv(
  here::here("Data", "Standardised-Processed-WEOC-Absorbance-Data.csv")) 
#eliminate all absorbance values above 600 nm
d <- d %>% filter(`Wavelength (nm)` <= 600)

# Fit a two-component exponential decay curve through the data for all our curves
fitted <- d %>%
  nest(-`Sample ID`) %>%
  mutate(
    fit = map(data, ~nls(`Standardized Absorbance (L mg-1 cm-1)` ~ SSasymp(`Wavelength (nm)`, yf, y0, log_alpha), data = .)),
    tidied = map(fit, tidy),
    augmented = map(fit, augment),
  )
# Produce a table of fit parameters: y0, yf, alpha
table <- fitted %>% 
  unnest(tidied) %>% 
  select(`Sample ID`, term, estimate) %>% 
  spread(term, estimate) %>% 
  mutate(alpha = exp(log_alpha))
#display table of fit parameters
table

#plot each absorbance curve along with the line of best fit
augmented <- fitted %>% 
  unnest(augmented)
qplot(`Wavelength (nm)`, `Standardized Absorbance (L mg-1 cm-1)`, data = augmented, geom = 'point', colour = `Sample ID`) +
  geom_line(aes(y=.fitted))

#now we have extracted the paramters of our lines of best fit, what do we do now?!?!?!
table$Habitat <- c(rep("Grassland",10), rep("Heathland",10),rep("Woodland",10))
table$Vegetation <- c(rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5))

#boxplot the alpha, which describes the curve ie how quicky we go from low wavelength (high mass C compounds) to high wavelength (low mass C compounds)
#plot carbon:nitrogen ratio
alpha_bxp <- ggboxplot(table, x = "Habitat", y = 'alpha', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = "alpha parameter") + theme(
    #remove x axis label
    axis.title.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#SUVA data
#list of wavelengths of interest
wavelength_of_interest <- list(254)#list(250, 254, 260, 265, 272, 280, 285, 300, 340, 350, 365, 400, 436, 465)
#filter data to extract absorbance at wavelength of interest
abs <- d %>%
  filter(`Wavelength (nm)` == wavelength_of_interest) %>% #filter for specific wavelength
  select(`Sample ID`, `Standardized Absorbance (L mg-1 cm-1)`) #select the relevant columns

#add in habitat and vegetation factors
abs$Habitat <- c(rep("Grassland",10), rep("Heathland",10),rep("Woodland",10))
abs$Vegetation <- c(rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5),rep("Bracken Present",5), rep("Bracken Absent", 5))

abs <- as.data.frame(abs)

#rename Nonbracken so that it is plotted first, and bracken second
abs$Vegetation[abs$Vegetation == 'Non-Bracken'] <- 'Bracken Absent'
abs$Vegetation[abs$Vegetation == 'Bracken'] <- 'Bracken Present'
#create string for y axis label
yaxis_label <- paste("Absorbance at ", wavelength_of_interest, "nm")


#plot absorbance
abs_bxp <- ggboxplot(abs, x = "Habitat", y = "`Standardized Absorbance (L mg-1 cm-1)`", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = yaxis_label) + theme(
    #remove x axis label
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
    
  ) 

#refresh the dataframe d with the correct one
d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#plot nitrogen
tsn_bxp <- ggboxplot(d, x = "Habitat", y = 'TotalNitrogen', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(x = "Habitat x Vegetation",
       y = expression("Total Soil Nitrogen (g kg"^-1*")")) + theme(
         #remove x axis label, tickes, labels
         axis.title.x=element_blank(),
         # Remove panel border
         panel.border = element_blank(),  
         # Remove panel grid lines
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         # Remove panel background
         panel.background = element_blank(),
         # Add axis line
         axis.line = element_line(colour = "black", linewidth = 0.5),
         #change colour and thickness of axis ticks
         axis.ticks = element_line(colour = "black", linewidth = 0.5),
         #change axis labels colour
         axis.title.y = element_text(colour = "black"),
         #change tick labels colour
         axis.text.y = element_text(colour = "black"),
       ) 
#create the panel figure
supp_panel_bxp <- ggarrange(abs_bxp, alpha_bxp, tsn_bxp, 
                            labels = c("A", "B", "C"),
                            ncol = 2, nrow = 2,
                            common.legend = TRUE, legend="top")

show(supp_panel_bxp)
ggsave(path = "figures", paste0(Sys.Date(), "_Supp-2-Figure.svg"), supp_panel_bxp)

#### cations (Supp. Figure 3) ----
d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
)
#graphing boxplots. Units are ppm aka mg per kg
Al_bxp <- ggboxplot(d, x = "Habitat", y = "Al_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Al (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Zn_bxp <- ggboxplot(d, x = "Habitat", y = "Zn_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Zn (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Mn_bxp <- ggboxplot(d, x = "Habitat", y = "Mn_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Mn (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Mg_bxp <- ggboxplot(d, x = "Habitat", y = "Mg_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Mg (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

#graphing boxplots. Units are ppm aka mg per kg
Ca_bxp <- ggboxplot(d, x = "Habitat", y = "Ca_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Ca (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 


#graphing boxplots. Units are ppm aka mg per kg
Na_bxp <- ggboxplot(d, x = "Habitat", y = "Na_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("Na (mg g"^-1*")")) + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 



#graphing boxplots. Units are ppm aka mg per kg
K_bxp <- ggboxplot(d, x = "Habitat", y = "K_mgperg", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y =  expression("K (mg g"^-1*")"))+ theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

panel_bxp <- ggarrange(Na_bxp, K_bxp, Ca_bxp, Mg_bxp, Zn_bxp, Mn_bxp, Al_bxp,
                       labels = c("A", "B", "C", "D", "E", "F", "G"),
                       ncol = 2, nrow = 4,
                       common.legend = TRUE, legend="top")

show(panel_bxp)
ggsave(path = "figures", paste0(Sys.Date(), "_Supp-3-Figure.svg"), panel_bxp)

#### pH vs Al (Supp. Figure 4) ----
# Clean the data, removing non 0s
pHAlmodel <- na.omit(d)

                           #linear model of the relationship
model <- lm(pHAlmodel$Al_mgperg ~ pHAlmodel$pH)
summary(model)


#plot the relationship
plot <- ggplot(pHAlmodel, aes(x = pH, y = Al_mgperg)) +
  geom_point(color = "black", size = 2) + # Add points
  geom_smooth(method = "lm", color = "red", se = TRUE) + # Add regression line
  # Line of best fit
  annotate("text", x = 5.3, y = 0.2, 
           label = paste("y =", 
                         round(coef(model)[2], 1), "x +", 
                         round(coef(model)[1], 2), 
                         "\nR² =", 
                         round(summary(model)$r.squared, 2)), 
           hjust = 0, vjust = 1, size = 4, color = "black") + 
  labs(
    x = "pH",
    y = expression("[Al] (mg g"^-1*")")
  ) +  theme(
    
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 

show(plot)

ggsave(path = "Figures", paste0(Sys.Date(), "_Supp-Figure-4.svg"), plot)


####
#### morphotype abundance, evenness, shannon, simpson (Figure 5, Supp. Figure 5)----
d <- readr::read_csv(here::here("Data", "Week-1+2-mites-springtails-morphotypes.csv"), show_col_types = FALSE) 
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0

#rename Nonbracken so that it is plotted first, and bracken second
d$Vegetation[d$Vegetation == 'NonBracken'] <- 'Bracken Absent'
d$Vegetation[d$Vegetation == 'Bracken'] <- 'Bracken Present'
# 
# #just the morphotypes
# spe <- d[,-(1:3)]
# #summary of mesofauna proportions
# d <- d %>%
#   mutate(Mites = rowSums(select(spe, contains("Mite") | contains("Mstg") | contains("Noth") | contains("Astig") | contains("Box") | contains("Turt")| contains("Tank")| contains("Pter")| contains("Grp")) > 0))
# #number of collembola morphospecies
# d <- d %>%
#   mutate(Collembola = rowSums(select(spe, contains("Sprg") | contains("Neel") | contains("Smin") ) > 0))
# 
# ##analyse mites, springtails 
# 
# # Create a list of group codes (or patterns to match)
# codes <- c("GN", "GB", "UN", "UB", "WN", "WB")
# # Initialize an empty list to store results
# group_summary <- list()
# # Loop through each code and calculate mean and SD for rows where 'ID' contains the code
# for(code in codes) {
#   # Filter rows where ID contains the pattern (code)
#   group_data <- d %>% filter(grepl(code, `Sample ID`))
#   
#   # Calculate mean and standard deviation for the target_column
#   summary_stats <- group_data %>%
#     summarise(  #change first argument to column of choice
#       mean_value = mean(`Mites`, na.rm = TRUE),
#       sd_value = sd(`Mites`, na.rm = TRUE)
#     )
#   
#   # Store the result in the list
#   group_summary[[code]] <- summary_stats
# }
# # Combine the results into a single data frame
# df_summary <- bind_rows(group_summary, .id = "group_code")
# # View the result
# print(df_summary)
# 
# 
# 
# #plot number of individuals per 100g dry soil for each sample
# d_abundance <- d

d_abundance$Total <- rowSums(d_abundance[, (4:355)])
#if we assume each soil core was fully filled to 10 cm depth in a 5 cm diamter cylinder, we can standardize by soil volume and scale up to per m2
volume <- 0.1 * (pi*(0.025^2))
d_abundance$CoreVolume <- volume
#number of mesofauna per m3, to a depth of 10cm
d_abundance$individualsperm2to10cmdepth <- ((d_abundance$Total/d_abundance$CoreVolume)*0.1)/1000

# displays 100,000 as 100,000, not as 1e5
library(scales)

#plot carbon
abundance_bxp <- ggboxplot(d_abundance, x = "Habitat", y = 'individualsperm2to10cmdepth', color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = expression("1000 individuals per m"^2*" soil (to 10 cm depth)")) + theme( #remove x axis label
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) + scale_y_continuous(labels = label_comma())



#just the morphotypes
spe <- d[,-(1:3)]
#species richness
d$richness <- apply(spe[,]>0,1,sum)
#calculate diversity
d$shannon <- diversity(spe[,], "shannon")
d$simpson <- diversity(spe[,], "simpson")
d$evenness <- d$shannon / log(d$richness)

#graphing abundances
evenness_bxp <- ggboxplot(d, x = "Habitat", y = "evenness", color = "Vegetation", ylab = "Pielou's Evenness", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#boxplot of grassland vs heathland, species richness
richness_bxp <- ggboxplot(d, x = "Habitat", y = "richness", color = "Vegetation", ylab = "Morphospecies \n Richness", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#graphing boxplots with bracken split from nonbracken
shannon_bxp <-ggboxplot(d, x = "Habitat", y = "shannon", color = "Vegetation", ylab = "Shannon Diversity", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#graphing boxplots with bracken split from nonbracken
simpson_bxp <- ggboxplot(d, x = "Habitat", y = "simpson", color = "Vegetation", ylab = "Simpson Diversity", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15)) + theme(
  #remove x axis label
  axis.title.x=element_blank(),
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "black", linewidth = 0.5),
  #change colour and thickness of axis ticks
  axis.ticks = element_line(colour = "black", linewidth = 0.5),
  #change axis labels colour
  axis.title.y = element_text(colour = "black"),
  #change tick labels colour
  axis.text.y = element_text(colour = "black"),
) 
#
nem <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
)

#plot nematodes
nem_bxp <- ggboxplot(nem, x = "Habitat", y = "Nematodespergdrysoil", color = "Vegetation", palette = c("black", "limegreen"), lwd = 0.75, add = "jitter",  add.params = list(size = 1.5,alpha = 1, width = 0.15))  +
  labs(y = expression("Nematodes (individuals g dry soil "^-1*")"))  + theme(
    #remove x axis label, tickes, labels
    axis.title.x=element_blank(),
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.y = element_text(colour = "black"),
  ) 



all_bxp <- ggarrange(abundance_bxp, evenness_bxp, shannon_bxp, simpson_bxp, nem_bxp, 
                     labels = c("A", "B", "C", "D", "E"),
                     ncol = 2, nrow = 3,
                     common.legend = TRUE, legend="top")
show(all_bxp)
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), '_Figure-5.svg'), width = 9, height = 12, all_bxp)

ggsave(path = "Figures", paste0(Sys.Date(), '_Supp-Figure-5.svg'), width = 7, height = 4, richness_bxp)



#### NMDS of morphotype communities week 1 (Figure 6 LHS) ----
#csv containing only mites and springtails
d <- readr::read_csv(here::here("Data", "Week-1-mites-springtails-morphotypes.csv"), show_col_types = FALSE) 
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#remove all empty rows
d <- d[1:30,]
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0
#replace row index with sample names
rownames(d) <- d[,1]
#just the morphospecies counts
spe <- d[,-(1:3)]
spe <- as.matrix(spe)
#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#plot the NMDS
plot(example_NMDS, col = "white")
#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
#colors =c(rep("#44AA99",5),rep("#117733",5), rep("#88CCEE",5),rep("#332288",5), rep("#AA4499", 5), rep("#882255", 5)) 
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#shapes for point codes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5))
#display the stress if using all morphospecies
#text(-0.6,2.2, paste("Stress = ", round(example_NMDS$stress, 3)))
#display the stress if using only mites and springtails
text(-1.7,1.5, paste("Stress = ", round(example_NMDS$stress, 3)))
#visualise the points and ellipses
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment 
    
    #ellipses are of the kind se = standard error (or sd = standard deviation) at a 95% confidence
    ordiellipse(example_NMDS$point[grep(i,treat),], kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }


#save the file using Export -> Save As Image -> Width = 655, Height = 500 

#save the legend as its own plot, to be modified in Inkscape later
#plot the NMDS
plot(example_NMDS, col = "white")
#legend for text codes
#legend(-3.5,-1, legend = c("Grassland Bracken", "Grassland Non-bracken", "Heathland Bracken", "Heathland Non-bracken", "Woodland Bracken", "Woodland Non-Bracken"), fill = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), bty = "n")
#legend for point codes
legend(-1.7,0.9, legend=c("Grassland Bracken Present", "Grassland Bracken Absent", "Heathland Bracken Present", "Heathland Bracken Absent", "Woodland Bracken Present", "Woodland Bracken Absent"), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2"), pch = c(15, 0,16,1,17,2))

#ggsave(path = "Figures", paste0(Sys.Date(), '_morphospecies-nmds.svg'), last_plot(), width = 8, height = 6, dpi = 300)

# do PERMANOVA analysis


#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- d[,(2:3)]
#run the permanova
morph_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
morph_permanova


#run an ANOSIM. The ANOSIM test is similar to an ANOVA hypothesis test, but it uses a dissimilarity matrix as input instead of raw data. It is also non-parametric, meaning it doesn’t assume much about your data (like normal distribution etc), so it’s a good bet for often-skewed microbial abundance data. As a non-parametric test, ANOSIM uses ranked dissimilarities instead of actual distances, and in this way it’s a very nice complement to an NMDS plot. The main point of the ANOSIM test is to determine if the differences between two or more groups are significant.
#run an anosim - when grouping by habitat
ano = anosim(as.matrix(spe), grouping = d$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(spe), grouping = d$Vegetation, permutations = 9999, distance = "bray")
# When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. A Significance value less than 0.05 is generally considered to be statistically significant, and means the null hypothesis can be rejected. “The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups” (GUSTAME). In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition.
ano
plot(ano)

#### NMDS of morphotype communities week 2 (Figure 6 RHS) ----
#csv containing only mites and springtails

d <- readr::read_csv(here::here("Data", "Week-1+2-mites-springtails-morphotypes.csv"), show_col_types = FALSE) 
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#remove all empty rows
d <- d[1:30,]
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0
#replace row index with sample names
rownames(d) <- d[,1]
#just the morphospecies counts
spe <- d[,-(1:3)]
spe <- as.matrix(spe)

#all morphospecies
#spe <- all_data[,51:436]
#replace row index with sample names
rownames(spe) <- all_data[,1]
spe <- as.matrix(spe)

#k is the number of reduced dimensions
#trymax sets the default number of iterations
example_NMDS <- metaMDS(spe, distance = "bray", k = 2, maxit = 999, trymax = 500)
#Shephard plot shows scatter around the regession between the interpoint distances in the final configuration (i.e. the distances between each pair of communities) against their original dissimilarities.  Large scatter around the line suggests the original dissimilarities are not well preserved in the reduced number of dimensions
stressplot(example_NMDS)

#plot the NMDS
plot(example_NMDS, col = "white")


#assign the treatments to relevant rows of the dataframe
treat=c(rep("Grassland Bracken Present",5),rep("Grassland Bracken Absent",5), rep("Heathland Bracken Present",5),rep("Heathland Bracken Absent",5), rep("Woodland Bracken Present", 5), rep("Woodland Bracken Absent", 5))
#set the colour for each treatment
#colors =c(rep("#44AA99",5),rep("#117733",5), rep("#88CCEE",5),rep("#332288",5), rep("#AA4499", 5), rep("#882255", 5)) 
colors =c(rep("#999999",5),rep("#E69F00",5), rep("#56B4E9",5),rep("#009E73",5), rep("#CC79A7", 5), rep("#0072B2", 5)) 
#shapes for point codes
pchs<- c(rep(15, 5), rep(0, 5), rep(16, 5), rep(1, 5), rep(17, 5), rep(2, 5))
#display the stress for all morphotypes
#text(-0.8,1.4, paste("Stress = ", round(example_NMDS$stress, 3)))
#display the stress for only mites and springtails
text(-2,1.3, paste("Stress = ", round(example_NMDS$stress, 3)))
#visualise the points and ellipses
for(i in unique(treat)) {
  #we have added an if statement so we can chose which points and ellipses to plot at a time e.g. i == "Grassland Bracken".  If we want to plot all ellipses simultaneously, set i == i
  if(i == i){
    #plot the sample IDs on the NMDS, with the colour specific to the treatment
    # orditorp(example_NMDS$point[grep(i,treat),],display="sites", col=colors[grep(i,treat)], cex=0.7, air=0.01)
    #plot point codes for each site
    points(example_NMDS$point[grep(i,treat),], pch = pchs[grep(i,treat)], col = colors[grep(i,treat)], cex = 0.7)
    #plots ellipse with ellipse centered on the centroid of the samples from the same treatment (and thus encapsulating 95% of the variance)
    ordiellipse(example_NMDS$point[grep(i,treat),],kind = "se", conf = 0.95, draw="polygon",
                groups=treat[treat==i],col=colors[grep(i,treat)],label=F) } }



#save the file using Export -> Save As Image -> Width = 655, Height = 500 

# do PERMANOVA analysis
#data frame containing the independent variables (Habitat, Vegetation) we shall be using in our PERMANOVA
idvs <- d[,(2:3)]
#run the permanova
morph_permanova <- adonis2(spe ~ Habitat*Vegetation, idvs, permutations = 999, method = "bray", by = "terms")
morph_permanova


#run an ANOSIM. The ANOSIM test is similar to an ANOVA hypothesis test, but it uses a dissimilarity matrix as input instead of raw data. It is also non-parametric, meaning it doesn’t assume much about your data (like normal distribution etc), so it’s a good bet for often-skewed microbial abundance data. As a non-parametric test, ANOSIM uses ranked dissimilarities instead of actual distances, and in this way it’s a very nice complement to an NMDS plot. The main point of the ANOSIM test is to determine if the differences between two or more groups are significant.
#run an anosim - when grouping by habitat
ano = anosim(as.matrix(spe), grouping = d$Habitat, permutations = 9999, distance = "bray")
#check output of anosim
ano
plot(ano)
#run an anosim - when grouping by vegetation
ano = anosim(as.matrix(spe), grouping = d$Vegetation, permutations = 9999, distance = "bray")
# When interpreting these results you want to look at the ANOSIM statistic R and the Significance values. A Significance value less than 0.05 is generally considered to be statistically significant, and means the null hypothesis can be rejected. “The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups” (GUSTAME). In other words, the higher the R value, the more dissimilar your groups are in terms of microbial community composition.
ano
plot(ano)



#### nestedness and turnover analysis ----
d <- readr::read_csv(here::here("Data", "Week-1+2-mites-springtails-morphotypes.csv"), show_col_types = FALSE) 
#order samples by ID alphabetically
d <- arrange(d, d["Sample ID"])
d <- as.data.frame(d)
#remove all empty rows
d <- d[1:30,]
#replace null (empty excell cell) with "0"
d[is.na(d)] <- 0
#replace row index with sample names
rownames(d) <- d[,1]
#just the morphospecies counts
spe <- d[,-(1:3)]
spe <- as.matrix(spe)

#just the morphospecies counts
species_matrix <- spe
rownames(species_matrix) <-d[,1]
species_matrix <- as.matrix(species_matrix)
# Convert the matrix to a presence-absence (binary) matrix if not already
species_binary <- ifelse(species_matrix > 0, 1, 0)

#run nestedness/turnover analysis on each of the treatments
#grassland bracken present
gb <- species_binary[(1:5),]
#grassland bracken absent
gn <- species_binary[(6:10),]
#heathland bracken present
hb <- species_binary[(11:15),]
#heathland bracken absent
hn <- species_binary[(16:20),]
#woodland bracken present
wb <- species_binary[(21:25),]
#woodland bracken absent
wn <- species_binary[(26:30),]

#get betapart object
gb.core <- betapart.core(gb)
gn.core <- betapart.core(gn)
hb.core <- betapart.core(hb)
hn.core <- betapart.core(hn)
wb.core <- betapart.core(wb)
wn.core <- betapart.core(wn)

# multiple site measures
gb.multi <- beta.multi(gb.core, index.family = "jaccard")
gn.multi <- beta.multi(gn.core, index.family = "jaccard")
hb.multi <- beta.multi(hb.core, index.family = "jaccard")
hn.multi <- beta.multi(hn.core, index.family = "jaccard")
wb.multi <- beta.multi(wb.core, index.family = "jaccard")
wn.multi <- beta.multi(wn.core, index.family = "jaccard")

#compile the measures into a table. JTU = value of the turnover component, measured as turnover fraction of Jaccard dissimilarity.  JNE = value of the nestedness component, measured as nestedness-resultant faction of Jaccard dissimilarity.  .JAC = value of the overall beta diversity, measured as the Jaccard dissimilarity

#rownames
habitats <- c("Grassland", "Grassland", "Heathland", "Heathland", "Woodland", "Woodland")
bracken <- c("Present", "Absent","Present", "Absent","Present", "Absent")
#put all the dissimilarity measures into a table
# Combine lists into a data frame, each list becomes a row
beta_table <- as.data.frame(rbind(gb.multi, gn.multi, hb.multi, hn.multi, wb.multi, wn.multi))
#ensure all data is numeric
beta_table <- as.data.frame(lapply(beta_table, function(x) as.numeric(as.character(x))))
# Combine with new columns on the left
beta_table <- cbind(Habitat = habitats, Bracken = bracken, beta_table)
beta_table <- as.data.frame(beta_table)



#### model soil moisture against nematode abundances ----
d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 

#for modelling effects of moisture on nematode abundanve
library(Hmisc)
p <- ggplot(d,aes(`Water content as percentage wet soil mass`, Nematodespergdrysoil, colour = Vegetation)) +
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm')  +
  scale_colour_manual(name = "Vegetation", 
                      values = c("black", "limegreen")) + labs(x = "Water Content (%)", y = "Nematodes per g dry soil") +  theme(
                        
                        # Remove panel border
                        panel.border = element_blank(),  
                        # Remove panel grid lines
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        # Remove panel background
                        panel.background = element_blank(),
                        # Add axis line
                        axis.line = element_line(colour = "black", linewidth = 0.5),
                        #change colour and thickness of axis ticks
                        axis.ticks = element_line(colour = "black", linewidth = 0.5),
                        #change axis labels colour
                        axis.title.y = element_text(colour = "black"),
                        #change tick labels colour
                        axis.text.y = element_text(colour = "black"),
                      ) 

show(p)
#nested anova
anova <- aov(d$Nematodespergdrysoil ~ d$`Water content as percentage wet soil mass` + factor(d$Vegetation))
summary(anova)


#### ANOVA/stats tests ----
indvs <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 
#indvs <- na.omit(indvs)

anova <- aov(indvs$individualsperm2to10cmdepth ~ indvs$Habitat * indvs$Vegetation)
summary(anova)

#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
print(tukey)

#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

plot(anova, 1)
#levene test.  if p value < 0.05, there is eidence to suggest that the variance across groups is statistically significantly different.
leveneTest(indvs$individualsperm2to10cmdepth ~ indvs$Habitat * indvs$Vegetation)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residules
#extract the residuals
nestaov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = nestaov_residuals)





#### Moran's I ----
d <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 

#the coordinates of sample sites
coords <- d[c(4,5)]
dists <- as.matrix(coords)
#create spatial neighbours
nb <- knn2nb(knearneigh(coords, k = 4))
#nb <- dnearneigh(coords, d1 = 0, d2 = 30)  # adjust distance!
hist(dists)
#create spatial weights
lw <- nb2listw(nb, style = "W")
#test an environment variable
moran.test(d$individualsperm2to10cmdepth, lw)



moran.plot(d$Morphotype_shannon, lw)
#now try on the residuals from community dataset
#run NMDS first in NMDS tab (veg or morphotype), then:
site_scores <- scores(example_NMDS, display = "sites")  # 2 columns for 2 axes
head(site_scores)
#run Moran on NMDS axes
# NMDS axis 1
moran_axis1 <- moran.test(site_scores[,1], lw, zero.policy = TRUE)
moran_axis1

# NMDS axis 2
moran_axis2 <- moran.test(site_scores[,2], lw, zero.policy = TRUE)
moran_axis2

#d <- na.omit(d)

#the coordinates of sample sites
coords <- d[c(4,5)]
dists <- as.matrix(coords)
#create spatial neighbours
nb <- knn2nb(knearneigh(coords, k = 4))
#nb <- dnearneigh(coords, d1 = 0, d2 = 30)  # adjust distance!
hist(dists)
#create spatial weights
lw <- nb2listw(nb, style = "W")
#test an environment variable
moran.test(d$`K_mg/g`, lw)


#### Effect of moisture on inverts extracted ----
#extract data from master dataframe
moist <- readr::read_csv(
  here::here("Data", "Soil-parameters-measured.csv")
) 

#remove "non bracken" so that "bracken" comes alphabetically second, and so bracken box plots will be plotted to the right of non bracken box plots
moist$Vegetation[moist$Vegetation == 'Non-bracken'] <- 'Bracken Absent'
moist$Vegetation[moist$Vegetation == 'Bracken'] <- 'Bracken Present'


#week 1 extracts
wk1 <- readr::read_csv(
  here::here("Data", "Week-1-mites-springtails-morphotypes.csv"), show_col_types = FALSE
) 
#order samples by ID alphabetically
wk1 <- arrange(wk1, wk1["Sample ID"])
wk1 <- as.data.frame(wk1)
#remove all empty rows
wk1 <- wk1[1:30,]
#replace null (empty excell cell) with "0"
wk1[is.na(wk1)] <- 0
#replace row index with sample names
rownames(wk1) <- c("GB1", "GB2", "GB3", "GB4", "GB5", "GN10", "GN3", "GN4", "GN5","GN9", "UB1","UB10","UB4","UB6", "UB9", "UN1", "UN10", "UN3", "UN7", "UN8", "WB3", "WB4", "WB6", "WB7", "WB8", "WN2", "WN4", "WN5", "WN6", "WN8")
#make sure our variables are coded as factors
wk1$Habitat <- factor(wk1$Habitat, levels = c("Grassland", "Heathland", "Woodland"), labels = c("Grassland", "Heathland", "Woodland"))
#just the morphospecies counts
wk1spe <- wk1[,-(1:3)]
wk1spe <- as.matrix(wk1spe)
#number of morphospecies extracted after 7 days
#species richness
wk1$Richness_1 <- apply(wk1spe[,]>0,1,sum)
#number of individuals extracted after 7 days
wk1$Individuals_1 <- rowSums(wk1spe)


#csv containing week 1+2 extracts

wk2 <- readr::read_csv(
  here::here("Data", "Week-1+2-mites-springtails-morphotypes.csv"), show_col_types = FALSE
) 
#order samples by ID alphabetically
wk2 <- arrange(wk2, wk2["Sample ID"])
wk2 <- as.data.frame(wk2)
#remove all empty rows
wk2 <- wk2[1:30,]
#replace null (empty excell cell) with "0"
wk2[is.na(wk2)] <- 0
#replace row index with sample names
rownames(wk2) <- c("GB1", "GB2", "GB3", "GB4", "GB5", "GN10", "GN3", "GN4", "GN5","GN9", "UB1","UB10","UB4","UB6", "UB9", "UN1", "UN10", "UN3", "UN7", "UN8", "WB3", "WB4", "WB6", "WB7", "WB8", "WN2", "WN4", "WN5", "WN6", "WN8")
#make sure our variables are coded as factors
wk2$Habitat <- factor(wk2$Habitat, levels = c("Grassland", "Heathland", "Woodland"), labels = c("Grassland", "Heathland", "Woodland"))
#just the morphospecies counts
wk2spe <- wk2[,-(1:3)]
wk2spe <- as.matrix(wk2spe)
#number of morphospecies extracted after 7 days
#species richness
wk2$Richness_2 <-apply(wk2spe[,]>0,1,sum)
#number of individuals extracted after 7 days
wk2$Individuals_2 <- rowSums(wk2spe)


moist <- moist[c(1,2,3,23)]
#create dataframe containing moisture, and richness/individuals from 7 and 14 day extracts
moist <- merge(moist, wk1[, c("Sample ID", "Richness_1", "Individuals_1")], by = "Sample ID", all.x = TRUE)
#create dataframe containing moisture, and richness/individuals from 7 and 14 day extracts
moist <- merge(moist, wk2[, c("Sample ID", "Richness_2", "Individuals_2")], by = "Sample ID", all.x = TRUE)

#reshape the data into long format
# Example of reshaping (assuming your columns are named 'Richness_1' and 'Richness_2' for species richness at 1 and 2 weeks)
moist_long <- moist %>%
  pivot_longer(cols = c("Richness_1", "Richness_2", "Individuals_1", "Individuals_2"),
               names_to = c(".value", "Week"),
               names_pattern = "(.*)_(.*)")

names(moist_long)[4] <- "Soil Moisture"
#number of morphospecies and individuals gained from an additional week of extraction
moist$RichnessGained <- moist$Richness_2 - moist$Richness_1
moist$IndividualsGained <- moist$Individuals_2 - moist$Individuals_1



#see if moisture, x habitat, deretmines catch of mesofauna
#for modelling effects of moisture/extract wewek on mesofauna catch.richness
library(glmmTMB)
library(ggeffects)
model <- glmmTMB(Richness ~ Week*Habitat,
                 data = moist_long,
                 family = gaussian,
)
summary(model)

