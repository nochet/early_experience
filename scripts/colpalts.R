# Hex colorss
number_of_colors <- 8
pallete_names <- rownames(RColorBrewer::brewer.pal.info)
pallete_color_generator <- RColorBrewer::brewer.pal

# make color hexcodes:
palletes<-lapply(pallete_names,
                 pallete_color_generator,
                 n = number_of_colors) 

# turn into matrix:
palletes<-do.call(rbind,palletes)

rownames(palletes)<-pallete_names

# dump hard-coded R code to create the object to the console:
dump('palletes','')


# Colorblind with > 8
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)


colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pie(rep(1, 8), col = colorBlindBlack8)


# Colorblind with 9
palette.colors(palette = "Okabe-Ito")

# For sequential or diverging colormaps
library(rcartocolor)
display_carto_all(colorblind_friendly = TRUE)
