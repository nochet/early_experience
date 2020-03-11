# Functions for color mapping

# cbPalette - color blind-friendly
pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")

tx_fill_map <- function(){
  cmap <- scale_fill_manual(
    name = "Larval-Adult Diet",
    values = c("C_C" = pal[1],
               "C_DR" = pal[2],
               "DR_C" = pal[3],
               "DR_DR" = pal[4]))
  return(cmap)
}

tx_color_map <- function(){
  cmap <- scale_color_manual(
    name = "Larval-Adult Diet",
    values = c("C_C" = pal[1],
               "C_DR" = pal[2],
               "DR_C" = pal[3],
               "DR_DR" = pal[4]))
  return(cmap)
}
