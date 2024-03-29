# Functions for color mapping

pal <- c("#999999", "#E69F00")
#pal <- c("#999999", "#E69F00", "#56B4E9")

tx_fill_map <- function(){
  cmap <- scale_fill_OkabeIto(
    name = "Diet",
    values = c("C" = pal[1],
               "DR" = pal[2]))
               #"HS" = pal[3]))
  return(cmap)
}

tx_color_map <- function(){
  cmap <- scale_color_OkabeIto(
    name = "Diet",
    values = c("C" = pal[1],
               "DR" = pal[2]))
               #"HS" = pal[3]))
  return(cmap)
}
