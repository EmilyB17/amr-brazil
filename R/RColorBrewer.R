### QUICK CODE FOR RCOLORBREWER

require(RColorBrewer)

# show all palettes
display.brewer.all(colorblindFriendly = TRUE)

# pick one palette
display.brewer.pal(n = 10, name = "Paired")

# return the hex code
hex <- brewer.pal(n = 10, name = "Paired")

# save the desired colors to your variables
farmcols <- c("farm1" = hex[1], "farm2" = hex[2])
sitecols <- c("feces" = hex[3], "rumen" = hex[4], "SNP" = hex[5])

# accent colors
accentcols <- c(hex[6:10])


