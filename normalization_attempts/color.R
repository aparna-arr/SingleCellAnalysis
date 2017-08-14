colrange <- 100

getcolor<-function(data) {
    colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
    cols <- colfunc(colrange)

    min_val <- min(data)
    max_val <- max(data)
    
    max_range <- max_val - min_val
    inc <- max_range / colrange

    bins<-ceiling((data - min_val) / inc)

    return(cols[bins])
}
