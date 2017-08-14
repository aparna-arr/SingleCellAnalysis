library(fpc)
library(cluster)

widths<-NULL
for (i in 1:10) {
    print("On i",i)
    p<-pam(t(cd_mtec), i)
    s<-silhouette(p)
    widths<-c(widths,summary(s)$avg.width)
}
