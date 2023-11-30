library(plotly)

data <- read.table("Spatial_distribution_1.txt",sep=";")
names(data) = c("indiv","x","y","z","y1","parent","diameter","species")

fig <- plot_ly(data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'lines',
               opacity = 1, line = list(width = 6, everscale = FALSE))
fig

#compare the distances, are they moved by the same amount each time?
d=rep(0,10)
for (i in 1:10){
  d[i] = dist(data[i:(i+1),2:4])
print(d[i])
}
plot(d,ylim=c(0,max(d))) #now identical

unique(data$diameter)