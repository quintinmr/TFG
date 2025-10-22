library("readxl")
library("ggdendro")
library("cluster")
library("tidyverse")
library("factoextra")

# Datos
genotipos <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
A <- c(2.5, 5.1, 2.8, 3.0, 8.3, 4.6, 9.3, 5.3)
B <- c(2.7, 4.9, 3.2, 2.4, 3.3, 6.1, 4.6, 3.8)
datos <- data.frame(A, B, row.names = genotipos)

centroides <- data.frame(
  C1 = c(2.766666666, 2.7666666),
  C2 = c(8.80, 3.95),
  C3 = c(5.0, 4.9333)
)

distancia <- function(punto, centroide) {
  sqrt((punto[1] - centroide[1])^2 + (punto[2] - centroide[2])^2)
}

distancias <- t(apply(datos[, c("A", "B")], 1, function(punto) {
  c(
    dist_C1 = distancia(punto, centroides$C1),
    dist_C2 = distancia(punto, centroides$C2),
    dist_C3 = distancia(punto, centroides$C3)
  )
}))

distancias
# Aplicar k-means con k=3
set.seed(123)  # Para reproducibilidad
kmeans_result <- kmeans(datos, centers = 3)

# Añadir los clústeres al data frame
datos$cluster <- as.factor(kmeans_result$cluster)

library(ggplot2)

# Graficar los clusters
ggplot(datos, aes(x = A, y = B, color = cluster, label = rownames(datos))) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 4) +
  labs(title = "Clustering k-means (k=3)", x = "A", y = "B") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

df<-as.data.frame(scale(datos))
k2 <- kmeans(datos, centers = 3, nstart = 25)
fviz_cluster(k2,data=datos)

set.seed(123)
fviz_nbclust(datos, kmeans, method = "wss")


# Instala la librería si aún no la tienes
# install.packages("factoextra")

# Cargar la librería
library(factoextra)

# Crear el dataframe
genotipos <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")
A <- c(2.5, 5.1, 2.8, 3.0, 8.3, 4.6, 9.3, 5.3)
B <- c(2.7, 4.9, 3.2, 2.4, 3.3, 6.1, 4.6, 3.8)
datos <- data.frame(A, B, row.names = genotipos)

# Método del codo con k.max = 8 (máximo permitido)
fviz_nbclust(datos, kmeans, method = "wss", k.max = 7) +
  ggtitle("Método del Codo para el ejemplo del k-means")


fviz_nbclust(datos, kmeans, method = "silhouette", k.max = 7) +
  ggtitle("Método de la Silueta para el ejemplo del k-means")

