# Importar las librerías necesarias
library(DynamicCancerDriverKM) # Libreria del paquete solicitado por el profesor 
library(e1071)  # Librería para SVM
library(caret) # Se utiliza para entrenar y evaluar modelos predictivos
library(dplyr) # Manipular y transformar datos 
library(pROC) # se utilizada para evaluar el rendimiento del modelo (curva (AUC) y curvas ROC)  
library(tidyverse) # manipulación y visualización de datos 

# Cargar los datos (BRCA_normal Y BRCA_PT) y unirlos en un solo dataset
merged_data <- rbind(DynamicCancerDriverKM::BRCA_normal, DynamicCancerDriverKM::BRCA_PT)


# Eliminar las columnas no relevantes para el análisis
columnas_no_relevantes <- c("barcode", "bcr_patient_barcode", "bcr_sample_barcode", "vital_status", "days_to_death", "treatments_radiation_treatment_or_therapy")
merged_data <- merged_data[, !(names(merged_data) %in% columnas_no_relevantes)]

# Verificar si hay datos nulos en el dataset
any(is.na(merged_data))



# VERIFICAR EL UMBRAL DE EXPRESION DE LAS MUESTRAS PARA CADA GEN 

# Obtener la matriz de muestras (sin la columna 'sample_type')
matriz_muestras <- as.matrix(merged_data[, -1])

# Calcular el umbral como el valor máximo en la matriz de muestras
umbral <- 0.0001 * max(matriz_muestras)

# Vector lógico que indica si cada muestra está activa para cada gen
muestras_activas <- matriz_muestras > umbral

# Contar la cantidad de TRUE para cada gen
verdaderos_por_gen <- colSums(muestras_activas)


# FILTRAR LOS GENES QUE SE EXPRESAN EN MAS DEL 20% DE LA MUESTRAS

# Calcular el umbral para conservar la columna con (20% del total de la muestra)
umbral_eliminar_columna <- nrow(matriz_muestras) * 0.2

# Encontrar las columnas a conservar (que tienen menos del 20% de TRUE)
columnas_a_conservar <- which(verdaderos_por_gen >= umbral_eliminar_columna)

# Filtrar el DataFrame original para conservar solo las columnas necesarias
filtered_data <- merged_data[, c(1, columnas_a_conservar + 1)] 


# CAMBIAR EL NOMBRE DE LAS COLUMNAS "Ensembl.ID" A "HGNC.symbol"

# Listar los nombres de los genes de "Ensembl.ID"
gen_names <- colnames(filtered_data)[-1]

# Burcar los nombres de los genes "Ensembl.ID"  con "HGNC.symbol"
genes_names <- AMCBGeneUtils::changeGeneId(gen_names, from = "Ensembl.ID")

# Cambiar los nombre el el DataFrame merged_filtrado de "Ensembl.ID" a "HGNC.symbol"                
colnames(filtered_data)[-1] <- genes_names$HGNC.symbol                 



# FILTRAR LOS GENES QUE ESTAN PRESENTES EN "filtered_data" y "la red PPI"


# Obtener los nombres de genes en PPI
genes_ppi <- DynamicCancerDriverKM::PPI$`Input-node Gene Symbol`

# Obtener los nombres de genes en filtered_data
genes_merged <- colnames(filtered_data)[-1]  # Excluir la columna "sample_type"

# Encontrar los genes comunes
genes_comunes <- intersect(genes_ppi, genes_merged)

# Filtrar el DataFrame PPI para incluir solo las filas con genes comunes
genes_comunes <- PPI[genes_ppi %in% genes_comunes, ]


# FRECUENCIA DE LOS GENES


# Contar la frecuencia de genes en 'Input-node Gene Symbol'
input_gene_counts <- genes_comunes %>%
  group_by(`Input-node Gene Symbol`) %>%
  summarize(input_count = n())

# Contar la frecuencia de genes en 'Output-node Gene Symbol'
output_gene_counts <- genes_comunes %>%
  group_by(`Output-node Gene Symbol`) %>%
  summarize(output_count = n())

# Unir las dos tablas por los nombres de los genes
merged_counts <- full_join(input_gene_counts, output_gene_counts, 
                           by = c("Input-node Gene Symbol" = "Output-node Gene Symbol"))

# Rellenar NA con 0
merged_counts[is.na(merged_counts)] <- 0

# Sumar las columnas 'input_count' y 'output_count' para obtener el total
merged_counts$total_count <- rowSums(merged_counts[, c("input_count", "output_count")])

# Ordenar por 'total_count' de mayor a menor
merged_counts <- merged_counts %>%
  arrange(desc(total_count))

# Filtrar los primeros 100 genes
top_genes <- head(merged_counts, 100)

# Obtener los nombres de los 100 genes de top_genes
top_100_genes <- top_genes$`Input-node Gene Symbol`


# Obtener la variable de respuesta 'y'
y <- filtered_data$sample_type

# Filtrar las columnas de tu dataframe 'data' para mantener solo los genes seleccionados
X <- filtered_data[, top_100_genes]

# Convertir la variable de respuesta a factor
y <- as.factor(y)

# Dividir los datos en conjuntos de entrenamiento y prueba
set.seed(123)
train_indices <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# Ajustar el modelo SVM
model <- svm(y_train ~ ., data = cbind(X_train, y_train), kernel = "linear")

# Realizar predicciones en el conjunto de prueba
predictions <- predict(model, newdata = cbind(X_test, y_test))

# Evaluar el rendimiento del modelo
confusionMatrix(predictions, y_test)

precision <- posPredValue(predictions, y_test)
roc_curve <- roc(y_test, as.numeric(predictions == "Tumor"), levels = c("Primary Tumor", "Solid Tissue Normal"))
roc_auc <- auc(roc_curve)
cat("Precision:", precision,  "AUC:", roc_auc)




  
  
  
  
  