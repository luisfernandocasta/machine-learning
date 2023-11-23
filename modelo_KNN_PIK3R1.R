# Importar las librerías necesarias
library(DynamicCancerDriverKM) # Libreria del paquete solicitado por el profesor
library(class)  # Librería para KNN
library(dplyr) # Manipular y transformar datos
library(caret) # se utiliza para entrenar y evaluar modelos de aprendizaje automático
library(ROCR) # Se utiliza para evaluar el rendimiento de modelos de clasificación
library(tidyverse) # manipulación y visualización de datos 

# Cargar los datos (BRCA_normal y BRCA_PT) y unirlos en un solo dataframe
merged_data <- rbind(DynamicCancerDriverKM::BRCA_normal, DynamicCancerDriverKM::BRCA_PT)

load("C:\\Users\\Santi\\Downloads\\geneScore.rdata")

# Eliminar las columnas no relevantes para el análisis
columnas_no_relevantes <- c("barcode", "bcr_patient_barcode", "bcr_sample_barcode", "vital_status", "days_to_death", "treatments_radiation_treatment_or_therapy")
merged_data <- merged_data[, !(names(merged_data) %in% columnas_no_relevantes)]

# Verificar si hay datos nulos en el dataset
any(is.na(merged_data))

# Calcular el porcentaje de expresión para cada gen
expression_percentage <- colMeans(merged_data > 0)
umbral <- 0.2
filtered_data <- merged_data %>%
  select(names(expression_percentage[expression_percentage >= umbral]))

# Ordenar el data frame por el campo 'score' de manera descendente
df_sorted <-  prub %>% arrange(desc(score))
# Seleccionar los primeros 100 genes con mayor score
top_genes <- df_sorted[1:100, ]


# Transponer top_genes para que los nombres de los genes estén en columnas
transposed_genes <- as.data.frame(t(top_genes))
# Asignar los nombres de las columnas usando la primera fila (antes de la transposición)
colnames(transposed_genes) <-  transposed_genes[1, ]
# Eliminar la primera fila, ya que ahora es la cabecera de las columnas
transposed_genes <- transposed_genes[-1, ]

# Obtener la variable de respuesta 'y'
y <- filtered_data$sample_type

# Filtrar las columnas de tu dataframe 'data' para mantener solo los genes seleccionados
X <- filtered_data[, colnames(transposed_genes)]

# Convertir la variable de respuesta a factor
y <- as.factor(y)

# Normalizar datos
X <- scale(X)

# Dividir los datos en conjuntos de entrenamiento y prueba con validación cruzada
set.seed(123)
indices <- createDataPartition(y, p = 0.7, list = FALSE)
train_indices <- indices
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# Ajustar el modelo KNN con validación cruzada
ctrl <- trainControl(method = "cv", number = 5)  # 5-fold cross-validation
k_values <- c(1, 3, 5, 7, 9)
tune_grid <- expand.grid(k = k_values)
model <- train(X_train, y_train, method = "knn", tuneGrid = tune_grid, trControl = ctrl)

# Realizar predicciones en el conjunto de prueba
predictions <- predict(model, newdata = X_test)

# Evaluar el rendimiento del modelo
confusionMatrix(predictions, y_test)

# Curva ROC y AUC
predictions <- as.numeric(predictions)
prediction_obj <- prediction(predictions, y_test)
perf <- performance(prediction_obj, "tpr", "fpr")
plot(perf, main = "ROC Curve", col = "blue", lwd = 2)

