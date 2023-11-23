# Importar las librerías necesarias
library(DynamicCancerDriverKM) # Libreria del paquete solicitado por el profesor 
library(glmnet) # Se utiliza para ajustar modelos de regresión lineal generalizada con regularización tipo L1 (lasso) y L2 (ridge
library(caret) # Se utiliza para entrenar y evaluar modelos predictivos
library(dplyr) # Manipular y transformar datos 
library(pROC) # se utilizada para evaluar el rendimiento del modelo (curva (AUC) y curvas ROC)  
library(tidyverse) # manipulación y visualización de datos 

# Cargar los datos (BRCA_normal Y BRCA_PT) y unirlos en un solo dataset
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

# Filtrar las columnas del dataframe para mantener solo los genes seleccionados
X <- filtered_data[, colnames(transposed_genes)]

# Convertir la variable de respuesta a factor
y <- as.factor(y)


# Dividir los datos en conjuntos de entrenamiento y prueba
set.seed(123)
train_indices <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# Ajustar un modelo de regresión logística con regularización L1 (lasso)
model <- cv.glmnet(as.matrix(X_train), y_train, family = "binomial", alpha = 1)

# Realizar predicciones en el conjunto de prueba
predictions <- predict(model, newx = as.matrix(X_test), s = "lambda.min", type = "response")

# Convertir las probabilidades a clases
predictions <- as.factor(ifelse(predictions > 0.5, levels(y)[2], levels(y)[1]))

# Evaluar el rendimiento del modelo
confusionMatrix(predictions, y_test)

precision <- posPredValue(predictions, y_test)
roc_curve <- roc(y_test, as.numeric(predictions == "Tumor"), levels = c("Primary Tumor", "Solid Tissue Normal"))
roc_auc <- auc(roc_curve)
cat("Precision:", precision,  "AUC:", roc_auc)
