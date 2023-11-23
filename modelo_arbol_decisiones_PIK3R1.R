# Importar las librerías necesarias
library(DynamicCancerDriverKM) # Libreria del paquete solicitado por el profesor 
library(rpart)  # Librería para árboles de decisiones
library(caret) # Se utiliza para entrenar y evaluar modelos predictivos
library(dplyr) # Manipular y transformar datos 
library(pROC) # se utilizada para evaluar el rendimiento del modelo (curva (AUC) y curvas ROC)  
library(tidyverse) # manipulación y visualización de datos 

# Cargar datos
merged_data <- rbind(DynamicCancerDriverKM::BRCA_normal, DynamicCancerDriverKM::BRCA_PT)

load("C:\\Users\\Santi\\Downloads\\geneScore.rdata")

# Eliminar las columnas no relevantes
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

# Filtrar las columnas de 'filtered_data' para mantener solo los genes seleccionados
X <- filtered_data[, colnames(transposed_genes)]

# Convertir la variable de respuesta a factor
y <- as.factor(y)

# Definir la validación cruzada (k-fold cross-validation)
control <- trainControl(method = "cv", number = 5) # El número de folds (k) se puede cambiar según las necesidades

# Ajustar el modelo de árbol de decisiones con validación cruzada
model <- train(y ~ ., data = data.frame(X, y), method = "rpart", trControl = control)

# Realizar predicciones
predictions <- predict(model, newdata = data.frame(X), type = "raw")  

# Convertir a factor
predictions <- as.factor(predictions)

# Evaluar el rendimiento del modelo
confusionMatrix(predictions, y)

# Calcular la curva ROC y el área bajo la curva (AUC)
precision <- posPredValue(predictions, y)
roc_curve <- roc(y, as.numeric(predictions == "Tumor"), levels = c("Primary Tumor", "Solid Tissue Normal"))
roc_auc <- auc(roc_curve)
cat("Precision:", precision, "AUC:", roc_auc)
