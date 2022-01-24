
library(ggpubr)
library(ggplot2)
library(cowplot)
library(glmnet)
library(caret)
library(rms)
library(ROCR)
library(MASS)
library(pROC)


outDir <- "output"
if(!dir.exists(outDir)) dir.create(outDir, recursive = T)

### read-in data
dat.retro <- read.delim('data/retrospective_cohort.txt', stringsAsFactors = F, check.names = F)
dat.pro <- read.delim('data/prospective_cohort.txt', stringsAsFactors = F, check.names = F)


### tmp variable
dat.run <- dat.retro
dat.run$PEF <- NULL
x_valid <- dat.pro

### set random seed
set.seed(123)

### convert data type of category variables to factor
dat.run$site.of.onset <- factor(dat.run$site.of.onset)
dat.run$gender <- factor(dat.run$gender, levels = c(1, 2), labels = c('Male', 'Female'))
dat.run$diagnostic.level <- factor(dat.run$diagnostic.level)
dat.run$BP_cates <- factor(dat.run$BP_cates)

### split sample to training set and test set
train_ind <- sample(nrow(dat.run), floor(nrow(dat.run) * 0.7))
x_train <- dat.run[train_ind,]
y_train <- as.double(dat.run$response[train_ind])
x_test <- dat.run[-train_ind,]
y_test <- as.double(dat.run$response[-train_ind])


### build prediction model
full.model <- glm(response ~., data = x_train, family = binomial)
# Make predictions
probabilities <- full.model %>% predict(x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Prediction accuracy
observed.classes <- x_test$response
print(mean(predicted.classes == observed.classes))

### final model selected by backward searched with AIC
step.model <- full.model %>% stepAIC(direction = 'backward', trace = FALSE)
coef(step.model)

### summary model
model_summay <- as.data.frame(summary(step.model)$coefficients, stringsAsFactors = F)
write.table(data.frame(variable=rownames(model_summay ), model_summay ), col.names = T, row.names = F, sep = "\t", quote = F,
            file = paste0('output/PEF.final_model.summary.tsv'))

# Make predictions
probabilities <- predict(step.model, x_test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Prediction accuracy
observed.classes <- x_test$response
print(mean(predicted.classes == observed.classes))

pred <- prediction(probabilities, observed.classes)
perf <- performance(pred, 'tpr', 'fpr')

roc_test <- roc(x_test$response, probabilities, algorithm = 1)
print(roc_test$auc)

mycoords <- coords(roc_test, "all")
best.coords <- coords(roc_test, "best", best.method="youden")

#calculated CI, SE and Z value
auc_ci <- ci(roc_test)
se <- sqrt(var(roc_test))
z <- (roc_test$auc - 0.5) / se
auc.p <- 2 * pt(-abs(z), df = Inf)

### AUC plot
pdf(paste0('output/PEF.AUC.Test.pdf'), width = 6.5, height = 6.5)
plot(perf, legacy.axes = TRUE, predictor = 'probabilities',
     main="ROC Curve (Test set n=95)", col="#1c61b6", xlim=c(0,1), ylim=c(0,1), extends = c(0,0))
abline(a = 0, b = 1, lty=2)
legend("bottomright", legend=c(paste0('AUC = ', signif(roc_test$auc, digits = 3))),
       col=c("#1c61b6"), lwd=2)
dev.off()

###### valid data
x_valid_1 <- x_valid[, colnames(x_test)]
x_valid_1$gender <- factor(x_valid_1$gender, levels = c(1,2), labels = c('Male', 'Female'))
x_valid_1$site.of.onset <- factor(x_valid_1$site.of.onset)
x_valid_1$BP_cates <- factor(x_valid_1$BP_cates)
x_valid_1$diagnostic.level <- factor(x_valid_1$diagnostic.level)
probabilities <- predict(step.model, x_valid_1, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Prediction accuracy
observed.classes <- x_valid_1$response
print(mean(predicted.classes == observed.classes))

pred <- prediction(probabilities, observed.classes)
perf <- performance(pred, 'tpr', 'fpr')

roc_test <- roc(x_valid_1$response, probabilities, algorithm = 1)

mycoords <- coords(roc_test, "all")
best.coords <- coords(roc_test, "best", best.method="youden")

auc_ci <- ci(roc_test)

se <- sqrt(var(roc_test))
z <- (roc_test$auc - 0.5) / se
auc.p <- 2 * pt(-abs(z), df = Inf)


print(roc_test$auc)
pdf(paste0('output/PEF.AUC.Validation.pdf'), width = 6.5, height = 6.5)
plot(perf, legacy.axes = TRUE, predictor = 'probabilities',
     main="ROC Curve (prospective data set n=97)", col="#1c61b6")
abline(a = 0, b = 1, lty=2)
legend("bottomright", legend=c(paste0('AUC = ', signif(roc_test$auc, digits = 3))),
       col=c("#1c61b6"), lwd=2)
dev.off()


### calibrate plot
coef(step.model)
model.1 <- lrm(response ~ gender + site.of.onset + site.of.onset + onset.age + CREA + `HDL-C` + 
                 BMI + PLT +`BASO#` + URIC + TG + CHOL, data = x_train, x = T, y = T)
cali <- calibrate(model.1, method = 'boot', B = 1000)
pdf(file = paste0('output/PEF.calibrate_curve.pdf'), width = 6.5, height = 6.5)
plot(cali, xlim=c(0,1), ylim = c(0,1))
dev.off()



####### heatmap plot
require(ComplexHeatmap)
df_plot <- dat.retro
df_plot <- df_plot[order(df_plot$PEF),]
df_plot$gender <- factor(df_plot$gender, levels = c(1, 2), labels = c('Male', 'Female'))
df_plot$site.of.onset <- factor(df_plot$site.of.onset)
ha <- HeatmapAnnotation(
  PEF = anno_barplot(df_plot$PEF, gp = gpar(fill = '#D3D3D3', col = '#D3D3D3'), 
                     axis_param = list(side = 'right',
                                       at = c(0, 5, 12),
                                       labels = c('0', '5', '12')),
                     annotation_name_gp= gpar(fontsize = 20),
                     height = unit(4, "cm"), border = FALSE),
  annotation_name_side = "right"
) %v% 
  HeatmapAnnotation(
    BMI = df_plot$BMI,
    gender = df_plot$gender,
    URIC = df_plot$URIC,
    onset.age = df_plot$onset.age,
    CREA = df_plot$CREA,
    'HDL-C' = df_plot$`HDL-C`,
    site.of.onset = df_plot$site.of.onset,
    PLT = df_plot$PLT,
    'BASO#' = df_plot$`BASO#`,
    TG = df_plot$TG,
    CHOL = df_plot$CHOL,
    annotation_name_side = "right",
    #annotation_label = c(rep(" ", 5)),
    col = list(
      gender = c('Male' = '#B7D19C', 'Female' = '#BF82CA'),
      site.of.onset = c('3' = '#e6194b', "1" = 'steelblue', '2' = 'lightsalmon')
      
    )
  )
pdf('output/heatmap.PEF.pdf',
    width = 11, height = 5.5, onefile = FALSE)
draw(ha,heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()








