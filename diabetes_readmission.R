# readmission prediction in diabetes patients

# The dataset represents 10 years (1999-2008) of clinical care at 130 US hospitals 
# and integrated delivery networks. It includes over 50 features representing 
# patient and hospital outcomes.
# 
# Beata Strack, Jonathan P. DeShazo, Chris Gennings, Juan L. Olmo, 
# Sebastian Ventura, Krzysztof J. Cios, and John N. Clore, 
# ???Impact of HbA1c Measurement on Hospital Readmission Rates: Analysis of 70,000 
# Clinical Database Patient Records,??? BioMed Research International, 
# vol. 2014, Article ID 781670, 11 pages, 2014. 
# 
# https://archive.ics.uci.edu/ml/datasets/Diabetes+130-US+hospitals+for+years+1999-2008

library(dplyr)
library(GGally)
library(ggplot2)
library(corrplot)
library(psych)
library(caret)
library(rpart)
library(randomForest)
library(nnet)
library(e1071)
library(ROCR)
library(CORElearn)
library(lasso2)
#library(pscl)

# read in data
filename <- '/Users/1Air/documents/614/dataset_diabetes/diabetic_data.csv'
data <- read.table(filename, sep = ",", header = T, na.strings = "?")
head(data)

#load(file = "data2.rdata")

#########################################################

# PREPROCESSING, CLEANING

# remove unnecessary columns patient/encounter IDs, weight(97% missing), drugs
data <- select(data, -encounter_id, -patient_nbr, -weight,-(25:41),-(43:47))

# exploratory analysis and plots
summary(data)
ggpairs(data)
c <- cor(data, use= "pairwise.complete.obs")
corrplot(c)
# time-in-hospital is positively correlated with number of lab procedures,
# number of non-lab procedures, number of medications and number of diagnoses
# number of emergency visits correlates with number of inpatient visits

# fix some missing values
data$race[is.na(data$race)] <- "Other"
any(is.na(data$race)) # false

# PLOTS

# variable distributions
plot(data$age, main = "age distribution") # age: mode 70-80yrs normal distribution, right skewed
plot(data$gender, main = "gender distribution") # gender: female 53% male 47%
plot(data$A1Cresult, main = "A1C") # A1Cresult: 84% no A1c results, 8% >8
plot(data$readmitted, main = "readmissions") # readmission: >50% no readmission
plot(data2$admission_source, main = "admission source") # emergency 60%
plot(data2$discharged_to, main = "readmissions") # transferred to another facility 70%
# race: 75% caucasian
# admission source: emergency >50%
# time in hospital: mode 3 days
# max_glu_serum: none in >90%

g <- ggplot(data2, aes(x=age, y=time_in_hospital))
g + geom_boxplot(aes(fill=readmitted))
# patients with <30 day readmissions in their 70s-80s had longer time in hospital
# patients in their 30s-40s with <30 day readmission spent longer in the hospital

g <- ggplot(data2,aes(x=A1Cresult, y=num_medications))
g + geom_boxplot(aes(color=A1Cresult)) 
# not much difference in distribution across groups

g <- ggplot(data2,aes(x=A1Cresult, y=time_in_hospital))
g + geom_boxplot(aes(fill=diabetesMed)) + facet_grid(. ~ readmitted)
# patients with no readmission had generally had less time in hospital
# for those not taking diabetes medication

g <- ggplot(data2,aes(x=age, y=num_medications))
g + geom_boxplot(aes(fill=age))
# number of medications was highest in 60-70yr olds

g <- ggplot(data2,aes(x=diag2, y=time_in_hospital))
g + geom_boxplot(aes(fill=diag2))
# respiratory and injury diagnosis 2 stayed longer in hospital


#################################################################

# FEATURE EXTRACTION

data2 <- data

data2$diag_1 <- as.numeric(levels(data2$diag_1)[data2$diag_1])
data2$diag_2 <- as.numeric(levels(data2$diag_2)[data2$diag_2])
data2$diag_3 <- as.numeric(levels(data2$diag_3)[data2$diag_3])

# diagnosis1
data2$diagnosis_group <- factor( rep("other",nrow(data2)),ordered = F, 
                                 levels = c("circulatory","respiratory","Digestive","Diabetes","Injury",
                                            "Musculoskeletal","Genitourinary","Neoplasms","other"))
data2$diagnosis_group[data2$diag_1>=390 & data2$diag_1 <= 459 | data2$diag_1==785] <- "circulatory"
data2$diagnosis_group[data2$diag_1>=460 & data2$diag_1 <= 519 | data2$diag_1==786] <- "respiratory"
data2$diagnosis_group[data2$diag_1>=520 & data2$diag_1 <= 579 | data2$diag_1==787] <- "Digestive"
data2$diagnosis_group[data2$diag_1>=250 & data2$diag_1 < 251] <- "Diabetes"
data2$diagnosis_group[data2$diag_1>800 & data2$diag_1 <= 999] <- "Injury"
data2$diagnosis_group[data2$diag_1>=710 & data2$diag_1 <= 739] <- "Musculoskeletal"
data2$diagnosis_group[data2$diag_1>=580 & data2$diag_1 <= 629 | data2$diag_1==788] <- "Genitourinary"
data2$diagnosis_group[data2$diag_1>=140 & data2$diag_1 <= 239 | data2$diag_1>=790 & 
                        data2$diag_1 <= 799 | data2$diag_1==780 | data2$diag_1>=240 & data2$diag_1 < 250 |
                        data2$diag_1>=251 & data2$diag_1 <= 279 | data2$diag_1>=680 & data2$diag_1 <= 709 |
                        data2$diag_1>=001 & data2$diag_1 <= 139 | data2$diag_1==781 |
                      data2$diag_1==782 | data2$diag_1==784] <- "Neoplasms"

# diagnosis_2
data2$diagnosis_2 <- factor( rep("other",nrow(data2)),ordered = F, 
                                 levels = c("circulatory","respiratory","Digestive","Diabetes","Injury",
                                            "Musculoskeletal","Genitourinary","Neoplasms","other"))
data2$diagnosis_2[data2$diag_2>=390 & data2$diag_2 <= 459 | data2$diag_2==785] <- "circulatory"
data2$diagnosis_2[data2$diag_2>=460 & data2$diag_2 <= 519 | data2$diag_2==786] <- "respiratory"
data2$diagnosis_2[data2$diag_2>=520 & data2$diag_2 <= 579 | data2$diag_2==787] <- "Digestive"
data2$diagnosis_2[data2$diag_2>=250 & data2$diag_2 < 251] <- "Diabetes"
data2$diagnosis_2[data2$diag_2>800 & data2$diag_2 <= 999] <- "Injury"
data2$diagnosis_2[data2$diag_2>=710 & data2$diag_2 <= 739] <- "Musculoskeletal"
data2$diagnosis_2[data2$diag_2>=580 & data2$diag_2 <= 629 | data2$diag_2==788] <- "Genitourinary"
data2$diagnosis_2[data2$diag_2>=140 & data2$diag_2 <= 239 | data2$diag_2>=790 & 
                        data2$diag_2 <= 799 | data2$diag_2==780 | data2$diag_2>=240 & data2$diag_2 < 250 |
                        data2$diag_2>=251 & data2$diag_2 <= 279 | data2$diag_2>=680 & data2$diag_2 <= 709 |
                        data2$diag_2>=001 & data2$diag_2 <= 139 | data2$diag_2==781 |
                        data2$diag_2==782 | data2$diag_2==784] <- "Neoplasms"

# diagnosis_3
data2$diagnosis_3 <- factor( rep("other",nrow(data2)),ordered = F, 
                                 levels = c("circulatory","respiratory","Digestive","Diabetes","Injury",
                                            "Musculoskeletal","Genitourinary","Neoplasms","other"))
data2$diagnosis_3[data2$diag_3>=390 & data2$diag_3 <= 459 | data2$diag_3==785] <- "circulatory"
data2$diagnosis_3[data2$diag_3>=460 & data2$diag_3 <= 519 | data2$diag_3==786] <- "respiratory"
data2$diagnosis_3[data2$diag_3>=520 & data2$diag_3 <= 579 | data2$diag_3==787] <- "Digestive"
data2$diagnosis_3[data2$diag_3>=250 & data2$diag_3 < 251] <- "Diabetes"
data2$diagnosis_3[data2$diag_3>800 & data2$diag_3 <= 999] <- "Injury"
data2$diagnosis_3[data2$diag_3>=710 & data2$diag_3 <= 739] <- "Musculoskeletal"
data2$diagnosis_3[data2$diag_3>=580 & data2$diag_3 <= 629 | data2$diag_3==788] <- "Genitourinary"
data2$diagnosis_3[data2$diag_3>=140 & data2$diag_3 <= 239 | data2$diag_3>=790 & 
                        data2$diag_3 <= 799 | data2$diag_3==780 | data2$diag_3>=240 & data2$diag_3 < 250 |
                        data2$diag_3>=251 & data2$diag_3 <= 279 | data2$diag_3>=680 & data2$diag_3 <= 709 |
                        data2$diag_3>=001 & data2$diag_3 <= 139 | data2$diag_3==781 |
                        data2$diag_3==782 | data2$diag_3==784] <- "Neoplasms"
# admission_source
data2$admission_source <- factor( rep("other",nrow(data2)),ordered = F, 
                             levels = c("clinic_referral", "emergency","other"))
data2$admission_source[data2$admission_source_id==c(1,2,3)]<- "clinic_referral"
data2$admission_source[data2$admission_source_id==7]<- "emergency"

# discharged_to
data2$discharged_to <- factor( rep("transferred",nrow(data2)),ordered = F, 
                                  levels = c("home", "transferred","left_AMA"))
data2$discharged_to[data2$discharge_disposition_id==c(1,6,8)]<- "home"
data2$discharged_to[data2$discharge_disposition_id==7]<- "left_AMA"

data2 <- select(data2, -diag_1, -diag_2, -diag_3, -admission_type_id, -discharge_disposition_id)
data2 <- select(data2, -medical_specialty)
data2 <- rename(data2, diag1 = diagnosis_group, diag2=diagnosis_2, diag3 = diagnosis_3)

# payer_code
data2$payer_code2 <- factor( rep("other",nrow(data2)),ordered = F, 
                               levels = c("other", "self_pay"))
data2$payer_code2[data2$payer_code=="SP"]<- "self_pay"
data2 <- select(data2, -payer_code)
data2 <- select(data2, -admission_source_id)
data2 <- rename(data2, payer_code=payer_code2)

#######################################################

# QUICK PCA with numeric variables

y <- select(data2, readmitted)
X <- select(data2, time_in_hospital, num_lab_procedures, num_procedures, num_medications, 
            number_outpatient, number_emergency, number_inpatient, number_diagnoses)
# no rotation
pca_noRot <- principal(X, nfactors = 5, rotate = "none")
rotation2_noRot <- data.frame(cbind(pca_noRot$score, y))
head(rotation2_noRot)
pca_noRot$loadings
# cumulative 77%
# pc1 number of medications and time in hosptial
# pc2 number of in-patient visits and emergency
# pc3 number of procedures
# pc4 number of out-patient visits
# pc5 number of diagnoses

# linear model of class as a function of PCs
linModel_noRot <- glm(readmitted ~ PC1 + PC2 + PC3 + PC4 + PC5, data = rotation2_noRot, family = binomial)
summary(linModel_noRot)
# all PCs are significant ***

# PCA with varimax rotation
pca2 <- principal(X, nfactors = 5, rotate = "varimax") 
rotation2 <- data.frame(cbind(pca2$score, y))
pca2$loadings
# RC1 lab_procedures, time in hospital
# RC3 num_procedures, medications
# RC2 emergency visits, in-patient
# RC5 number diagnoses
# RC4 number_outpatient

plot(pca2)
summary(rotation2)
plot(rotation2)

###########################################################

# SPLIT DATA INTO TRAINING AND TESTING SET

set.seed(123)
inTrain <- createDataPartition(y = data2$readmitted, p = .66,list = FALSE)
train <- data2[ inTrain,]
test <- data2[-inTrain,]
nrow(train) # 67167
nrow(test) # 3459


########################################################################

# LOGISTIC REGRESSION

fit_all <- glm(readmitted ~., data=train, family=binomial)
summary(fit_all)
# significant predictors     p-value
#______________________________________
# discharged_to transferred 4.69e-15
# discharged_toleft_AMA     0.000358
# admission_sourceother     0.029632
# diag3Genitourinary        0.000647
# diag3Diabetes             0.002519
# diag2Neoplasms            0.008182
# diag2Diabetes             0.000185
# diag2respiratory          0.005380
# diag1Neoplasms            3.12e-05
# diag1Genitourinary        0.002069
# diag1Musculoskeletal      0.028831
# diag1Digestive            0.000808
# diag1respiratory          4.13e-14
# diabetesMedYes            1.68e-06
# changeNo                  0.001678
# insulinUp                 0.027916
# insulinSteady             7.34e-07
# insulinNo                 3.02e-07
# age[20-30)                0.015920  
# age[30-40)                0.015453  
# age[40-50)                0.018261  
# age[50-60)                0.022707  
# age[60-70)                0.010500  
# age[70-80)                0.007974 
# age[80-90)                0.007353 
# age[90-100)               0.011584
# time_in_hospital          7.34e-07
# num_procedures            1.25e-06
# num_medications           3.91e-05
# number_emergency          7.90e-05
# number_inpatient          < 2e-16 
# number_diagnoses          1.79e-12


# pseudo R-squared for logistic regression model
logisticPseudoR2s <- function(LogModel) {
  dev <- LogModel$deviance 
  nullDev <- LogModel$null.deviance 
  modelN <-  length(LogModel$fitted.values)
  R.l <-  1 -  dev / nullDev
  R.cs <- 1- exp ( -(nullDev - dev) / modelN)
  R.n <- R.cs / ( 1 - ( exp (-(nullDev / modelN))))
  cat("Pseudo R^2 for logistic regression\n")
  cat("Hosmer and Lemeshow R^2  ", round(R.l, 3), "\n")
  cat("Cox and Snell R^2        ", round(R.cs, 3), "\n")
  cat("Nagelkerke R^2           ", round(R.n, 3),    "\n")
}

logisticPseudoR2s(fit_all)
#Pseudo R^2 for logistic regression
#Hosmer and Lemeshow R^2   0.038 
#Cox and Snell R^2         0.027 
#Nagelkerke R^2            0.053 

# not very impressive predictor


# main effects, with A1C result
fit_a1c <- glm(readmitted ~ race+age+discharged_to+time_in_hospital+
             num_lab_procedures+num_procedures+num_medications+number_outpatient+
             number_emergency+number_inpatient+number_diagnoses+
             insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult, 
           data=train, family = binomial)
summary(fit_a1c)
# results not very different from fit_all
logisticPseudoR2s(fit_a1c)
pR2(fit_a1c)
# adjusted R-squared mostly same as fit_all
anova(fit_a1c, test="Chisq")

# predictors with largest effect on Residual deviance
#time_in_hospital    1   103.55     67148      46825 < 2.2e-16 ***
#number_inpatient    1  1036.04     67142      45538 < 2.2e-16 ***
#number_emergency    1   147.03     67143      46574 < 2.2e-16 ***
#age                 9    82.35     67151      46991 5.506e-14 ***
#discharged_to       2    62.78     67149      46929 2.332e-14 ***



######################################################################

# RPART DECISION TREES

rpart_tree <- rpart(formula = readmitted ~ age+discharged_to+time_in_hospital+
                      num_lab_procedures+num_procedures+num_medications+number_outpatient+
                      number_emergency+number_inpatient+number_diagnoses+
                      insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult, 
                    data=train, method = 'class')
summary(rpart_tree)

#Variable importance
#number_inpatient  number_emergency number_outpatient 
#     93                 6                 1 

# prediction

test$pred_readmit <- predict(rpart_tree, test, type="class")
table(predict(rpart_tree, test, type="class"), test$readmitted)
prop.table(table(test$readmitted, test$pred_readmit),1)

#         <30        >30         NO
# <30 0.00000000 0.28288617 0.71711383
# >30 0.00000000 0.19943727 0.80056273
# NO  0.00000000 0.08020131 0.91979869

# predicts no readmission with 92% accuracy

confusionMatrix(test$pred_readmit, test$readmitted)
# Accuracy : 0.5662  
#                   Class: <30 Class: >30 Class: NO
#Sensitivity               0.000    0.19944    0.9198
#Specificity               1.000    0.88524    0.2196

######################################################################

# RANDOM FOREST
Rf_fit<-randomForest(formula=readmitted ~ age+discharged_to+time_in_hospital+
                       num_lab_procedures+num_procedures+num_medications+number_outpatient+
                       number_emergency+number_inpatient+number_diagnoses+
                       insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult,
                     data=train)
print(Rf_fit)

test$pred_readmit <- predict(Rf_fit, test, type = "response")
table(test$readmitted, test$pred_readmit)
prop.table(table(test$readmitted, test$pred_readmit),1)

        #<30         >30          NO
#<30 0.013986014 0.374255374 0.611758612
#>30 0.003227141 0.341001241 0.655771618
#NO  0.001072214 0.163137297 0.835790489

#predicts no readmission with 84% accuracy

importance(Rf_fit)
#                       MeanDecreaseGini
#num_lab_procedures        5454.0726
#num_medications           4582.3072
#time_in_hospital          3278.6999
#age                       2731.0087
#diag3                     3147.4519
#diag2                     3052.8240
#diag1                     2884.9768
#number_diagnoses          2025.1397
#num_procedures            2021.3074
#number_inpatient          1768.9447
#insulin                   1440.8706
#A1Cresult                  963.3717
#number_outpatient          844.8839
#discharged_to              776.9524
#number_emergency           650.8291
#change                     606.1124
#diabetesMed                362.3747



##############################################

# NEURAL NET

nnet_model <- nnet(formula = readmitted ~ race+age+discharged_to+time_in_hospital+
                     num_lab_procedures+num_procedures+num_medications+number_outpatient+
                     number_emergency+number_inpatient+number_diagnoses+
                     insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult, 
                   data=train, size = 10, maxit = 100)

test$pred_readmit <- predict(nnet_model, test, type = "class")
prop.table(table(test$readmitted, test$pred_readmit),1)
# 10 hidden nodes, max iterations 100
#_________________________________________
# correctly predicts NO readmission 93%
# correctly predicts >30 readmission 17%
# correctly predicts <30 readmission 2%

# 10 hidden nodes, max iterations 200
#_________________________________________
# correctly predicts NO readmission 84%
# correctly predicts >30 readmission 34%
# correctly predicts <30 readmission 0.2%

confusionMatrix(test$pred_readmit, test$readmitted)
#Overall Accuracy : 0.5669  
#                       Class: <30 Class: >30 Class: NO
#Sensitivity            0.027872    0.17370    0.9321
#Specificity            0.995319    0.90980    0.1954

##############################################

# SUPPORT VECTOR MACHINES

SVMmodel <- svm(readmitted ~ age+discharged_to+time_in_hospital+
                  num_lab_procedures+num_procedures+num_medications+number_outpatient+
                  number_emergency+number_inpatient+number_diagnoses+
                  insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult,
                data=train, kernel = "linear")
                  #kernel = "rbf", gamma = 0.1, cost = 1)
print(SVMmodel)
summary(SVMmodel)
x <- select(test, -readmitted)
y <- select(test, readmitted)
pred <- predict(SVMmodel, x)
test$pred_readmit <- pred
prop.table(table(test$readmitted, test$pred_readmit),1)
confusionMatrix(test$pred_readmit, test$readmitted)

# visualize (classes by color, SV by crosses):
plot(cmdscale(dist(x)),
     col = as.integer(y),
     pch = c("o","+")[1:(nrow(data2)) %in% SVMmodel$index + 1])


##############################################

# KNN
# build kNN kernel regressor by preventing tree splitting

modelKernel <- CoreModel(readmitted ~ age+discharged_to+time_in_hospital+
                           num_lab_procedures+num_procedures+num_medications+number_outpatient+
                           number_emergency+number_inpatient+number_diagnoses+
                           insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult, 
                         train, model="regTree", modelTypeReg=6)
print(modelKernel)

# prediction on test set
pred <- predict(knnModel, test, type="class")
mEval <- modelEval(knnModel, test$readmitted, pred)
print(mEval)

test$pred_readmit <- pred
prop.table(table(test$readmitted, test$pred_readmit),1)
confusionMatrix(test$pred_readmit, test$readmitted)
# Overall Accuracy : 0.5674  
# correctly predicts NO readmission 93%
# correctly predicts >30 readmission 18%
# correctly predicts <30 readmission 0%

#                       Class: <30 Class: >30 Class: NO
#Sensitivity              0.0000    0.18577    0.9321
#Specificity              1.0000    0.89971    0.2030

##############################################

# KNN
knnModel <- CoreModel(readmitted ~ age+discharged_to+time_in_hospital+
                        num_lab_procedures+num_procedures+num_medications+number_outpatient+
                        number_emergency+number_inpatient+number_diagnoses+
                        insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult,
                      data=train, model="knn", kInNN = 10)
print(knnModel)

# prediction on test set
pred <- predict(knnModel, test, type="class")
mEval <- modelEval(knnModel, test$readmitted, pred)
print(mEval) # evaluation of the model
#$accuracy 0.5122691
#$AUC 0.524898

test$pred_readmit <- pred
prop.table(table(test$readmitted, test$pred_readmit),1)
confusionMatrix(test$pred_readmit, test$readmitted)

#       <30        >30         NO
#<30 0.02874903 0.37762238 0.59362859
#>30 0.02085230 0.35953662 0.61961109
#NO  0.01393878 0.27475473 0.71130649

# correctly predicts NO readmission 71%
# correctly predicts >30 readmission 35%
# correctly predicts <30 readmission 2%


##############################################

# NAIVE BAYES 
# CoreModel implementation

# build decision tree with naive Bayes in the leaves
nbModel <- CoreModel("readmitted", train, model="tree", modelType=4)
print(nbModel)

# prediction on test set
pred <- predict(nbModel, test, type="class")
mEval <- modelEval(nbModel, test$readmitted, pred)
print(mEval) # evaluation of the model

test$pred_readmit <- pred
prop.table(table(test$readmitted, test$pred_readmit),1)
confusionMatrix(test$pred_readmit, test$readmitted)
#accuracy 0.5571548
#AUC 0.5475671
# correctly predicts NO readmission 79%
# correctly predicts >30 readmission 36%
# correctly predicts <30 readmission 3%

#                       Class: <30 Class: >30 Class: NO
#Sensitivity            0.037296     0.3614    0.7916
#Specificity            0.984547     0.7719    0.3910

##############################################

# NAIVE BAYES 
# e1071 implementation

nbayesmodel <- naiveBayes(readmitted ~ age+discharged_to+time_in_hospital+
                            num_lab_procedures+num_procedures+num_medications+number_outpatient+
                            number_emergency+number_inpatient+number_diagnoses+
                            insulin+change+diabetesMed+diag1+diag2+diag3+A1Cresult, 
                          data = train)

pred <- predict(nbayesmodel, test, type = "class")
test$pred_readmit <- pred
prop.table(table(test$readmitted, test$pred_readmit),1)
confusionMatrix(test$pred_readmit, test$readmitted)
# Overall Accuracy : 0.5597  
# correctly predicts NO readmission 91%
# correctly predicts >30 readmission 16%
# correctly predicts <30 readmission 9%

#                       Class: <30 Class: >30 Class: NO
#Sensitivity             0.08987    0.16475    0.9129
#Specificity             0.97521    0.90952    0.2202

##############################################

write.csv(data2, file = "/Users/1Air/documents/614/dataset_diabetes/processed_diabetes.csv", sep=",", na="?", row.names = F)

# fastest
# naive bayes
# logistic regression
# nnet

# kNN
# randomForest
# rpart decision tree

# slow
# kNN
# svm
