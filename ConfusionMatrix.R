#I will start by reading in a data set about breast cancer
breastcancer = read.csv(file.choose(),header=TRUE)

#here we see a list of the names in the data set and the subsequent data tied
#to the names
names(breastcancer)
head(breastcancer)

#Now I am writing a function that will install a group of packages.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

packages <- c("Hmisc","corrplot","PerformanceAnalytics","dominanceanalysis","pscl")

ipak(packages)

library(Hmisc)
library(corrplot)
library(PerformanceAnalytics)
library(pscl)
library(dominanceanalysis)



#Writing the logistic regression
#NOTE: Error "Error in eval(family$initialize) : y values must be 0 <= y <= 1" 
#popped up, so I had to change the BenignMalignant values from 2,4 to 0,1
breastcancer$BenignMalignant <- as.factor(breastcancer$BenignMalignant)

modBM <- glm(BenignMalignant ~ ClumpThickness+CellSizeUniformity+CellShapeUniformity
              +MarginalAdhesion+SingleEpithelialCell+BareNuclei+BlandChromatin+
                NormalNucleoli+Mitoses, data = 
                              breastcancer, family = binomial) 

summary(modBM)
anova(modBM, test="Chisq")

#calculating pseudo R squared
#The McFadden value is near 1, so a good model
pR2(modBM)

da.glm.fit()("names")
daBM<-dominanceAnalysis(modBM)

#impact of addition of a variable (col) to a circumstance (row)
getFits(daBM,"r2.m")

dominanceMatrix(daBM, type="complete",fit.functions = "r2.m" )
contributionByLevel(daBM,fit.functions="r2.m")
dominanceMatrix(daBM, type="conditional",fit.functions = "r2.m")
averageContribution(daBM,fit.functions = "r2.m")
dominanceMatrix(daBM, type="general",fit.functions = "r2.m")

#Preparing the data to plug into the Confusion Matrix equations
#Generating the prediction score (probability) value
breastcancer$prediction2=exp(predict(modBM,breastcancer))/1+exp(predict(modBM,breastcancer))
head(breastcancer)

#Creating a "for" loop to replace the predictive values with 2,4 (benign,malignant)
for(i in 1:699){
  if(breastcancer$prediction2[i]>0.5){
    breastcancer$prediction2[i]=4
  }
  else(breastcancer$prediction2[i] = 2)
}

#Here, we see the changes created by lines 74-79
head(breastcancer)

#Developing a table to compare errors in prediction vs actual 
#We see a resulting 97% success rate 
table(as.factor(breastcancer$BenignMalignant),as.factor(breastcancer$prediction2))

#Storing this table as an object
confmat=table(as.factor(breastcancer$BenignMalignant),as.factor(breastcancer$prediction2))
(confmat)


#Assigning Confusion Matrix variables.

#We can assume that a value of 4 (malignant/cancerous) is a positive test 
#result/prediction, and that a value of 2 (benign/non-cancerous) is a negative
#test result/prediction

#Assigning what variables I can from the confmat table
TP=confmat[2,2]
TN=confmat[1,1]
FP=confmat[1,2]
FN=confmat[2,1]

#Using vec_count from the vctrs library to find the (P) Positive and (N) Negative 
#values in BenignMalignant

library(vctrs)
BMposneg = vec_count(breastcancer$BenignMalignant)
P=BMposneg[2,2]
N=BMposneg[1,2]

#Using vec_count to find the (PP) Predicted Positive and (PN) Predicted Negative
#values in BenignMalignant
BMpredposneg = vec_count(breastcancer$prediction2)
PP=BMpredposneg[2,2]
PN=BMpredposneg[1,2]


#Plugging in Confusion Matrix variables into equations

#Total Population
TotP=P+N
#Prevalence
Prev=(P/(P+N))
#Accuracy
ACC=((TP+TN)/(P+N))
#Positive Predictive Value
PPV=TP/PP
#False Discovery Rate
FDR=FP/PP
#False Omission Rate
FOR=FN/PN
#Negative Predictive Vale
NPV=TN/PN
#True Positive Rate
TPR=TP/P
#False Negative Rate
FNR=FN/P
#False Positive Rate
FPR=FP/N
#True Negative Rate
TNR=TN/N
#Bookmaker Informedness
BM=TPR+TNR-1
#Prevalence Threshold
PT=((sqrt(TPR*FPR))-FPR)/(TPR-FPR)
#Positive Likelihood Ratio
LRp=TPR/FPR
#Negative Likelhood Ratio
LRn=FNR/TNR
#Markedness
MK=PPV+NPV-1
#Diagnostic Odds Ratio
DOR=LRp/LRn
#Balanced Accuracy
BA=(TPR+TNR)/2
#F1 Score
F1=(2*TP)/(2*TP+FP+FN)
#Fowlkes-Mallows Index
FM=sqrt(PPV*TPR)
#Matthews Correlation Coefficient
MCC=(sqrt(TPR*TNR*PPV*NPV))-(sqrt(FNR*FPR*FOR*FDR))
#Critical Success Index
CSI=(TP)/(TP+FN+FP)
