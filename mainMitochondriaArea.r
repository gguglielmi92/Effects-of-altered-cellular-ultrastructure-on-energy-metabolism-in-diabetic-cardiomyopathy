################################################################################
#Project: ageingHeart
#File: mainMitochondriaArea.r
#Author: Giovanni Guglielmi
#Description:
# -main script
################################################################################

################################################################################
#########################WORK SPACE#############################################
################################################################################

setwd("C:/Users/gguglielmi/OneDrive - The University of Melbourne/uniPhD/projects/MitochindriaAreaShourya")

################################################################################
######################### PACKAGES #############################################
################################################################################

################################ install #######################################

install.packages("nlme")
install.packages("car")

################################### call #######################################

require("nlme")
require("car")


############################### function for FPR ###############################
funcFPR <- function(t_value, estimation, std_err, dFreedom){
  center1 <- estimation/std_err
  y0 <- dt(x = t_value, df = dFreedom, ncp = 0)
  y1 <- dt(x = t_value, df = dFreedom, ncp = center1)
  lrt <- y1/(2*y0)
  fpr <- 1/(1+lrt)
  return(fpr)
}

################################################################################
######################## READ THE DATA #########################################
################################################################################

control <- read.csv(file = "controlDataset.csv", sep = ",", header = T)

diabete <- read.csv(file = "diabeteDataset.csv", sep = ",", header = T)

#set the column labels
myLabels <- c("mitoAreaFrac", "MyoAreaFrac", "glycAreaFrac", "tTubAreaFrac", 
              "ratioMitoArea","cellCrossSecSize", 
              "madArea_mitoAndMyo")


colnames(control) <- myLabels
colnames(diabete) <- myLabels

#control row names
obsControl1 <- 1:5
obsControl2 <- 6:13
obsControl3 <- 14:19

rownames(control)[obsControl1] <- paste0("control1_", 1:length(obsControl1)) 
rownames(control)[obsControl2] <- paste0("control2_", 1:length(obsControl2))
rownames(control)[obsControl3] <- paste0("control3_", 1:length(obsControl3))

#diabete row names
obsDiabete1 <- 1:8
obsDiabete2 <- 9:16
obsDiabete3 <- 17:21

rownames(diabete)[obsDiabete1] <- paste0("diabete1_", 1:length(obsDiabete1)) 
rownames(diabete)[obsDiabete2] <- paste0("diabete2_", 1:length(obsDiabete2))
rownames(diabete)[obsDiabete3] <- paste0("diabete3_", 1:length(obsDiabete3))

myData_raw <- rbind(control, diabete)
myData_raw$condition <- factor(x = c(rep("control", times = nrow(control)), 
                                     rep("diabete", times = nrow(diabete))),  
                               levels = c("control", "diabete"))
myData_raw$technical <- factor(x = c(rep("control1", times = length(obsControl1)), 
                                     rep("control2", times = length(obsControl2)), 
                                     rep("control3", times = length(obsControl3)), 
                                     rep("diabete1", times = length(obsDiabete1)), 
                                     rep("diabete2", times = length(obsDiabete2)), 
                                     rep("diabete3", times = length(obsDiabete3))), 
                               levels = c("control1", "control2", "control3", 
                                          "diabete1", "diabete2", "diabete3"))
#summary
summary(myData_raw)
class(myData_raw$condition)


rm(obsControl1, obsControl2, obsControl3, 
   obsDiabete1, obsDiabete2, obsDiabete3, 
   myLabels, 
   control, diabete)

######################## set the data ######################################

myData <- data.frame(mitoAreaFrac = double(), 
                     MyoAreaFrac = double(), 
                     glycAreaFrac = double(),
                     tTubAreaFrac = double(), 
                     ratioMitoArea = double(), 
                     cellCrossSecSize = double(), 
                     madArea_mitoAndMyo = double())

#get the control biological samples
myData[1, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "control1"), 1:7], 
                               MARGIN = 2, FUN = mean))
myData[2, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "control2"), 1:7],
                               MARGIN = 2, FUN = mean))
myData[3, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "control3"), 1:7],
                               MARGIN = 2, FUN = mean))

#get the diabete biological sample
myData[4, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "diabete1"), 1:7], 
                               MARGIN = 2, FUN = mean))
myData[5, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "diabete2"), 1:7],
                               MARGIN = 2, FUN = mean))
myData[6, ] <- as.vector(apply(X = myData_raw[which(myData_raw$technical == "diabete3"), 1:7],
                               MARGIN = 2, FUN = mean))

#now attach the categories
myData$condition <- factor(x = c(rep("control", times = 3), 
                                 rep("diabete", times = 3)), 
                           levels = c("control", "diabete"))

#rownames
rownames(myData) <- c(paste0("control", 1:3), paste0("diabete", 1:3))

########################### create dataset of results ##########################
myResults <- data.frame(QUESTION	= NA, 
                        STATISTICS = NA,
                        ESTIMATION = NA, 
                        STD_ERR = NA ,
                        CI_LOW = NA, 
                        CI_UP = NA,
                        t_VALUE = NA, 
                        p_VALUE	= NA, 
                        FPR = NA)

################################################################################
#q1: Is there any correlation between area fractions of mitochondria vs        #
#myofibrils. In each cross section                                             #
################################################################################

#correlation
myCorr <- cor.test(x = myData$mitoAreaFrac, y = myData$MyoAreaFrac, method = "pearson")

myResults$QUESTION[1] <- "Correlation area fractions of mitochondria vs myofibrils"
myResults$STATISTICS[1] <- "t_test df = 4"
myResults$ESTIMATION[1] <- myCorr$estimate
myResults$STD_ERR[1] <- sqrt((1 - myCorr$estimate^2)/(nrow(myData) - 2))
myResults$CI_LOW[1] <- myCorr$conf.int[1]
myResults$CI_UP[1] <- myCorr$conf.int[2]
myResults$t_VALUE[1] <- myCorr$statistic
myResults$p_VALUE[1] <- myCorr$p.value
myResults$FPR[1] <- funcFPR(t_value = myResults$t_VALUE[1], 
                            estimation = myResults$ESTIMATION[1], 
                            std_err = myResults$STD_ERR[1],
                            dFreedom = (nrow(myData)-2))

myCorr
rm(myCorr)

################################################################################
#q2: Is there any correlation between area fractions of mitochondria vs        #
#glycogen in each cross section ?                                              #
################################################################################

#correlation
myCorr <- cor.test(x = myData$mitoAreaFrac, y = myData$glycAreaFrac, method = "pearson")

myResults[2, ] <- NA 
myResults$QUESTION[2] <- "correlation between area fractions of mitochondria vs glycogen"
myResults$STATISTICS[2] <- "t_test df = 4"
myResults$ESTIMATION[2] <- myCorr$estimate
myResults$STD_ERR[2] <- sqrt((1 - myCorr$estimate^2)/(nrow(myData) - 2))
myResults$CI_LOW[2] <- myCorr$conf.int[1]
myResults$CI_UP[2] <- myCorr$conf.int[2]
myResults$t_VALUE[2] <- myCorr$statistic
myResults$p_VALUE[2] <- myCorr$p.value
myResults$FPR[2] <- funcFPR(t_value = myResults$t_VALUE[2], 
                            estimation = myResults$ESTIMATION[2], 
                            std_err = myResults$STD_ERR[2],
                            dFreedom = (nrow(myData)-2))

myCorr
rm(myCorr)
################################################################################
#q3: Is there any correlation between area fractions of mitochondria vs        #
#T-tubules in each cross section ?                                             #
################################################################################

#correlation
myCorr <- cor.test(x = myData$mitoAreaFrac, y = myData$tTubAreaFrac, method = "pearson")

myResults[3, ] <- NA 
myResults$QUESTION[3] <- "correlation between area fractions of mitochondria vs T-tubules"
myResults$STATISTICS[3] <- "t_test df = 4"
myResults$ESTIMATION[3] <- myCorr$estimate
myResults$STD_ERR[3] <- sqrt((1 - myCorr$estimate^2)/(nrow(myData) - 2))
myResults$CI_LOW[3] <- myCorr$conf.int[1]
myResults$CI_UP[3] <- myCorr$conf.int[2]
myResults$t_VALUE[3] <- myCorr$statistic
myResults$p_VALUE[3] <- myCorr$p.value
myResults$FPR[3] <- funcFPR(t_value = myResults$t_VALUE[3], 
                            estimation = myResults$ESTIMATION[3], 
                            std_err = myResults$STD_ERR[3],
                            dFreedom = (nrow(myData)-2))

myCorr
rm(myCorr)

################### LMM APPROACH ###############################################

################################################################################
#q4: Is the overall mitochondria / (myofibril + mitochondria) ratio different  #
#    between control and diabetic groups ? I mean is there any significant     #
#    difference between column E of both file ?                                #
################################################################################
myLme <- lme(fixed = ratioMitoArea ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[4, ] <- NA 
myResults$QUESTION[4] <- "mitochondria / (myofibril + mitochondria) different between control and diabetic"
myResults$STATISTICS[4] <- "t_test df = 4"
myResults$ESTIMATION[4] <- mySumm$tTable[2, 1]
myResults$STD_ERR[4] <- mySumm$tTable[2, 2]
myResults$CI_LOW[4] <-  Confint(myLme)[2, 2]
myResults$CI_UP[4] <- Confint(myLme)[2, 3]
myResults$t_VALUE[4] <- mySumm$tTable[2, 4]
myResults$p_VALUE[4] <- mySumm$tTable[2, 5]
myResults$FPR[4] <- funcFPR(t_value = myResults$t_VALUE[4], 
                            estimation = myResults$ESTIMATION[4], 
                            std_err = myResults$STD_ERR[4],
                            dFreedom = mySumm$tTable[2, 3])


#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)

################################################################################
#q5: Is the MAD of mitochondrial area distribution different between control   #
#    and diabetic groups ? i.e.,  is there any significant difference between  #
#    column G of both file ??                                                  #
################################################################################

myLme <- lme(fixed = madArea_mitoAndMyo ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[5, ] <- NA 
myResults$QUESTION[5] <- "MAD of mitochondrial area distribution different between control and diabetic groups"
myResults$STATISTICS[5] <- "t_test df = 4"
myResults$ESTIMATION[5] <- mySumm$tTable[2, 1]
myResults$STD_ERR[5] <- mySumm$tTable[2, 2]
myResults$CI_LOW[5] <-  Confint(myLme)[2, 2]
myResults$CI_UP[5] <- Confint(myLme)[2, 3]
myResults$t_VALUE[5] <- mySumm$tTable[2, 4]
myResults$p_VALUE[5] <- mySumm$tTable[2, 5]
myResults$FPR[5] <- funcFPR(t_value = myResults$t_VALUE[5], 
                            estimation = myResults$ESTIMATION[5], 
                            std_err = myResults$STD_ERR[5],
                            dFreedom = mySumm$tTable[2, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)
################################################################################
#q6 : Is there any significant difference between column A                     #
#    (Mitochondria area fraction) of both file ?                               #
################################################################################

myLme <- lme(fixed = mitoAreaFrac ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[6, ] <- NA 
myResults$QUESTION[6] <- "difference mitochondria area fraction between control and diabetic groups"
myResults$STATISTICS[6] <- "t_test df = 4"
myResults$ESTIMATION[6] <- mySumm$tTable[2, 1]
myResults$STD_ERR[6] <- mySumm$tTable[2, 2]
myResults$CI_LOW[6] <-  Confint(myLme)[2, 2]
myResults$CI_UP[6] <- Confint(myLme)[2, 3]
myResults$t_VALUE[6] <- mySumm$tTable[2, 4]
myResults$p_VALUE[6] <- mySumm$tTable[2, 5]
myResults$FPR[6] <- funcFPR(t_value = myResults$t_VALUE[6], 
                            estimation = myResults$ESTIMATION[6], 
                            std_err = myResults$STD_ERR[6],
                            dFreedom = mySumm$tTable[2, 3])


#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)

################################################################################
#q7 : Is there any significant difference between column B                     #
#     (Myofibrils area fraction) of both file                                  #
################################################################################

myLme <- lme(fixed = MyoAreaFrac ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[7, ] <- NA 
myResults$QUESTION[7] <- "difference Myofibrils area fraction between control and diabetic groups"
myResults$STATISTICS[7] <- "t_test df = 4"
myResults$ESTIMATION[7] <- mySumm$tTable[2, 1]
myResults$STD_ERR[7] <- mySumm$tTable[2, 2]
myResults$CI_LOW[7] <-  Confint(myLme)[2, 2]
myResults$CI_UP[7] <- Confint(myLme)[2, 3]
myResults$t_VALUE[7] <- mySumm$tTable[2, 4]
myResults$p_VALUE[7] <- mySumm$tTable[2, 5]
myResults$FPR[7] <- funcFPR(t_value = myResults$t_VALUE[7], 
                            estimation = myResults$ESTIMATION[7], 
                            std_err = myResults$STD_ERR[7],
                            dFreedom = mySumm$tTable[2, 3])


#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)


boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)
################################################################################
#q8 : Is there any significant difference between column C (Glycogen area      #
#     fraction) of both file                                                   #
################################################################################

myLme <- lme(fixed = glycAreaFrac ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[8, ] <- NA 
myResults$QUESTION[8] <- "difference Glycogen area fraction between control and diabetic groups"
myResults$STATISTICS[8] <- "t_test df = 4"
myResults$ESTIMATION[8] <- mySumm$tTable[2, 1]
myResults$STD_ERR[8] <- mySumm$tTable[2, 2]
myResults$CI_LOW[8] <-  Confint(myLme)[2, 2]
myResults$CI_UP[8] <- Confint(myLme)[2, 3]
myResults$t_VALUE[8] <- mySumm$tTable[2, 4]
myResults$p_VALUE[8] <- mySumm$tTable[2, 5]
myResults$FPR[8] <- funcFPR(t_value = myResults$t_VALUE[8], 
                            estimation = myResults$ESTIMATION[8], 
                            std_err = myResults$STD_ERR[8],
                            dFreedom = mySumm$tTable[2, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)

################################################################################
#q9 : Is there any significant difference between column D (T-tubules area     #
#     fraction) of both file                                                   #
################################################################################

myLme <- lme(tTubAreaFrac ~ condition, 
             random = ~1 | technical, 
             data = myData_raw, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[9, ] <- NA 
myResults$QUESTION[9] <- "difference T-tubules area fraction between control and diabetic groups"
myResults$STATISTICS[9] <- "t_test df = 4"
myResults$ESTIMATION[9] <- mySumm$tTable[2, 1]
myResults$STD_ERR[9] <- mySumm$tTable[2, 2]
myResults$CI_LOW[9] <-  Confint(myLme)[2, 2]
myResults$CI_UP[9] <- Confint(myLme)[2, 3]
myResults$t_VALUE[9] <- mySumm$tTable[2, 4]
myResults$p_VALUE[9] <- mySumm$tTable[2, 5]
myResults$FPR[9] <- funcFPR(t_value = myResults$t_VALUE[9], 
                            estimation = myResults$ESTIMATION[9], 
                            std_err = myResults$STD_ERR[9],
                            dFreedom = mySumm$tTable[2, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~myData_raw$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(myData_raw$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(myData_raw$technical)){
  points(x = as.numeric(myData_raw$condition[myData_raw$technical == levels(myData_raw$technical)[i]])[1], 
         y = fitSub[myData_raw$technical == levels(myData_raw$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

#remove all the previos stuff
rm(e_hat, fitOverall, fitSub, i, myLme, mySumm)
rm(myData, myData_raw)

#                                                                              #
#                                                                              #
#                                                                              #
############################### NEW SET OF QUESTIONS ###########################
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################
######################## READ THE DATA #########################################
################################################################################

############################ conditions datasets ###############################

cc <- read.csv("simCC.csv")
cd <- read.csv("simCD.csv")
dd <- read.csv("simDD.csv")

#compact the three conditions datasets

condData <- rbind(cc, cd, dd) 
condData$technical <- factor(x = condData$technical, 
                             levels = c("cc1", "cc2", "cc3", 
                                        "cd1", "cd2", "cd3", 
                                        "dd1", "dd2", "dd3"))
condData$condition <- factor(x = condData$condition, 
                             levels = c("cc", "cd", "dd"))

summary(condData)

#remove the datasets of the different conditions
rm(cc, cd, dd)

########################### parameters dataset #################################

paramControl <- read.csv("simControlStructureParameters.csv") 
paramDiabete <- read.csv("simDiabeteStructureParameters.csv")

########### compact the parameters
paramData <- rbind(paramControl, paramDiabete)
paramData$technical <- factor(paramData$technical, 
                              levels = c("control1", "control2", "control3", 
                                         "diabete1", "diabete2", "diabete3"))
paramData$condition <- factor(x = paramData$condition, 
                              levels = c("control", "diabete"))

summary(paramData)

#remove the dataset of the structures
rm(paramControl, paramDiabete)

################################################################################
#q1 : Is M/M MEAN significantly different between control and diabetes?        #
################################################################################
myLme <- lme(MM_MEAN ~ condition, 
             random = ~1|technical, 
             data = paramData, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[10, ] <- NA 
myResults$QUESTION[10] <- "difference M/M MEAN between control and diabetic groups"
myResults$STATISTICS[10] <- "t_test df = 4"
myResults$ESTIMATION[10] <- mySumm$tTable[2, 1]
myResults$STD_ERR[10] <- mySumm$tTable[2, 2]
myResults$CI_LOW[10] <-  Confint(myLme)[2, 2]
myResults$CI_UP[10] <- Confint(myLme)[2, 3]
myResults$t_VALUE[10] <- mySumm$tTable[2, 4]
myResults$p_VALUE[10] <- mySumm$tTable[2, 5]
myResults$FPR[10] <- funcFPR(t_value = myResults$t_VALUE[10], 
                             estimation = myResults$ESTIMATION[10], 
                             std_err = myResults$STD_ERR[10],
                             dFreedom = mySumm$tTable[2, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~paramData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(paramData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(paramData$technical)){
  points(x = as.numeric(paramData$condition[paramData$technical == levels(paramData$technical)[i]])[1], 
         y = fitSub[paramData$technical == levels(paramData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

#remove all the previos stuff
rm(e_hat, fitOverall, fitSub, i, myLme, mySumm)

################################################################################
#q2: Is M/M MAD significantly different between control and diabetes?          #
################################################################################

myLme <- lme(MM_MAD ~ condition, 
             random = ~1|technical, 
             data = paramData, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[11, ] <- NA 
myResults$QUESTION[11] <- "difference M/M MAD between control and diabetic groups"
myResults$STATISTICS[11] <- "t_test df = 4"
myResults$ESTIMATION[11] <- mySumm$tTable[2, 1]
myResults$STD_ERR[11] <- mySumm$tTable[2, 2]
myResults$CI_LOW[11] <-  Confint(myLme)[2, 2]
myResults$CI_UP[11] <- Confint(myLme)[2, 3]
myResults$t_VALUE[11] <- mySumm$tTable[2, 4]
myResults$p_VALUE[11] <- mySumm$tTable[2, 5]
myResults$FPR[11] <- funcFPR(t_value = myResults$t_VALUE[11], 
                             estimation = myResults$ESTIMATION[11], 
                             std_err = myResults$STD_ERR[11],
                             dFreedom = mySumm$tTable[2, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~paramData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(paramData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(paramData$technical)){
  points(x = as.numeric(paramData$condition[paramData$technical == levels(paramData$technical)[i]])[1], 
         y = fitSub[paramData$technical == levels(paramData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

#remove all the previos stuff
rm(e_hat, fitOverall, fitSub, i, myLme,mySumm)


################################################################################
#q3: Is MEAN of ADP/ATP significantly different amongst CC, CD and DD?         #
################################################################################

myLme <- lme(mean_AdpAtp ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[12, ] <- NA 
myResults$QUESTION[12] <- "MEAN ADP/ATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[12] <- "CD:  t_test df = 6"
myResults$ESTIMATION[12] <- mySumm$tTable[2, 1]
myResults$STD_ERR[12] <- mySumm$tTable[2, 2]
myResults$CI_LOW[12] <-  Confint(myLme)[2, 2]
myResults$CI_UP[12] <- Confint(myLme)[2, 3]
myResults$t_VALUE[12] <- mySumm$tTable[2, 4]
myResults$p_VALUE[12] <- mySumm$tTable[2, 5]
myResults$FPR[12] <- funcFPR(t_value = myResults$t_VALUE[12], 
                             estimation = myResults$ESTIMATION[12], 
                             std_err = myResults$STD_ERR[12],
                             dFreedom = mySumm$tTable[2, 3])

myResults[13, ] <- NA 
myResults$QUESTION[13] <- "MEAN ADP/ATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[13] <- "DD:  t_test df = 6"
myResults$ESTIMATION[13] <- mySumm$tTable[3, 1]
myResults$STD_ERR[13] <- mySumm$tTable[3, 2]
myResults$CI_LOW[13] <-  Confint(myLme)[3, 2]
myResults$CI_UP[13] <- Confint(myLme)[3, 3]
myResults$t_VALUE[13] <- mySumm$tTable[3, 4]
myResults$p_VALUE[13] <- mySumm$tTable[3, 5]
myResults$FPR[13] <- funcFPR(t_value = myResults$t_VALUE[13], 
                             estimation = myResults$ESTIMATION[13], 
                             std_err = myResults$STD_ERR[13],
                             dFreedom = mySumm$tTable[3, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm, e_hat, fitOverall, fitSub, i)
################################################################################
#q4: Is MEAN of Pi significantly different amongst CC, CD and DD?              #
################################################################################

myLme <- lme(mean_Pi ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[14, ] <- NA 
myResults$QUESTION[14] <- "MEAN Pi significantly different amongst CC, CD and DD"
myResults$STATISTICS[14] <- "CD:  t_test df = 6"
myResults$ESTIMATION[14] <- mySumm$tTable[2, 1]
myResults$STD_ERR[14] <- mySumm$tTable[2, 2]
myResults$CI_LOW[14] <-  Confint(myLme)[2, 2]
myResults$CI_UP[14] <- Confint(myLme)[2, 3]
myResults$t_VALUE[14] <- mySumm$tTable[2, 4]
myResults$p_VALUE[14] <- mySumm$tTable[2, 5]
myResults$FPR[14] <- funcFPR(t_value = myResults$t_VALUE[14], 
                             estimation = myResults$ESTIMATION[14], 
                             std_err = myResults$STD_ERR[14],
                             dFreedom = mySumm$tTable[2, 3])

myResults[15, ] <- NA 
myResults$QUESTION[15] <- "MEAN Pi significantly different amongst CC, CD and DD"
myResults$STATISTICS[15] <- "DD:  t_test df = 6"
myResults$ESTIMATION[15] <- mySumm$tTable[3, 1]
myResults$STD_ERR[15] <- mySumm$tTable[3, 2]
myResults$CI_LOW[15] <-  Confint(myLme)[3, 2]
myResults$CI_UP[15] <- Confint(myLme)[3, 3]
myResults$t_VALUE[15] <- mySumm$tTable[3, 4]
myResults$p_VALUE[15] <- mySumm$tTable[3, 5]
myResults$FPR[15] <- funcFPR(t_value = myResults$t_VALUE[15], 
                             estimation = myResults$ESTIMATION[15], 
                             std_err = myResults$STD_ERR[15],
                             dFreedom = mySumm$tTable[3, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm, e_hat, fitOverall, fitSub, i)

################################################################################
#q5: Is MEAN of XATP significantly different amongst CC, CD and DD?            #
################################################################################

myLme <- lme(mean_XATP ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[16, ] <- NA 
myResults$QUESTION[16] <- "MEAN of XATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[16] <- "CD:  t_test df = 6"
myResults$ESTIMATION[16] <- mySumm$tTable[2, 1]
myResults$STD_ERR[16] <- mySumm$tTable[2, 2]
myResults$CI_LOW[16] <-  Confint(myLme)[2, 2]
myResults$CI_UP[16] <- Confint(myLme)[2, 3]
myResults$t_VALUE[16] <- mySumm$tTable[2, 4]
myResults$p_VALUE[16] <- mySumm$tTable[2, 5]
myResults$FPR[16] <- funcFPR(t_value = myResults$t_VALUE[16], 
                             estimation = myResults$ESTIMATION[16], 
                             std_err = myResults$STD_ERR[16],
                             dFreedom = mySumm$tTable[2, 3])

myResults[17, ] <- NA 
myResults$QUESTION[17] <- "MEAN of XATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[17] <- "DD:  t_test df = 6"
myResults$ESTIMATION[17] <- mySumm$tTable[3, 1]
myResults$STD_ERR[17] <- mySumm$tTable[3, 2]
myResults$CI_LOW[17] <-  Confint(myLme)[3, 2]
myResults$CI_UP[17] <- Confint(myLme)[3, 3]
myResults$t_VALUE[17] <- mySumm$tTable[3, 4]
myResults$p_VALUE[17] <- mySumm$tTable[3, 5]
myResults$FPR[17] <- funcFPR(t_value = myResults$t_VALUE[17], 
                             estimation = myResults$ESTIMATION[17], 
                             std_err = myResults$STD_ERR[17],
                             dFreedom = mySumm$tTable[3, 3])


#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-1, max(fitSub)))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(myLme, mySumm)
################################################################################
#q6: Is MAD of ADP/ATP significantly different amongst CC, CD and DD?          #
################################################################################

myLme <- lme(mad_AdpAtp ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")

mySumm <- summary(myLme)
mySumm

myResults[18, ] <- NA 
myResults$QUESTION[18] <- "MAD of ADP/ATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[18] <- "CD:  t_test df = 6"
myResults$ESTIMATION[18] <- mySumm$tTable[2, 1]
myResults$STD_ERR[18] <- mySumm$tTable[2, 2]
myResults$CI_LOW[18] <-  Confint(myLme)[2, 2]
myResults$CI_UP[18] <- Confint(myLme)[2, 3]
myResults$t_VALUE[18] <- mySumm$tTable[2, 4]
myResults$p_VALUE[18] <- mySumm$tTable[2, 5]
myResults$FPR[18] <- funcFPR(t_value = myResults$t_VALUE[18], 
                             estimation = myResults$ESTIMATION[18], 
                             std_err = myResults$STD_ERR[18],
                             dFreedom = mySumm$tTable[2, 3])

myResults[19, ] <- NA 
myResults$QUESTION[19] <- "MAD of ADP/ATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[19] <- "DD:  t_test df = 6"
myResults$ESTIMATION[19] <- mySumm$tTable[3, 1]
myResults$STD_ERR[19] <- mySumm$tTable[3, 2]
myResults$CI_LOW[19] <-  Confint(myLme)[3, 2]
myResults$CI_UP[19] <- Confint(myLme)[3, 3]
myResults$t_VALUE[19] <- mySumm$tTable[3, 4]
myResults$p_VALUE[19] <- mySumm$tTable[3, 5]
myResults$FPR[19] <- funcFPR(t_value = myResults$t_VALUE[19], 
                             estimation = myResults$ESTIMATION[19], 
                             std_err = myResults$STD_ERR[19],
                             dFreedom = mySumm$tTable[3, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-0.05, max(fitSub)+0.05))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(mySumm, myLme)
################################################################################
#q7: Is MAD of Pi significantly different amongst CC, CD and DD?               #
################################################################################

myLme <- lme(mad_Pi ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[20, ] <- NA 
myResults$QUESTION[20] <- "MAD of Pi significantly different amongst CC, CD and DD"
myResults$STATISTICS[20] <- "CD:  t_test df = 6"
myResults$ESTIMATION[20] <- mySumm$tTable[2, 1]
myResults$STD_ERR[20] <- mySumm$tTable[2, 2]
myResults$CI_LOW[20] <-  Confint(myLme)[2, 2]
myResults$CI_UP[20] <- Confint(myLme)[2, 3]
myResults$t_VALUE[20] <- mySumm$tTable[2, 4]
myResults$p_VALUE[20] <- mySumm$tTable[2, 5]
myResults$FPR[20] <- funcFPR(t_value = myResults$t_VALUE[20], 
                             estimation = myResults$ESTIMATION[20], 
                             std_err = myResults$STD_ERR[20],
                             dFreedom = mySumm$tTable[2, 3])

myResults[21, ] <- NA 
myResults$QUESTION[21] <- "MAD of Pi significantly different amongst CC, CD and DD"
myResults$STATISTICS[21] <- "DD:  t_test df = 6"
myResults$ESTIMATION[21] <- mySumm$tTable[3, 1]
myResults$STD_ERR[21] <- mySumm$tTable[3, 2]
myResults$CI_LOW[21] <-  Confint(myLme)[3, 2]
myResults$CI_UP[21] <- Confint(myLme)[3, 3]
myResults$t_VALUE[21] <- mySumm$tTable[3, 4]
myResults$p_VALUE[21] <- mySumm$tTable[3, 5]
myResults$FPR[21] <- funcFPR(t_value = myResults$t_VALUE[21], 
                             estimation = myResults$ESTIMATION[21], 
                             std_err = myResults$STD_ERR[21],
                             dFreedom = mySumm$tTable[3, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-0.05, max(fitSub)+0.05))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}

rm(mySumm, myLme)
################################################################################
#q8: Is MAD of XATP significantly different amongst CC, CD and DD?             #
################################################################################

myLme <- lme(mad_XATP ~ condition, 
             random = ~1|technical, 
             data = condData, method = "REML")
mySumm <- summary(myLme)
mySumm

myResults[22, ] <- NA 
myResults$QUESTION[22] <- "MAD of XATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[22] <- "CD:  t_test df = 6"
myResults$ESTIMATION[22] <- mySumm$tTable[2, 1]
myResults$STD_ERR[22] <- mySumm$tTable[2, 2]
myResults$CI_LOW[22] <-  Confint(myLme)[2, 2]
myResults$CI_UP[22] <- Confint(myLme)[2, 3]
myResults$t_VALUE[22] <- mySumm$tTable[2, 4]
myResults$p_VALUE[22] <- mySumm$tTable[2, 5]
myResults$FPR[22] <- funcFPR(t_value = myResults$t_VALUE[22], 
                             estimation = myResults$ESTIMATION[22], 
                             std_err = myResults$STD_ERR[22],
                             dFreedom = mySumm$tTable[2, 3])

myResults[23, ] <- NA 
myResults$QUESTION[23] <- "MAD of XATP significantly different amongst CC, CD and DD"
myResults$STATISTICS[23] <- "DD:  t_test df = 6"
myResults$ESTIMATION[23] <- mySumm$tTable[3, 1]
myResults$STD_ERR[23] <- mySumm$tTable[3, 2]
myResults$CI_LOW[23] <-  Confint(myLme)[3, 2]
myResults$CI_UP[23] <- Confint(myLme)[3, 3]
myResults$t_VALUE[23] <- mySumm$tTable[3, 4]
myResults$p_VALUE[23] <- mySumm$tTable[3, 5]
myResults$FPR[23] <- funcFPR(t_value = myResults$t_VALUE[23], 
                             estimation = myResults$ESTIMATION[23], 
                             std_err = myResults$STD_ERR[23],
                             dFreedom = mySumm$tTable[3, 3])

#diagnostic
e_hat <- resid(myLme, type = "normalized")
qqnorm(e_hat)
qqline(e_hat)
hist(e_hat)
shapiro.test(x = e_hat)

boxplot(e_hat~condData$condition)
plot(predict(myLme), e_hat)

#prediction
fitOverall <- fitted(myLme, level = 0)
fitSub <- fitted(myLme, level = 1)
plot(x = as.numeric(condData$condition), y = fitOverall, 
     col = "red", type = "p", cex = 3, pch = 16, 
     ylim = c(min(fitSub)-0.05, max(fitSub)+0.05))
for (i in 1:nlevels(condData$technical)){
  points(x = as.numeric(condData$condition[condData$technical == levels(condData$technical)[i]])[1], 
         y = fitSub[condData$technical == levels(condData$technical)[i]][1], 
         cex = 2, col = "blue", pch = 16)
}
rm(myLme, mySumm)
#                                                                              #
#                                                                              #
#                                                                              #
################################################################################
############################## correlation analysis ############################
################################################################################
#                                                                              #
#                                                                              #
#                                                                              #
######################## prepare average data for condition ####################

variableToExclude <- c("technical", "condition")

condData_avg = data.frame(mean_AdpAtp = numeric(), 
                          mean_Pi = numeric(), 
                          mean_XATP = numeric(), 
                          mad_AdpAtp = numeric(), 
                          mad_Pi = numeric(), 
                          mad_XATP = numeric())
              
#cc
condData_avg[1, ] <- as.vector(apply(X = condData[condData$technical=="cc1", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[2, ] <- as.vector(apply(X = condData[condData$technical=="cc2", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[3, ] <- as.vector(apply(X = condData[condData$technical=="cc3", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))

#cd
condData_avg[4, ] <- as.vector(apply(X = condData[condData$technical=="cd1", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[5, ] <- as.vector(apply(X = condData[condData$technical=="cd2", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[6, ] <- as.vector(apply(X = condData[condData$technical=="cd3", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
#dd
condData_avg[7, ] <- as.vector(apply(X = condData[condData$technical=="dd1", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[8, ] <- as.vector(apply(X = condData[condData$technical=="dd2", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
condData_avg[9, ] <- as.vector(apply(X = condData[condData$technical=="dd3", 
                                                  colnames(condData)!= variableToExclude], 
                                     MARGIN = 2, FUN = mean))
#add the condition
condData_avg$condition <- factor(c(rep(x = "cc", times = 3), 
                                   rep(x = "cd", times = 3), 
                                   rep(x = "dd", times = 3)), 
                                 levels = c("cc", "cd", "dd"))
rm(variableToExclude)
################################ create parameter data #########################
paramData_avg <- data.frame(MM_MEAN = numeric() ,
                            MM_MAD = numeric())

#control
paramData_avg[1, ] <- as.vector(apply(X = paramData[paramData$technical == "control1", 
                                                c("MM_MEAN", "MM_MAD")], 
                                MARGIN = 2, FUN = mean))

paramData_avg[2, ] <- as.vector(apply(X = paramData[paramData$technical == "control2", 
                                                    c("MM_MEAN", "MM_MAD")], 
                                      MARGIN = 2, FUN = mean))

paramData_avg[3, ] <- as.vector(apply(X = paramData[paramData$technical == "control3", 
                                                    c("MM_MEAN", "MM_MAD")], 
                                      MARGIN = 2, FUN = mean))
#diabete
paramData_avg[4, ] <- as.vector(apply(X = paramData[paramData$technical == "diabete1", 
                                                    c("MM_MEAN", "MM_MAD")], 
                                      MARGIN = 2, FUN = mean))
paramData_avg[5, ] <- as.vector(apply(X = paramData[paramData$technical == "diabete2", 
                                                    c("MM_MEAN", "MM_MAD")], 
                                      MARGIN = 2, FUN = mean))
paramData_avg[6, ] <- as.vector(apply(X = paramData[paramData$technical == "diabete3", 
                                                    c("MM_MEAN", "MM_MAD")], 
                                      MARGIN = 2, FUN = mean))

paramData_avg$condition <- factor(c(rep("control", times = 3),
                                    rep("diabete", times = 3)), 
                                  levels = c("control", "diabete"))

################################################################################
#q9-11: Does M/M MEAN correlate with correlate  variables in CC, CD and DD      #
################################################################################

#################### CC ########################################################
myApply <- apply(X = condData[condData$condition == "cc", 2:4], MARGIN = 2, 
                 FUN = function(a){
                        cor.test(x = a, 
                        y = paramData[paramData$condition == "control", c("MM_MEAN")], 
                        method = "pearson")
                        })

myResults[24:26, ] <- NA

myApply$mean_AdpAtp

myResults$QUESTION[24] <- "Correlation M/M MEAN with mean_AdpAtp in control - cc"
myResults$STATISTICS[24] <- "t_test df = 17"
myResults$ESTIMATION[24] <- myApply$mean_AdpAtp$estimate
myResults$STD_ERR[24] <- sqrt((1 - myResults$ESTIMATION[24]^2)/myApply$mean_AdpAtp$parameter)
myResults$CI_LOW[24] <- myResults$ESTIMATION[24] + qt(p = 0.025, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[24]
myResults$CI_UP[24] <- myResults$ESTIMATION[24] + qt(p = 0.975, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[24]
myResults$t_VALUE[24] <- myApply$mean_AdpAtp$statistic
myResults$p_VALUE[24] <- myApply$mean_AdpAtp$p.value
myResults$FPR[24] <- funcFPR(t_value = myResults$t_VALUE[24], 
                             estimation = myResults$ESTIMATION[24], 
                             std_err = myResults$STD_ERR[24],
                             dFreedom = myApply$mean_AdpAtp$parameter)

myResults$QUESTION[25] <- "Correlation M/M MEAN with mean_Pi in control - cc"
myResults$STATISTICS[25] <- "t_test df = 17"
myResults$ESTIMATION[25] <- myApply$mean_Pi$estimate
myResults$STD_ERR[25] <- sqrt((1 - myResults$ESTIMATION[25]^2)/myApply$mean_Pi$parameter)
myResults$CI_LOW[25] <- myResults$ESTIMATION[25] + qt(p = 0.025, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[25]
myResults$CI_UP[25] <- myResults$ESTIMATION[25] + qt(p = 0.975, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[25]
myResults$t_VALUE[25] <- myApply$mean_Pi$statistic
myResults$p_VALUE[25] <- myApply$mean_Pi$p.value
myResults$FPR[25] <- funcFPR(t_value = myResults$t_VALUE[25], 
                             estimation = myResults$ESTIMATION[25], 
                             std_err = myResults$STD_ERR[25],
                             dFreedom = myApply$mean_Pi$parameter)

myResults$QUESTION[26] <- "Correlation M/M MEAN with mean_XATP in control - cc"
myResults$STATISTICS[26] <- "t_test df = 17"
myResults$ESTIMATION[26] <- myApply$mean_XATP$estimate
myResults$STD_ERR[26] <- sqrt((1 - myResults$ESTIMATION[26]^2)/myApply$mean_XATP$parameter)
myResults$CI_LOW[26] <- myResults$ESTIMATION[26] + qt(p = 0.025, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[26]
myResults$CI_UP[26] <- myResults$ESTIMATION[26] + qt(p = 0.975, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[26]
myResults$t_VALUE[26] <- myApply$mean_XATP$statistic
myResults$p_VALUE[26] <- myApply$mean_XATP$p.value
myResults$FPR[26] <- funcFPR(t_value = myResults$t_VALUE[26], 
                             estimation = myResults$ESTIMATION[26], 
                             std_err = myResults$STD_ERR[26],
                             dFreedom = myApply$mean_XATP$parameter)

rm(myApply)
#################### CD ########################################################
myApply <- apply(X = condData[condData$condition == "cd", 2:4], MARGIN = 2, 
                 FUN = function(a){
                        cor.test(x = a, 
                                 y = paramData[paramData$condition == "control", c("MM_MEAN")], 
                                 method = "pearson")
                        })

myResults[27:29, ] <- NA

myResults$QUESTION[27] <- "Correlation M/M MEAN with mean_AdpAtp in control - cd"
myResults$STATISTICS[27] <- "t_test df = 17"
myResults$ESTIMATION[27] <- myApply$mean_AdpAtp$estimate
myResults$STD_ERR[27] <- sqrt((1 - myResults$ESTIMATION[27]^2)/(myApply$mean_AdpAtp$parameter))
myResults$CI_LOW[27] <- myResults$ESTIMATION[27] + qt(p = 0.025, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[27]
myResults$CI_UP[27] <- myResults$ESTIMATION[27] + qt(p = 0.975, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[27]
myResults$t_VALUE[27] <- myApply$mean_AdpAtp$statistic
myResults$p_VALUE[27] <- myApply$mean_AdpAtp$p.value
myResults$FPR[27] <- funcFPR(t_value = myResults$t_VALUE[27], 
                             estimation = myResults$ESTIMATION[27], 
                             std_err = myResults$STD_ERR[27],
                             dFreedom = myApply$mean_AdpAtp$parameter)

myResults$QUESTION[28] <- "Correlation M/M MEAN with mean_Pi in control - cd"
myResults$STATISTICS[28] <- "t_test df = 17"
myResults$ESTIMATION[28] <- myApply$mean_Pi$estimate
myResults$STD_ERR[28] <- sqrt((1 - myResults$ESTIMATION[28]^2)/myApply$mean_Pi$parameter)
myResults$CI_LOW[28] <- myResults$ESTIMATION[28] + qt(p = 0.025, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[28]
myResults$CI_UP[28] <- myResults$ESTIMATION[28] + qt(p = 0.975, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[28]
myResults$t_VALUE[28] <- myApply$mean_Pi$statistic
myResults$p_VALUE[28] <- myApply$mean_Pi$p.value
myResults$FPR[28] <- funcFPR(t_value = myResults$t_VALUE[28], 
                             estimation = myResults$ESTIMATION[28], 
                             std_err = myResults$STD_ERR[28],
                             dFreedom = myApply$mean_Pi$parameter)

myResults$QUESTION[29] <- "Correlation M/M MEAN with mean_XATP in control - cd"
myResults$STATISTICS[29] <- "t_test df = 17"
myResults$ESTIMATION[29] <- myApply$mean_XATP$estimate
myResults$STD_ERR[29] <- sqrt((1 - myResults$ESTIMATION[29]^2)/myApply$mean_XATP$parameter)
myResults$CI_LOW[29] <- myResults$ESTIMATION[29] + qt(p = 0.025, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[29]
myResults$CI_UP[29] <- myResults$ESTIMATION[29] + qt(p = 0.975, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[29]
myResults$t_VALUE[29] <- myApply$mean_XATP$statistic
myResults$p_VALUE[29] <- myApply$mean_XATP$p.value
myResults$FPR[29] <- funcFPR(t_value = myResults$t_VALUE[29], 
                             estimation = myResults$ESTIMATION[29], 
                             std_err = myResults$STD_ERR[29],
                             dFreedom = myApply$mean_XATP$parameter)

rm(myApply)


#################### DD ########################################################

myApply <- apply(X = condData[condData$condition == "dd", 2:4], MARGIN = 2, 
                 FUN = function(a){
                 cor.test(x = a, 
                          y = paramData[paramData$condition == "diabete", c("MM_MEAN")], 
                          method = "pearson")
                          })

myResults[30:32, ] <- NA

myResults$QUESTION[30] <- "Correlation M/M MEAN with mean_AdpAtp in diabete - dd"
myResults$STATISTICS[30] <- "t_test df = 19"
myResults$ESTIMATION[30] <- myApply$mean_AdpAtp$estimate
myResults$STD_ERR[30] <- sqrt((1 - myResults$ESTIMATION[30]^2)/myApply$mean_AdpAtp$parameter)
myResults$CI_LOW[30] <- myResults$ESTIMATION[30] + qt(p = 0.025, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[30]
myResults$CI_UP[30] <- myResults$ESTIMATION[30] + qt(p = 0.975, df = myApply$mean_AdpAtp$parameter)*myResults$STD_ERR[30]
myResults$t_VALUE[30] <- myApply$mean_AdpAtp$statistic
myResults$p_VALUE[30] <- myApply$mean_AdpAtp$p.value
myResults$FPR[30] <- funcFPR(t_value = myResults$t_VALUE[30], 
                             estimation = myResults$ESTIMATION[30], 
                             std_err = myResults$STD_ERR[30],
                             dFreedom = myApply$mean_AdpAtp$parameter)

myResults$QUESTION[31] <- "Correlation M/M MEAN with mean_Pi in diabete - dd"
myResults$STATISTICS[31] <- "t_test df = 19"
myResults$ESTIMATION[31] <- myApply$mean_Pi$estimate
myResults$STD_ERR[31] <- sqrt((1 - myResults$ESTIMATION[31]^2)/myApply$mean_Pi$parameter)
myResults$CI_LOW[31] <- myResults$ESTIMATION[31] + qt(p = 0.025, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[31]
myResults$CI_UP[31] <- myResults$ESTIMATION[31] + qt(p = 0.975, df = myApply$mean_Pi$parameter)*myResults$STD_ERR[31]
myResults$t_VALUE[31] <- myApply$mean_Pi$statistic
myResults$p_VALUE[31] <- myApply$mean_Pi$p.value
myResults$FPR[31] <- funcFPR(t_value = myResults$t_VALUE[31], 
                             estimation = myResults$ESTIMATION[31], 
                             std_err = myResults$STD_ERR[31],
                             dFreedom = myApply$mean_Pi$parameter)

myResults$QUESTION[32] <- "Correlation M/M MEAN with mean_XATP in diabete - dd"
myResults$STATISTICS[32] <- "t_test df = 19"
myResults$ESTIMATION[32] <- myApply$mean_XATP$estimate
myResults$STD_ERR[32] <- sqrt((1 - myResults$ESTIMATION[32]^2)/myApply$mean_XATP$parameter)
myResults$CI_LOW[32] <- myResults$ESTIMATION[32] + qt(p = 0.025, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[32]
myResults$CI_UP[32] <- myResults$ESTIMATION[32] + qt(p = 0.975, df = myApply$mean_XATP$parameter)*myResults$STD_ERR[32]
myResults$t_VALUE[32] <- myApply$mean_XATP$statistic
myResults$p_VALUE[32] <- myApply$mean_XATP$p.value
myResults$FPR[32] <- funcFPR(t_value = myResults$t_VALUE[32], 
                             estimation = myResults$ESTIMATION[32], 
                             std_err = myResults$STD_ERR[32],
                             dFreedom = myApply$mean_XATP$parameter)

rm(myApply)

################################################################################
#q11-14: Does M/M MAD correlate with correlate  variables in CC, CD and DD     #
################################################################################

#################### CC ########################################################
myApply <- apply(X = condData[condData$condition == "cc", 5:7], MARGIN = 2, 
                 FUN = function(a){
                        cor.test(x = a, 
                                 y = paramData[paramData$condition == "control", c("MM_MAD")], 
                                 method = "pearson")
                                 })

myResults[33:35, ] <- NA

myResults$QUESTION[33] <- "Correlation M/M MAD with mad_AdpAtp in control - cc"
myResults$STATISTICS[33] <- "t_test df = 17"
myResults$ESTIMATION[33] <- myApply$mad_AdpAtp$estimate
myResults$STD_ERR[33] <- sqrt((1 - myResults$ESTIMATION[33]^2)/myApply$mad_AdpAtp$parameter)
myResults$CI_LOW[33] <- myResults$ESTIMATION[33] + qt(p = 0.025, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[33]
myResults$CI_UP[33] <- myResults$ESTIMATION[33] + qt(p = 0.975, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[33]
myResults$t_VALUE[33] <- myApply$mad_AdpAtp$statistic
myResults$p_VALUE[33] <- myApply$mad_AdpAtp$p.value
myResults$FPR[33] <- funcFPR(t_value = myResults$t_VALUE[33], 
                             estimation = myResults$ESTIMATION[33], 
                             std_err = myResults$STD_ERR[33],
                             dFreedom = myApply$mad_AdpAtp$parameter)

myResults$QUESTION[34] <- "Correlation M/M MAD with mad_Pi in control - cc"
myResults$STATISTICS[34] <- "t_test df = 17"
myResults$ESTIMATION[34] <- myApply$mad_Pi$estimate
myResults$STD_ERR[34] <- sqrt((1 - myResults$ESTIMATION[34]^2)/myApply$mad_Pi$parameter)
myResults$CI_LOW[34] <- myResults$ESTIMATION[34] + qt(p = 0.025, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[34]
myResults$CI_UP[34] <- myResults$ESTIMATION[34] + qt(p = 0.975, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[34]
myResults$t_VALUE[34] <- myApply$mad_Pi$statistic
myResults$p_VALUE[34] <- myApply$mad_Pi$p.value
myResults$FPR[34] <- funcFPR(t_value = myResults$t_VALUE[34], 
                             estimation = myResults$ESTIMATION[34], 
                             std_err = myResults$STD_ERR[34],
                             dFreedom = myApply$mad_Pi$parameter)

myResults$QUESTION[35] <- "Correlation M/M MAD with mad_XATP in control - cc"
myResults$STATISTICS[35] <- "t_test df = 17"
myResults$ESTIMATION[35] <- myApply$mad_XATP$estimate
myResults$STD_ERR[35] <- sqrt((1 - myResults$ESTIMATION[35]^2)/myApply$mad_XATP$parameter)
myResults$CI_LOW[35] <- myResults$ESTIMATION[35] + qt(p = 0.025, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[35]
myResults$CI_UP[35] <- myResults$ESTIMATION[35] + qt(p = 0.975, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[35]
myResults$t_VALUE[35] <- myApply$mad_XATP$statistic
myResults$p_VALUE[35] <- myApply$mad_XATP$p.value
myResults$FPR[35] <- funcFPR(t_value = myResults$t_VALUE[35], 
                             estimation = myResults$ESTIMATION[35], 
                             std_err = myResults$STD_ERR[35],
                             dFreedom = myApply$mad_XATP$parameter)

rm(myApply)



#################### CD ########################################################
myApply <- apply(X = condData[condData$condition == "cd", 5:7], MARGIN = 2, 
                 FUN = function(a){
                        cor.test(x = a, 
                                 y = paramData[paramData$condition == "control", c("MM_MAD")], 
                                 method = "pearson")
      })

myResults[36:38, ] <- NA

myResults$QUESTION[36] <- "Correlation M/M MAD with mad_AdpAtp in control - cd"
myResults$STATISTICS[36] <- "t_test df = 17"
myResults$ESTIMATION[36] <- myApply$mad_AdpAtp$estimate
myResults$STD_ERR[36] <- sqrt((1 - myResults$ESTIMATION[36]^2)/myApply$mad_AdpAtp$parameter)
myResults$CI_LOW[36] <- myResults$ESTIMATION[36] + qt(p = 0.025, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[36]
myResults$CI_UP[36] <- myResults$ESTIMATION[36] + qt(p = 0.975, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[36]
myResults$t_VALUE[36] <- myApply$mad_AdpAtp$statistic
myResults$p_VALUE[36] <- myApply$mad_AdpAtp$p.value
myResults$FPR[36] <- funcFPR(t_value = myResults$t_VALUE[36], 
                             estimation = myResults$ESTIMATION[36], 
                             std_err = myResults$STD_ERR[36],
                             dFreedom = myApply$mad_AdpAtp$parameter)

myResults$QUESTION[37] <- "Correlation M/M MAD with mad_Pi in control - cd"
myResults$STATISTICS[37] <- "t_test df = 17"
myResults$ESTIMATION[37] <- myApply$mad_Pi$estimate
myResults$STD_ERR[37] <- sqrt((1 - myResults$ESTIMATION[37]^2)/myApply$mad_Pi$parameter)
myResults$CI_LOW[37] <- myResults$ESTIMATION[37] + qt(p = 0.025, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[37]
myResults$CI_UP[37] <- myResults$ESTIMATION[37] + qt(p = 0.975, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[37]
myResults$t_VALUE[37] <- myApply$mad_Pi$statistic
myResults$p_VALUE[37] <- myApply$mad_Pi$p.value
myResults$FPR[37] <- funcFPR(t_value = myResults$t_VALUE[37], 
                             estimation = myResults$ESTIMATION[37], 
                             std_err = myResults$STD_ERR[37],
                             dFreedom = myApply$mad_Pi$parameter)

myResults$QUESTION[38] <- "Correlation M/M MAD with mad_XATP in control - cd"
myResults$STATISTICS[38] <- "t_test df = 17"
myResults$ESTIMATION[38] <- myApply$mad_XATP$estimate
myResults$STD_ERR[38] <- sqrt((1 - myResults$ESTIMATION[38]^2)/myApply$mad_XATP$parameter)
myResults$CI_LOW[38] <- myResults$ESTIMATION[38] + qt(p = 0.025, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[38]
myResults$CI_UP[38] <- myResults$ESTIMATION[38] + qt(p = 0.975, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[38]
myResults$t_VALUE[38] <- myApply$mad_XATP$statistic
myResults$p_VALUE[38] <- myApply$mad_XATP$p.value
myResults$FPR[38] <- funcFPR(t_value = myResults$t_VALUE[38], 
                             estimation = myResults$ESTIMATION[38], 
                             std_err = myResults$STD_ERR[38],
                             dFreedom = myApply$mad_XATP$parameter)

rm(myApply)
#################### DD ########################################################

myApply <- apply(X = condData[condData$condition == "dd", 5:7], MARGIN = 2, 
                 FUN = function(a){
                        cor.test(x = a, 
                        y = paramData[paramData$condition == "diabete", c("MM_MAD")], 
                        method = "pearson")
                        })

myResults[39:41, ] <- NA

myResults$QUESTION[39] <- "Correlation M/M MAD with mad_AdpAtp in diabete - dd"
myResults$STATISTICS[39] <- "t_test df = 19"
myResults$ESTIMATION[39] <- myApply$mad_AdpAtp$estimate
myResults$STD_ERR[39] <- sqrt((1 - myResults$ESTIMATION[39]^2)/myApply$mad_AdpAtp$parameter)
myResults$CI_LOW[39] <- myResults$ESTIMATION[39] + qt(p = 0.025, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[39]
myResults$CI_UP[39] <- myResults$ESTIMATION[39] + qt(p = 0.975, df = myApply$mad_AdpAtp$parameter)*myResults$STD_ERR[39]
myResults$t_VALUE[39] <- myApply$mad_AdpAtp$statistic
myResults$p_VALUE[39] <- myApply$mad_AdpAtp$p.value
myResults$FPR[39] <- funcFPR(t_value = myResults$t_VALUE[39], 
                             estimation = myResults$ESTIMATION[39], 
                             std_err = myResults$STD_ERR[39],
                             dFreedom = myApply$mad_AdpAtp$parameter)

myResults$QUESTION[40] <- "Correlation M/M MAD with mad_Pi in diabete - dd"
myResults$STATISTICS[40] <- "t_test df = 19"
myResults$ESTIMATION[40] <- myApply$mad_Pi$estimate
myResults$STD_ERR[40] <- sqrt((1 - myResults$ESTIMATION[40]^2)/myApply$mad_Pi$parameter)
myResults$CI_LOW[40] <- myResults$ESTIMATION[40] + qt(p = 0.025, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[40]
myResults$CI_UP[40] <- myResults$ESTIMATION[40] + qt(p = 0.975, df = myApply$mad_Pi$parameter)*myResults$STD_ERR[40]
myResults$t_VALUE[40] <- myApply$mad_Pi$statistic
myResults$p_VALUE[40] <- myApply$mad_Pi$p.value
myResults$FPR[40] <- funcFPR(t_value = myResults$t_VALUE[40], 
                             estimation = myResults$ESTIMATION[40], 
                             std_err = myResults$STD_ERR[40],
                             dFreedom = myApply$mad_Pi$parameter)

myResults$QUESTION[41] <- "Correlation M/M MAD with mad_XATP in diabete - dd"
myResults$STATISTICS[41] <- "t_test df = 19"
myResults$ESTIMATION[41] <- myApply$mad_XATP$estimate
myResults$STD_ERR[41] <- sqrt((1 - myResults$ESTIMATION[41]^2)/myApply$mad_XATP$parameter)
myResults$CI_LOW[41] <- myResults$ESTIMATION[41] + qt(p = 0.025, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[41]
myResults$CI_UP[41] <- myResults$ESTIMATION[41] + qt(p = 0.975, df = myApply$mad_XATP$parameter)*myResults$STD_ERR[41]
myResults$t_VALUE[41] <- myApply$mad_XATP$statistic
myResults$p_VALUE[41] <- myApply$mad_XATP$p.value
myResults$FPR[41] <- funcFPR(t_value = myResults$t_VALUE[41], 
                             estimation = myResults$ESTIMATION[41], 
                             std_err = myResults$STD_ERR[41],
                             dFreedom = myApply$mad_XATP$parameter)

rm(myApply)


write.csv(x = myResults, file = "computationalResults.csv")
