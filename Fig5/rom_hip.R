require(nlme)
require(graphics)
require(car)
require(MASS)

library(multcomp)

setwd("./")

# empty environment
rm(list=ls())

# load data
data_all = read.csv("./rom.csv")

####################################################################################

# desired joints
desJointsName <- "hip"
desJoints     <- as.factor(desJointsName);

# desired days
desDays     <- as.factor(c("-1","1","9","45"));

# For hip
desRatID <- as.factor(c("rat22","rat24","rat25","J5","J9"))

# extract data
data_rom <- data_all[data_all$ratID %in% desRatID & data_all$joint %in% desJoints & data_all$des_days %in% desDays,]

data_rom$ratID    <- factor(data_rom$ratID)
data_rom$joint    <- factor(data_rom$joint)

data_rom$des_days <- as.factor(data_rom$des_days)

# Data transformation
data_rom_trf <- data_rom

boxplot(data_rom_trf$data)
hist(data_rom_trf$data,breaks=10)

############################### Fit model ##########################################

lCtr <- lmeControl(maxIter = 5000, niterEM = 500, msMaxIter=1000, msMaxEval=1000, msVerbose = FALSE, opt = 'optim')

# Days (and intercept) as a random variable (to account for the variability cross subjects)
rom.lm_nonNest_r  <- lme(data ~ des_days, random = ~des_days|ratID, data=data_rom_trf, na.action=na.omit, method = "ML", control=lCtr)

rom.lm <- rom.lm_nonNest_r
summary(rom.lm)

############################ assumptions ############################

#### Fitted values

plot(rom.lm, data~fitted(.), abline = c(0,1)) 
plot(rom.lm, data~fitted(.)|ratID, abline = c(0,1))
plot(rom.lm, data~fitted(.)|des_days, abline = c(0,1))
plot(rom.lm, data~fitted(.)|ratID*des_days, abline = c(0,1))


#### Residuals (should have constant variance and be normally distributed)

plot(rom.lm)                       # Residuals (should have constant variance)
qqp(residuals(rom.lm),"norm")      # Residuals (should be normally distributed)

boxplot(residuals(rom.lm))

myx <- seq(min(residuals(rom.lm)), max(residuals(rom.lm)), length.out= 100)
myhist <- hist(residuals(rom.lm),breaks=20)
multiplier <- myhist$counts / myhist$density

varMean <- rom.lm$sigma * mean(fitted(rom.lm)) ^ (2*rom.lm$modelStruct$varStruct)
normal <- dnorm(x = myx, mean = 0, sd = sqrt(varMean))

normal <- dnorm(x = myx, mean = 0, sd = rom.lm$sigma)

lines(myx, normal* multiplier[1], col = "blue", lwd = 2)

shapiro.test(residuals(rom.lm))


#### Random effects (should be normally distributed)

randBase <- ranef(rom.lm)

randEff <- randBase$`(Intercept)`
qqp( randEff )                         # Random effects (should be notmally distributed)

myx <- seq(min(randEff), max(randEff), length.out= 100)
myhist <- hist(randEff)
multiplier <- myhist$counts / myhist$density
normal <- dnorm(x = myx, mean = 0, sd = as.numeric(VarCorr(rom.lm)[1,2]) )
lines(myx, normal*multiplier[1], col = "blue", lwd = 2)

shapiro.test(randEff)

randEff <- randBase$des_days1
qqp( randEff )  

randEff <- randBase$des_days9
qqp( randEff )  

randEff <- randBase$des_days45
qqp( randEff )  

#############################################################################
################################ PREDICTIONS ################################
#############################################################################

# unique values of the factors
pred <- expand.grid(joint=unique(data_rom$joint),des_days=unique(data_rom$des_days))

# polulation level (level=0) predictions
y=predict(rom.lm, level=0, newdata=pred)
pred$pred_val = y

pred <- pred[order(pred$joint, pred$des_days),]

############ Compute menas ############
# remove non-used columns and NaN
dat_rom.small <- data_rom_trf[,c("ratID","joint","des_days","data")]
dat_rom.small <- dat_rom.small[!is.na(dat_rom.small$data),]

# mean over gait cycle and rats
data_rom.mean <- aggregate(dat_rom.small$data, by = list(dat_rom.small$joint, dat_rom.small$des_days), FUN = 'mean')
colnames(data_rom.mean) <- c("joint","des_days","data_mean")

data_rom.mean <- data_rom.mean[order(data_rom.mean$joint, data_rom.mean$des_days),]

pred$data_mean = data_rom.mean$data_mean
pred$diff      = pred$pred_val - pred$data_mean
pred$diff_perc = pred$diff/pred$data_mean*100

# mean over gait cycles for each rat
data_rom.meanSbj <- aggregate(dat_rom.small$data, by = list(dat_rom.small$ratID,dat_rom.small$joint, dat_rom.small$des_days), FUN = 'mean')
colnames(data_rom.meanSbj) <- c("ratID","joint","des_days","data_mean")

# Mean of the means
data_rom.meanSbjAll <- aggregate(data_rom.meanSbj$data_mean, by = list(data_rom.meanSbj$joint, data_rom.meanSbj$des_days), FUN = 'mean')
colnames(data_rom.meanSbjAll) <- c("joint","des_days","data_mean")

pred$data_mean2 = data_rom.meanSbjAll$data_mean

pred$diff2      = pred$pred_val - pred$data_mean2
pred$diff_perc2 = pred$diff2/pred$data_mean*100

###############################################################################
############################# PREDICTIONS SUBJECTS ############################
###############################################################################

# unique values of the factors
pred.sbj <- expand.grid(ratID=unique(data_rom$ratID), joint=unique(data_rom$joint), des_days=unique(data_rom$des_days))

# polulation level (level=0) predictions
y=predict(rom.lm, level=1, newdata=pred.sbj)
pred.sbj$pred_val = y

pred.sbj <- pred.sbj[order(pred.sbj$ratID, pred.sbj$joint, pred.sbj$des_days),]

############ Compute menas ############
# remove non-used columns and NaN
dat_rom.small <- data_rom_trf[,c("ratID","joint","des_days","data")]
dat_rom.small <- dat_rom.small[!is.na(dat_rom.small$data),]

# mean over gait cycles
data_rom.mean <- aggregate(dat_rom.small$data, by = list(dat_rom.small$ratID,dat_rom.small$joint, dat_rom.small$des_days), FUN = 'mean')
colnames(data_rom.mean) <- c("ratID","joint","des_days","data_mean")

data_rom.mean <- data_rom.mean[order(data_rom.mean$ratID, data_rom.mean$joint, data_rom.mean$des_days),]

pred.sbj$data_mean <- NA
for(i in 1:dim(pred.sbj)[1]){
  for(j in 1:dim(data_rom.mean)[1]){
    if( all(pred.sbj[i,1:3]==data_rom.mean[j,1:3]) ){
      pred.sbj[i,]$data_mean<-data_rom.mean[j,]$data_mean
    }
  }
}
pred.sbj$data_mean <- as.numeric(pred.sbj$data_mean)

pred.sbj$diff      = pred.sbj$pred_val - pred.sbj$data_mean
pred.sbj$diff_perc = pred.sbj$diff/pred.sbj$data_mean*100

############################ Post-hoc (time as factor) ############################

ph_conditional <- c("des_days1  = 0",
                    "des_days9  = 0",
                    "des_days45 = 0");

rom.ph <- glht(rom.lm, linfct = ph_conditional);
summary(rom.ph, test=adjusted("bonferroni"))

# RESULTS:

# - day1 significantly different from day-1
# - day9 significantly different from day-1
# - day45 NOT significantly different from day-1

