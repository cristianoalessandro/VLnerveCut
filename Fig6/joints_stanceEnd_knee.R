require(nlme)
require(graphics)
require(car)
require(MASS)

library(multcomp)

setwd("./")

# empty environment
rm(list=ls())

# load data
data_all = read.csv("./joints_stanceEnd.csv")

####################################################################################

# desired joints
desJointsName <- "knee"
desJoints     <- as.factor(desJointsName);

# desired days
desDays     <- as.factor(c("-1","1","14","48"));

# desired subjects
desRatID <- as.factor(c("rat21","rat22","rat24","rat25","J5","J9"))

# extract data
data_red <- data_all[data_all$ratID %in% desRatID & data_all$joint %in% desJoints & data_all$des_days %in% desDays,]

data_red$ratID    <- factor(data_red$ratID)
data_red$joint    <- factor(data_red$joint)

data_red$des_days <- as.factor(data_red$des_days)

# Data transformation
data_red_trf      <- data_red

boxplot(data_red_trf$data)
hist(data_red_trf$data,breaks=10)

############################### Fit model ##########################################

lCtr <- lmeControl(maxIter = 5000, niterEM = 500, msMaxIter=1000, msMaxEval=1000, msVerbose = FALSE, opt = 'optim')

# Days (and intercept) as a random variable (to account for the variability cross subjects)
lev.lm_nonNest_r  <- lme(data ~ des_days, random = ~des_days|ratID, data=data_red_trf, na.action=na.omit, method = "ML", control=lCtr)

lev.lm <- lev.lm_nonNest_r
summary(lev.lm)

############################ assumptions ############################

#### Fitted values

plot(lev.lm, data~fitted(.), abline = c(0,1)) 
plot(lev.lm, data~fitted(.)|ratID, abline = c(0,1))
plot(lev.lm, data~fitted(.)|des_days, abline = c(0,1))
plot(lev.lm, data~fitted(.)|ratID*des_days, abline = c(0,1))


#### Residuals (should have constant variance and be normally distributed)

plot(lev.lm)                       # Residuals (should have constant variance)
qqp(residuals(lev.lm),"norm")      # Residuals (should be normally distributed)

boxplot(residuals(lev.lm))

myx <- seq(min(residuals(lev.lm)), max(residuals(lev.lm)), length.out= 100)
myhist <- hist(residuals(lev.lm),breaks=20)
multiplier <- myhist$counts / myhist$density

varMean <- lev.lm$sigma * mean(fitted(lev.lm)) ^ (2*lev.lm$modelStruct$varStruct)
normal <- dnorm(x = myx, mean = 0, sd = sqrt(varMean))

normal <- dnorm(x = myx, mean = 0, sd = lev.lm$sigma)

lines(myx, normal* multiplier[1], col = "blue", lwd = 2)

shapiro.test(residuals(lev.lm))


#### Random effects (should be normally distributed)

randBase <- ranef(lev.lm)

randEff <- randBase$`(Intercept)`
qqp( randEff )                         # Random effects (should be notmally distributed)
myx <- seq(min(randEff), max(randEff), length.out= 100)
myhist <- hist(randEff)
multiplier <- myhist$counts / myhist$density
normal <- dnorm(x = myx, mean = 0, sd = as.numeric(VarCorr(lev.lm)[1,2]) )
lines(myx, normal*multiplier[1], col = "blue", lwd = 2)

shapiro.test(randEff)

randEff <- randBase$des_days1
qqp( randEff )  

randEff <- randBase$des_days14
qqp( randEff )  

randEff <- randBase$des_days48
qqp( randEff )  

#############################################################################
################################ PREDICTIONS ################################
#############################################################################

# unique values of the factors
pred <- expand.grid(joint=unique(data_red$joint),des_days=unique(data_red$des_days))

# polulation level (level=0) predictions
y=predict(lev.lm, level=0, newdata=pred)
pred$pred_val = y

pred <- pred[order(pred$joint, pred$des_days),]

############ Compute menas ############
# remove non-used columns and NaN
dat_lev.small <- data_red_trf[,c("ratID","joint","des_days","data")]
dat_lev.small <- dat_lev.small[!is.na(dat_lev.small$data),]

# mean over gait cycle and rats
data_red.mean <- aggregate(dat_lev.small$data, by = list(dat_lev.small$joint, dat_lev.small$des_days), FUN = 'mean')
colnames(data_red.mean) <- c("joint","des_days","data_mean")

data_red.mean <- data_red.mean[order(data_red.mean$joint, data_red.mean$des_days),]

pred$data_mean = data_red.mean$data_mean
pred$diff      = pred$pred_val - pred$data_mean
pred$diff_perc = pred$diff/pred$data_mean*100

# mean over gait cycles for each rat
data_red.meanSbj <- aggregate(dat_lev.small$data, by = list(dat_lev.small$ratID, dat_lev.small$joint, dat_lev.small$des_days), FUN = 'mean')
colnames(data_red.meanSbj) <- c("ratID","joint","des_days","data_mean")

# Mean of the means
data_red.meanSbjAll <- aggregate(data_red.meanSbj$data_mean, by = list(data_red.meanSbj$joint, data_red.meanSbj$des_days), FUN = 'mean')
colnames(data_red.meanSbjAll) <- c("joint","des_days","data_mean")

pred$data_mean2 = data_red.meanSbjAll$data_mean

pred$diff2      = pred$pred_val - pred$data_mean2
pred$diff_perc2 = pred$diff2/pred$data_mean2*100


###############################################################################
############################# PREDICTIONS SUBJECTS ############################
###############################################################################

# unique values of the factors
pred.sbj <- expand.grid(ratID=unique(data_red$ratID), joint=unique(data_red$joint), des_days=unique(data_red$des_days))

# polulation level (level=0) predictions
y=predict(lev.lm, level=1, newdata=pred.sbj)
pred.sbj$pred_val = y

pred.sbj <- pred.sbj[order(pred.sbj$ratID, pred.sbj$joint, pred.sbj$des_days),]

############ Compute menas ############
# remove non-used columns and NaN
dat_lev.small <- data_red_trf[,c("ratID","joint","des_days","data")]
dat_lev.small <- dat_lev.small[!is.na(dat_lev.small$data),]

# mean over gait cycles
data_red.mean <- aggregate(dat_lev.small$data, by = list(dat_lev.small$ratID, dat_lev.small$joint, dat_lev.small$des_days), FUN = 'mean')
colnames(data_red.mean) <- c("ratID","joint","des_days","data_mean")

data_red.mean <- data_red.mean[order(data_red.mean$ratID, data_red.mean$joint, data_red.mean$des_days),]

pred.sbj$data_mean <- NA
for(i in 1:dim(pred.sbj)[1]){
  for(j in 1:dim(data_red.mean)[1]){
    if( all(pred.sbj[i,1:3]==data_red.mean[j,1:3]) ){
      pred.sbj[i,]$data_mean<-data_red.mean[j,]$data_mean
    }
  }
}
pred.sbj$data_mean <- as.numeric(pred.sbj$data_mean)

pred.sbj$diff      = pred.sbj$pred_val - pred.sbj$data_mean
pred.sbj$diff_perc = pred.sbj$diff/pred.sbj$data_mean*100


############################ Post-hoc (time as factor) ############################

ph_conditional <- c("des_days1  = 0",
                    "des_days14 = 0",
                    "des_days48 = 0");

lev.ph <- glht(lev.lm, linfct = ph_conditional);
summary(lev.ph, test=adjusted("bonferroni"))
