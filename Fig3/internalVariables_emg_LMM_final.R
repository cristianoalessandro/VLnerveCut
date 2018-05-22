# Refer to file model_J5_log_randDay_nested.RData inside the model folder

require(nlme)
require(car)
require(graphics)

library(multcomp)
library(reshape2)

setwd("./")

# empty environment
rm(list=ls())

# load data
data_all = read.csv("./bursts_RF_VM.csv")

############################### DESIRED VALUES ###############################

# desired rats
desRat <- as.factor( c("rat21","rat22","rat24","rat25","J5","J9") )

# desired days
desDays <- as.factor( c("-1","1","9","45") )

#desired muscles
desMuscle <- as.factor( c("VM","RF") )

############################### EXTRACT DESIRED VALUES ###############################

data_emg <-data_all

# transform days into a factor
data_emg$des_days <- as.factor(data_emg$des_days)

# extract desired muscles, days and subjects
data_emg <- data_emg[data_emg$ratID %in% desRat & data_emg$des_days %in% desDays & data_emg$muscle %in% desMuscle,]

# remove undesired muscles, days and subjects from levels
data_emg$des_days <- droplevels(data_emg$des_days)
data_emg$ratID    <- droplevels(data_emg$ratID)
data_emg$muscle   <- droplevels(data_emg$muscle)


######################## FIT MODEL ########################

lCtr <- lmeControl(maxIter = 1000, msMaxIter=50, niterEM = 100, msMaxEval=200, msVerbose = FALSE, opt = 'optim')

data_emg_trf <- data_emg
data_emg_trf$data = log(data_emg_trf$data)

hist(data_emg_trf$data,breaks=10)
boxplot(data_emg_trf$data)

# Models
emg.lme <- lme(data ~ des_days*muscle, random = ~des_days|ratID/cycle, data=data_emg_trf, na.action=na.omit, method = "ML", control=lCtr )

######################## ...or load model already fitted

load('model_J5_log_randDay_nested.RData')
emg.lme = emg.lme8;

############################ assumptions ############################

emg.lm <- emg.lme
summary(emg.lm)

#### Fitted values

plot(emg.lm, data~fitted(.)|ratID*muscle, abline = c(0,1))
plot(emg.lm, data~fitted(.)|muscle*des_days, abline = c(0,1))
plot(emg.lm, data~fitted(.), abline = c(0,1)) 

#### Residuals (should have constant variance and be normally distributed)

plot(emg.lm)                       # Residuals (should have constant variance)
qqp(residuals(emg.lm),"norm")      # Residuals (should be normally distributed)

myx <- seq(min(residuals(emg.lm)), max(residuals(emg.lm)), length.out= 100)
myhist <- hist(residuals(emg.lm),100)
multiplier <- myhist$counts / myhist$density
normal <- dnorm(x = myx, mean = 0, sd = emg.lm$sigma)

lines(myx, normal* multiplier[1], col = "blue", lwd = 2)

#### Random effects (should be normally distributed)

randBase <- ranef(emg.lm)$ratID

# baseline
randEff <- randBase$`(Intercept)`
qqp( randEff )                         # Random effects (should be notmally distributed)

myx <- seq(min(randEff), max(randEff), length.out= 100)
myhist <- hist(randEff,5)
shapiro.test(randEff)

# day 1
randEff <- randBase$des_days1
qqp( randEff )  

# day 9
randEff <- randBase$des_days9
qqp( randEff )  

# day 45
randEff <- randBase$des_days45
qqp( randEff )  

################# cycle

randBase <- ranef(emg.lm)$cycle

# baseline
randEff <- randBase$`(Intercept)`
qqp( randEff )                         # Random effects (should be notmally distributed)
hist(randEff)

# day 1
randEff <- randBase$des_days1
qqp( randEff )  
hist(randEff,20)

# day 9
randEff <- randBase$des_days9
qqp( randEff )  
hist(randEff,20)

# day 45
randEff <- randBase$des_days45
qqp( randEff )  
hist(randEff,20)


#############################################################################
################################# Post-hoc ##################################
#############################################################################

############### two-sided tests for RF-vs-VM, together with single VM and RF tests, corrected with Bonferroni

ph_conditional <- c("des_days1  + des_days1:muscleVM  = 0",  # VM, day1=bas
                    "des_days9  + des_days9:muscleVM  = 0",  #     day9=bas
                    "des_days45 + des_days45:muscleVM = 0",  #     day45=bas
                    "des_days1 = 0",                         # RF, day1=bas
                    "des_days9 = 0",                         #     day9=bas   
                    "des_days45 = 0",                        #     day45=bas
                    "muscleVM   + des_days1:muscleVM  = 0",  # RF-vs-VM, day1
                    "muscleVM   + des_days9:muscleVM  = 0",  #           day9
                    "muscleVM   + des_days45:muscleVM = 0"); #           day45

emg.ph <- glht(emg.lm, linfct = ph_conditional);

summary(emg.ph, test = adjusted("bonferroni"))
plot(emg.ph)
