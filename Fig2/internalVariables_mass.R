#require(nlme)
require(car)
#require(graphics)
#library(multcomp)
#ibrary(reshape2)

setwd("./")

# empty environment
rm(list=ls())

# load data
data_all = read.csv("mass_quads.csv")

############################### DESIRED VALUES ###############################

# desired rats
desRat <- as.factor( c("Rat21","Rat22","Rat24","Rat25","J5","J9") )

#desired muscles
desMuscle <- as.factor( c("VL","VM","VI","RF") )

############################### EXTRACT DESIRED VALUES ###############################

data_mass <- data_all
data_mass <- data_mass[data_mass$ratId %in% desRat & data_mass$muscle %in% desMuscle,]

############################### VL ###############################
data_VL <- data_mass[data_mass$muscle=="VL",]
data_tmp <- data_VL

#desRat   <- as.factor( c("Rat21","Rat22","Rat24","Rat25","J5") )
#data_tmp <- data_tmp[data_tmp$ratId %in% desRat,]

# Check data
colMeans( cbind(data_tmp$ipsi,data_tmp$contra) )

# Are the differences normally distributed?
qqp          ( data_tmp$ipsi - data_tmp$contra )
shapiro.test ( data_tmp$ipsi - data_tmp$contra )

VL_t <- t.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
VL_t

############################### VM ###############################
data_VM  <- data_mass[data_mass$muscle=="VM",]
data_tmp <- data_VM

# Check data
colMeans( cbind(data_tmp$ipsi,data_tmp$contra) )

# Are the differences normally distributed?
qqp          ( data_tmp$ipsi - data_tmp$contra )
shapiro.test ( data_tmp$ipsi - data_tmp$contra )

VM_t <- t.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
VM_t

############################### VI ###############################
data_VI  <- data_mass[data_mass$muscle=="VI",]
data_tmp <- data_VI

# Remove outliers
#desRat   <- as.factor( c("Rat21","Rat22","Rat24","Rat25","J5") )
#data_tmp <- data_tmp[data_tmp$ratId %in% desRat,]

# Check data
colMeans( cbind(data_tmp$ipsi,data_tmp$contra) )

# Are the differences normally distributed?
qqp          ( data_tmp$ipsi - data_tmp$contra )
shapiro.test ( data_tmp$ipsi - data_tmp$contra )

VI_t <- t.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
#VI_t <- wilcox.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
VI_t

############################### RF ###############################
data_RF  <- data_mass[data_mass$muscle=="RF",]
data_tmp <- data_RF

# Remove outliers
#desRat   <- as.factor( c("Rat22","Rat24","Rat25","J5","J9") )
#data_tmp <- data_tmp[data_tmp$ratId %in% desRat,]

# Check data
colMeans( cbind(data_tmp$ipsi,data_tmp$contra) )

# Are the differences normally distributed?
qqp          ( data_tmp$ipsi - data_tmp$contra )
shapiro.test ( data_tmp$ipsi - data_tmp$contra )

RF_t <- t.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
#RF_t <- wilcox.test( data_tmp$ipsi, data_tmp$contra, paired=TRUE)
RF_t

################# Correction on p-values ######################

p_vec     <- c( VL_t$p.value , VM_t$p.value , VI_t$p.value , RF_t$p.value )
p_vec_adj <- p.adjust(p_vec , method="bonferroni")

muscle  <- c( "VL","VM","VI","RF" )
results <- data.frame( muscle, p_vec, p_vec_adj )

