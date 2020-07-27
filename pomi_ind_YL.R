#Impute the POMI+IND method using the MICE package in R 
# Description about the POMI+IND method is described in the manuscript entitled "Using multiple imputation to classify potential outcomes subgroups" 
# by Yun Li, Irina Bondarenko, Michael R. Elliott, Timothy P. Hofer and Jeremy M.G. Taylor  

#  c13: chemotherapy, the observed outcome Y
#  chemo0 = Y_0, chemo1 = Y_1  
#  rs: recurrence score assay testing status: 0, 1 
#  ghi_recurrence_scoren: recurrence score, a post-exposure variable 
#  ghi: rescurrence score assay results: low, intermediate, and high score 
#  age_survey: at survey 
#  comrb: the number of major comorbid conditions 
#  size: tumor size 
#  ncode: cancer stage 
#  grade: tumor grade 
#  race: race 
#  high-risk: genetic risk 
#  insgroup: insurance group 
#  education: education level 
#  m6: menstrual period status 
#  k5: income 
#  site27: site A or B status 

library(lattice)
library(sas7bdat)
library(mice)
library(haven)
utest0 <- read_sas("U:/Biostat/forR.sas7bdat", 
                 NULL)
#define categorical variables
summary(utest0);
utest0$m6 = factor(utest0$m6, levels = c(1,2,3));
utest0$grade = factor(utest0$grade, levels = c(1,2,3));
utest0$size = factor(utest0$size, levels = c(1,2,3));
utest0$ncode = factor(utest0$ncode, levels = c(0,1,2));
utest0$comrb = factor(utest0$comrb, levels = c(0,1,2))
utest0$race = factor(utest0$race, levels = c(1,2,3,4,5));
utest0$insgroup = factor(utest0$insgroup, levels =c(0,1,2,3,4));
utest0$education2 = factor(utest0$education2, levels = c(-1,0,1));
utest0$site27 = factor(utest0$site27, levels = 0:1);
utest0$high_risk = factor(utest0$high_risk, levels = 0:1);
utest0$chemo1 = factor(utest0$chemo1, levels = 0:1);
utest0$chemo0 = factor(utest0$chemo0, levels = 0:1);
utest0$rs = factor(utest0$rs, levels = 0:1);

summary(utest0)

myvars<-c("age_survey",   "comrb", "chemo0", "chemo1", "size", "ncode",
          "grade", "race" , "high_risk",  "insgroup", "education2", "m6", "k5",
          "ghi", "rs", "site27");
          
utest <- utest0[myvars]

summary(utest)

# prediction matrix
predT3=make.predictorMatrix(utest)   

#E.g., predT3["chemo1","chemo0"]=0 means that "chemo0" is excluded in the prediction of "chemo1" 
#The default is that all variables are included in the imputation model  

predT3["chemo1","chemo0"]=0
predT3["chemo0","chemo1"]=0
predT3["chemo1","rs"]=0
predT3["chemo0","rs"]=0
predT3["ghi","rs"]=0
predT3["chemo0","ghi"]=0
predT3["ghi","chemo0"]=0

predT3

tempdataT3<-mice(data = utest, m =5 , pred=predT3, maxit =20, seed = 500)
tempdataT3$warnings
tempdataT3$loggedEvents;
compT3<-complete(tempdataT3, action="long" )
dim(compT3)

compT1[] <- lapply(compT3, function(x) as.numeric(as.character(x)))
summary(compT3)
write.csv(compT3,file="U:/Biostat/fromR.csv")