#2013 Dec 9, check consistency of last week's gnls fitting results

rm(list=ls())
setwd("~/github/LOH_H2O2_2012-master/analysis")
files = list.files(path=".", pattern="*csv")
files

#2013 manual fitting results
#tbManu = read.csv("H2O2_Log_Plot_Summarized_data,2013May30.csv")
#tbManu = read.csv("../H2O2_Log_Plot_Summarized_data,2013Dec4_v2.csv") # new half-black results
#names(tbManu)=c("file", "Date","Strain","CvManu","CbManu", "Cb0.5Manu", "OD", "note13", "repeat")

#2013 auto fitting results with manual results
tb = read.csv("../_merged.tb.20131206_14,49pm.csv")
tb$file = as.character(tbAuto$files)


#compare Cv
tb$CvChanges = tb$CvManu / tb$Cv
summary(tb$CvChanges)
hist(tb$CvChanges,br=10)
summary(lm(tb$CvManu ~ tb$Cv)) #R2=0.14, ??
tb[tb$CvChanges>1.3, c("file", "CvManu", "Cv", "CvChanges")]
tb[tb$CvChanges<1/1.3, c("file", "CvManu", "Cv", "CvChanges")]

tb[grep("sod",tb$file), c("file", "CvManu", "Cv", "CvChanges") ] #sod, auto fit makes sense
tb[grep("SGU",tb$file), c("file", "CvManu", "Cv", "CvChanges") ] #SGU, auto fit seems better
tb[grep("M2-8",tb$file), c("file", "CvManu", "Cv", "CvChanges") ] #unchanged
tb[grep("M32",tb$file), c("file", "CvManu", "Cv", "CvChanges") ] #unchanged
#conclusions, batch Cv and manual Cv are similar

#compare Cb
summary(lm(tb$CvManu ~ tb$Cb)) #R2=0.79
tb$CbChanges = tb$CbManu / tb$Cb
summary(tb$CbChanges)
hist(tb$CbChanges, br=20)
tb[tb$CbChanges>2, c("file", "CvManu", "Cv", "CvChanges", "CbManu", "Cb", "CbChanges")]
tb[tb$CbChanges<1/2, c("file", "CvManu", "Cv", "CvChanges", "CbManu", "Cb", "CbChanges")]
#2013Dec4
# These differences likely occur because batch fitting did not take larger errors in some datapoints into account
# I should to use 'weights' in gnls. I could use the 1/ total_number_of_colonies as weights. 
# 

# choose the manual fitting for these expts
#tb$Cb[39] = tb$CbManu[39] #YPS163,20120503,YPS163-H2O2LOH.csv

tb$Cv.vs.Cb.Manu = tb$CvManu / tb$CbManu   
tb$Cv.vs.Cb.gnls = tb$Cv / tb$Cb

summary(lm(tb$Cv.vs.Cb.Manu ~ tb$Cv.vs.Cb.gnls)) #R=0.39
summary(lm(tb$CvManu/tb$Cb ~ tb$CvManu/tb$CbManu)) #R=0.32
summary(lm(tb$Cv/tb$Cb ~ tb$Cv/tb$CbManu)) #R=0.44

tb2= tb;
tb2$Cv = ifelse( tb2$CvChanges>1.3 |tb2$CvChanges< 1/1.3 , tb2$CvManu, tb2$Cv  )
tb2$Cb = ifelse( tb2$CbChanges>1.3 |tb2$CbChanges< 1/1.3 , tb2$CbManu, tb2$Cb )

summary(lm(tb2$Cv / tb2$Cb ~ tb2$CvManu / tb2$CbManu)) #R=0.46

write.csv(tb, "_merged.tb.20131209_9,15am.csv", row.names=F)
