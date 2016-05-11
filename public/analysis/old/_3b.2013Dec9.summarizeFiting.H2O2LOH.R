# 2013Dec 10, add Cb0.5.vs.Cv to output

# summarize H2O2-LOH data from _2a2.2013Dec5.stepBystep.gnslWeights.R

rm(list=ls())

setwd("~/github/LOH_H2O2_2012-master/analysis")
list.files(, pattern="_merged")

#read step-by-step fitting results, using _2a2.2013Dec5.stepBystep.gnslWeights.R
tbH2O2 = read.csv("_merged.tb.20131209_9am.csv")
str(tbH2O2)

#find strains, gsub
tbH2O2$strains = gsub("[\\,|\\.].*", "", tbH2O2$file)

strains2 = unique(tbH2O2$strains)
table(tbH2O2$strains)

tbH2O2$Cb.vs.CvManu = tbH2O2$CbManu / tbH2O2$CvManu
tbH2O2$Cb.vs.Cv = tbH2O2$Cb / tbH2O2$Cv
tbH2O2$Cv.vs.Cb = tbH2O2$Cv / tbH2O2$Cb
tbH2O2$Cb0.5.vs.Cv = tbH2O2$Cb0.5 / tbH2O2$Cv
tbH2O2$half.vs.full = tbH2O2$Cb0.5 / tbH2O2$Cb
str( tbH2O2  )

tbH2O2Report = data.frame(strains2)
names( tbH2O2Report ) = c("strains")
# strain = 'BY4743' i =1
for( i in 1:length(strains2) ) {  
  tmp = tbH2O2[tbH2O2$strains == strains2[i],  ]  
  
  tbH2O2Report$CvManu[i] = mean(tmp$CvManu, na.rm=T)
  tbH2O2Report$CvManuSTD[i] = sd(tmp$CvManu, na.rm=T)    
  tbH2O2Report$CvMean[i] = mean(tmp$Cv, na.rm=T)
  tbH2O2Report$CvSTD[i] = sd(tmp$Cv, na.rm=T)  
  
  tbH2O2Report$CbManu[i] = mean(tmp$CbManu, na.rm=T)
  tbH2O2Report$CbMean[i] = mean(tmp$Cb, na.rm=T)
  tbH2O2Report$CbSTD[i] = sd(tmp$Cb, na.rm=T)
  tbH2O2Report$bminH2O2[i] = mean(tmp$b.min, na.rm=T)
  tbH2O2Report$bmaxH2O2[i] = mean(tmp$b.max, na.rm=T)
      
  tmp$Cv.vs.Cb = tmp$Cv / tmp$CbManu   ###2013 Dec 9  change !!!!!!!!!
  tmp$Cb0.5.vs.Cv = tmp$Cb0.5 / tmp$Cv   ###2013 Dec 10  change  !!!!!!!!
  tbH2O2Report$Cv.vs.Cb[i] = mean(tmp$Cv.vs.Cb, na.rm=T)
  tbH2O2Report$Cv.vs.Cb.STD[i] = sd(tmp$Cv.vs.Cb, na.rm=T)
  tbH2O2Report$Cb.vs.Cv[i] =  1/ tbH2O2Report$Cv.vs.Cb[i]
  tbH2O2Report$Cb.vs.Cv.STD[i] = sd(tmp$Cb.vs.Cv, na.rm=T)
  
  tbH2O2Report$Cb0.5.vs.Cv[i] = mean(tmp$Cb0.5.vs.Cv)
  tbH2O2Report$Cb0.5.vs.Cv.STD[i] = sd(tmp$Cb0.5.vs.Cv, na.rm=T)
  #tbH2O2$Cb0.5.vs.Cv = tbH2O2$Cb0.5 / tbH2O2$Cv

  tbH2O2Report$Cb0.5Mean[i] = mean(tmp$Cb0.5, na.rm=T)
  tbH2O2Report$Cb.05STD[i] = sd(tmp$Cb0.5, na.rm=T)
  tbH2O2Report$b0.5minH2O2[i] = mean(tmp$b0.5.min, na.rm=T)
  tbH2O2Report$b0.5maxH2O2[i] = mean(tmp$b0.5.max, na.rm=T)
  
  #tbH2O2Report$Cb0.5.vs.Cv[i] = mean(tmp$Cb0.5.vs.Cv, na.rm=T) 
  #tbH2O2Report$Cb0.5.vs.Cv.STD[i] = sd(tmp$Cb0.5.vs.Cv, na.rm=T)
  #tbH2O2Report$half.vs.full[i] = sd(tmp$half.vs.full, na.rm=T)
    
  #tbH2O2Report$Cv.vs.Cb.byMeans[i] = tbH2O2Report$CvMean[i] / tbH2O2Report$CbMean[i]
  #tbH2O2Report$Cv.vs.Cb.byMeansSTD[i] = tbH2O2Report$Cv.vs.Cb.byMeans[i]*sqrt( (tbH2O2Report$CvSTD[i] / tbH2O2Report$CvMean[i])^2 + (tbH2O2Report$CbSTD[i] / tbH2O2Report$CbMean[i])^2 )

  #tbH2O2Report$Cb.vs.Cv.byMeans[i] = tbH2O2Report$CbMean[i] / tbH2O2Report$CvMean[i]    
  #tbH2O2Report$Cb.vs.Cv.byMeansSTD[i] = tbH2O2Report$Cb.vs.Cv.byMeans[i]*sqrt( (tbH2O2Report$CvSTD[i] / tbH2O2Report$CvMean[i])^2 + (tbH2O2Report$CbSTD[i] / tbH2O2Report$CbMean[i])^2 ) 
}
rownames(tbH2O2Report)= tbH2O2Report$strains

#write.csv(tbH2O2Report, file="output/LOHH2O2_averaged20131209_v2.csv", row.names=F)
write.csv(tbH2O2Report, file="output/LOHH2O2_averaged20131210_v1.csv", row.names=F)

#WildStrains=c('101S','M1-2', 'M13','M2-8','M32','M34','M5','M8','SGU57','YPS128','YPS163')
#tbH2O2ReportNat = tbH2O2Report[NatStrains,]

#write.csv(tbH2O2ReportNat, file="output/LOHH2O2_all_20131126.csv", quote=F, row.names=F)

#quit("no")

