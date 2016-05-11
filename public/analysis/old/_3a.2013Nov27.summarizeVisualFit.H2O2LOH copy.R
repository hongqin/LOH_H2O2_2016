# summarize H2O2-LOH data
#2013 Nov 20, calculate standard deviations

rm(list=ls())

setwd("~/github/LOH_H2O2_2012-master/analysis")
list.files(, pattern="csv")

#read individual fitting data, organized by experiments and strains
#2013 manual fitting results
tbH2O2 = read.csv("H2O2_Log_Plot_Summarized_data,2013May30.csv")
names(tbH2O2)=c("file", "Date","Strain","Cv","Cb", "OD", "note13", "repeat")
tbH2O2$file = as.character(tbManu$file)

#tbH2O2 = read.csv("_ctl.tb_out.20130530b.csv")
#tbH2O2$files = as.character( tbH2O2$files )
#str(tbH2O2[,1:4])

#find strains, gsub
tbH2O2$strains = gsub("[\\,|\\.].*", "", tbH2O2$file)

strains2 = unique(tbH2O2$strains)[c(1:2,3:8,10:15)]
table(tbH2O2$strains)

tbH2O2$Cb.vs.Cv = tbH2O2$Cb / tbH2O2$Cv
tbH2O2$Cv.vs.Cb = tbH2O2$Cv / tbH2O2$Cb

tbH2O2Report = data.frame(strains2)
names( tbH2O2Report ) = c("strains")
# strain = 'BY4743' i =1
for( i in 1:length(strains2) ) {  
  tmp = tbH2O2[tbH2O2$strains == strains2[i],  ]  
  tbH2O2Report$CvMean[i] = mean(tmp$Cv, na.rm=T)
  tbH2O2Report$CvSTD[i] = sd(tmp$Cv, na.rm=T)  
  tbH2O2Report$CbMean[i] = mean(tmp$Cb, na.rm=T)
  tbH2O2Report$CbSTD[i] = sd(tmp$Cb, na.rm=T)
#  tbH2O2Report$bminH2O2[i] = mean(tmp$bmin, na.rm=T)
#  tbH2O2Report$bmaxH2O2[i] = mean(tmp$bmax, na.rm=T)
  
  #Two ways to calculate the ratios
  tbH2O2Report$Cv.vs.Cb[i] = mean(tmp$Cv.vs.Cb, na.rm=T)
  tbH2O2Report$Cv.vs.Cb.STD[i] = sd(tmp$Cv.vs.Cb, na.rm=T)
  tbH2O2Report$Cb.vs.Cv[i] = mean(tmp$Cb.vs.Cv, na.rm=T) 
  tbH2O2Report$Cb.vs.Cv.STD[i] = sd(tmp$Cb.vs.Cv, na.rm=T)
  
  tbH2O2Report$Cv.vs.Cb.byMeans[i] = tbH2O2Report$CvMean[i] / tbH2O2Report$CbMean[i]
  tbH2O2Report$Cv.vs.Cb.byMeansSTD[i] = tbH2O2Report$Cv.vs.Cb.byMeans[i]*sqrt( (tbH2O2Report$CvSTD[i] / tbH2O2Report$CvMean[i])^2 + (tbH2O2Report$CbSTD[i] / tbH2O2Report$CbMean[i])^2 )

  tbH2O2Report$Cb.vs.Cv.byMeans[i] = tbH2O2Report$CbMean[i] / tbH2O2Report$CvMean[i]    
  tbH2O2Report$Cb.vs.Cv.byMeansSTD[i] = tbH2O2Report$Cb.vs.Cv.byMeans[i]*sqrt( (tbH2O2Report$CvSTD[i] / tbH2O2Report$CvMean[i])^2 + (tbH2O2Report$CbSTD[i] / tbH2O2Report$CbMean[i])^2 ) 
}
rownames(tbH2O2Report)= tbH2O2Report$strains

write.csv(tbH2O2Report, file="output/LOHH2O2_all_20131127.csv", quote=F, row.names=F)

#WildStrains=c('101S','M1-2', 'M13','M2-8','M32','M34','M5','M8','SGU57','YPS128','YPS163')
#tbH2O2ReportNat = tbH2O2Report[NatStrains,]

#write.csv(tbH2O2ReportNat, file="output/LOHH2O2_all_20131126.csv", quote=F, row.names=F)

#quit("no")

