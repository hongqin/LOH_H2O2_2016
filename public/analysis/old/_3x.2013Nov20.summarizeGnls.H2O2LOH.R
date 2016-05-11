# bactch automatical fitting is not as good as manual fitting results. 
# summarize H2O2-LOH data
#2013 Nov 20, calculate standard deviations

rm(list=ls())

setwd("~/github/LOH_H2O2_2012-master/analysis")
list.files(, pattern="csv")

#read individual fitting data, organized by experiments and strains
tbH2O2 = read.csv("_ctl.tb_out.20130530b.csv")
tbH2O2$files = as.character( tbH2O2$files )
str(tbH2O2[,1:4])

#find strains, gsub
tbH2O2$strains = gsub("[\\,|\\.].*", "", tbH2O2$files)

strains2 = unique(tbH2O2$strains)
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
  tbH2O2Report$bmin[i] = sd(tmp$Cb, na.rm=T)
  
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

write.csv(tbH2O2Report, file="output/LOHH2O2_all_20131125.csv", quote=F, row.names=F)

NatStrains=c('101S','M1-2', 'M13','M2-8','M32','M34','M5','M8','SGU57','YPS128','YPS163')
tbH2O2ReportNat = tbH2O2Report[NatStrains,]

write.csv(tbH2O2ReportNat, file="output/LOHH2O2_nat_20131125.csv", quote=F, row.names=F)

quit("no")

########### test on H2O2-LOH 
t.test( tbH2O2$Cv.vs.Cb, mu=1, alternative="greater") 
summary(tbH2O2$Cv.vs.Cb)
hist(tbH2O2$Cv.vs.Cb, br=20)

t.test( tbH2O2Report$Cv.vs.Cb, mu=1, alternative="greater") 
t.test( tbH2O2Report$Cv.vs.Cb.byMeans, mu=1, alternative="greater") 

summary( tbH2O2Report )
hist(tbH2O2Report$Cv.vs.Cb, br=20)
hist(tbH2O2Report$Cb.vs.Cv, br=10)

####### read previous data
#Load previous results
tb = read.table("summary.new.by.strain.csv", header=T, sep="\t");
tb.old = tb;
labels = names( tb.old );
tb = tb.old[c(1:11), c("strains","ARLS","R0","G","CLS","Tc", "Tg","Tmmax","Tbmax", "Td", "Tdmax","TLmax","Lmax", 
                       "b.max", "b.min", "strains", "L0.all", "L0.small" , "Pbt0","Pb0.5t0", "Pbt0.b") ];
tb$CLS.vs.Tc = tb$CLS / tb$Tc; 
tb$Tg.vs.Tc = tb$Tg / tb$Tc;
tb$strains = as.character(tb$strains)

#load qin06exg data
nat = read.table("062705.rls.cls.tab", sep="\t", header=T, colClasses=c("character", rep(NA,4)) );
# Qin06exg data is not completely included in Qin08 data
names(nat) = c("strains", "ARLS", "R0", "G", "CLS")

#check strains names, do they match? 
strains2 = unique(tbH2O2Report$strains)
intersect( strains2, tb$strains)

# merge tb and tbH2O2Report
tb1308 = merge(tb, tbH2O2Report)
tb1306 = merge(nat, tbH2O2Report)

############### regression Cv/Cb, R, G, RLS, CLS
summary( lm (Cv.vs.Cb ~ log(R0) + G, data=tb1308 ) ) #p=0.5
summary( lm (Cb.vs.Cv ~ log(R0) + G, data=tb1308 ) ) #p=0.7

summary( lm (Cb.vs.Cv ~ log(R0) + G + ARLS + CLS, data=tb1306 ) ) #p=0.015
summary( lm (Cb.vs.Cv ~ log(R0) + G         + CLS, data=tb1306 ) ) #p=0.05
summary( lm (Cb.vs.Cv ~ log(R0) + G + ARLS, data=tb1306 ) ) #p=0.008, but factors are not independent
summary( lm (Cb.vs.Cv ~  CLS + ARLS, data=tb1306 ) ) #p=0.008

summary( lm (Cv.vs.Cb ~ CLS, data=tb1308 ) ) #p=0.7
summary( lm (Cb.vs.Cv ~ ARLS, data=tb1308 ) ) #p=0.7

pv = 1:length(tb1308[1, ])
names(pv)=names(tb1308)[ 1:length(tb1308[1, ]) ]
for( jj in 1:length(tb1308[1, ])) {
  m = lm (tb1308$Cv.vs.Cb ~ tb1308[,jj] ) 
  sm = summary(m)
  pv[jj] = 1- pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
}

pv2 = 1:length(tb1306[1, ])
names(pv2)=names(tb1306)[ 1:length(tb1306[1, ]) ]
for( jj in 1:length(tb1306[1, ])) {
  m = lm (tb1306$Cv.vs.Cb ~ tb1306[,jj] ) 
  sm = summary(m)
  pv2[jj] = 1- pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
}



