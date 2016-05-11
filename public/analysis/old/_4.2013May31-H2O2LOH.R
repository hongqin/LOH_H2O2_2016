# use 2013 manual fitting results, and gnls reults

#reverse L0 axis in plot, 2012 July 30

# 2012Feb25, Tg.vs.Tc ~ ln(R0) + G
# Partial correlations are all negative, what does this mean? 
# Tg/Tc is a measure of ability to maintian recombiation rate during aging

 rm( list = ls() );
 #setwd("~/projects/LOH-oxidants2012.osX/analysis")
 setwd("~/github/LOH_H2O2_2012-master/analysis")
 
 #Load previous results
 tb = read.table("summary.new.by.strain.csv", header=T, sep="\t");
 tb.old = tb;
 labels = names( tb.old );
 tb = tb.old[, c("strain","ARLS","R0","G","CLS","Tc", "Tg","Tmmax","Tbmax", "Td", "Tdmax","TLmax","Lmax", 
 "b.max", "b.min", "strains", "L0.all", "L0.small" , "Pbt0","Pb0.5t0", "Pbt0.b") ];
 tb$CLS.vs.Tc = tb$CLS / tb$Tc; 
 tb$Tg.vs.Tc = tb$Tg / tb$Tc;
 tb$strain = as.character(tb$strain)

 #load exg06 data
 nat = read.table("062705.rls.cls.tab", sep="\t", header=T, colClasses=c("character", rep(NA,4)) );
 
 
###############
 #Load H2O2-LOH results
 list.files(pattern="csv")
 #tb2 = read.csv("H2O2_Log_Plot_Summarized_data,2012Jan24.csv")
 #names(tb2) = c("Date", "Strain", "Cv", "Cb", "OD600nm","Notes","repeat")
 tb2 = read.csv("H2O2_Log_Plot_Summarized_data,2013May30.csv") 
 names(tb2) = c("files","Date", "Strain", "Cv", "Cb", "OD600nm","Notes","repeat")

 tb2$Strain = as.character(tb2$Strain)
 tb2.old = tb2;
 tb2 = tb2[, c("Strain","Date","Cv","Cb")]
 tb2$Cv.vs.Cb = tb2$Cv / tb2$Cb
 tb2$Cb.vs.Cv = tb2$Cb / tb2$Cv
 head(tb2)
 
 tb3 = read.csv("_ctl.tb_out.20130530b.csv")
 tb3$strains = as.character( tb3$strains )
 strains3 = unique(tb3$strains)
 strains3 = strains3[ ! is.na(strains3)]
 
 #check strains names, do they match? 
 strains2 = unique(tb2$Strain)
 intersect( strains2, tb$strain)
 intersect( strains3, tb$strain)
 
 intersect( nat$strain, strains2)
 intersect( nat$strain, strains3)
###########
 
 # generate the means tb3m
 tb3m = data.frame(cbind(strains3 ))
 tb3m[,1] = as.character(tb3m[,1])          
 i=1;
 for( i in 1:length(strains3) ){
   sub = tb3[tb3$strains==strains3[i], ]
   tb3m$Cv.vs.Cb[i] = mean(sub$Cv.vs.Cb, na.rm=T)
   tb3m$Cv[i] = mean(sub$Cv, na.rm=T)
   tb3m$Cb[i] = mean(sub$Cb, na.rm=T)
 }
 tb3m$Cb.vs.Cv = 1 / tb3m$Cv.vs.Cb
 
 
#### analyze H2O2-LOH in tb2
 #raw values
 hist(tb2$Cv.vs.Cb, br=10)
 hist(log2(tb2$Cv.vs.Cb), br=10)
 summary(tb2) 
 hist(log2(1/tb2$Cv.vs.Cb), br=10)

 # generate the means
 tb2m = data.frame(cbind(strains2 ))
 tb2m[,1] = as.character(tb2m[,1])          
 i=1;
 for( i in 1:length(strains2) ){
   sub = tb2[tb2$Strain==strains2[i], ]
   tb2m$Cv.vs.Cb[i] = mean(sub$Cv.vs.Cb, na.rm=T)
   tb2m$Cb.vs.Cv[i] = mean(sub$Cb.vs.Cv, na.rm=T)
   tb2m$Cv[i] = mean(sub$Cv, na.rm=T)
   tb2m$Cb[i] = mean(sub$Cb, na.rm=T)
 }
 tb2m$Cv.vs.CbByMean = tb2m$Cv / tb2m$Cb ;   
 hist( 1/ tb2m$Cv.vs.Cb, br =10) 
 hist( 1/ tb2m$Cv.vs.CbByMean, br =10) 
 plot( tb2m$Cv.vs.Cb ~ tb2m$Cv.vs.CbByMean )

 # compare Cb/Cv and Tg/Tc
 t.test( 1 / tb2m$Cv.vs.Cb, mu=1, alternative="less") #p=0.074
 t.test( 1 / tb2m$Cv.vs.CbByMean, mu=1, alternative="less") #p=0.32
 
 t.test( log2(1 / tb2m$Cv.vs.Cb), mu=0, alternative="less") #p=0.020
 t.test( log2(1 / tb2m$Cv.vs.CbByMean), mu=0, alternative="less") #p=0.093
 
 wilcox.test( 1/ tb2m$Cv.vs.Cb, mu=1, alternative="less") #p=0.0863
 #wilcox.test( 1/ tb2m$Cv.vs.CbByMean, mu=1, alternative="less") #p=0.38
 # Cb/Cv < 1
 
 hist(tb$Tg.vs.Tc, br=10)  
 t.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00035
 wilcox.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00042
 #Tg/Tc > 1

 ### side by side bar-plots of Tg/Tc Cb/Cv
 mystep=0.2
 my.breaks = seq( 0.1,  round(max( c( 1/tb2m$Cv.vs.Cb[1:14], tb$Tg.vs.Tc ) + 0.1, 1)), by= mystep ); 
 h.H2O2  <- hist( 1/tb2m$Cv.vs.Cb, br= my.breaks, xlab = "Cb/Cv", ylab = "relative density", freq=F ) ; 
 h.aging <- hist(tb$Tg.vs.Tc, br= my.breaks, xlab = "Tg/Tc",  ylab = "relative density", freq=F ) ;
 
 #generate the comparison table
 bins <-  data.frame( rbind(h.H2O2$density,h.aging$density) )  ;
 my.mids = my.breaks[-length(my.breaks)] + mystep/2
 #my.mids
 names( bins ) <- my.mids
 row.names(bins) <- c( "H2O2", "Chronological Aging" )
 bins

 pdf("Figure_sideBYside.pdf", width=8, height=5)
 barplot( as.matrix(bins), beside=T, col=c("black","gray"), ylab="Relative Frequency", xlab="Ratios",
          legend= c( "Cb/Cv H2O2", "Tg/Tc Aging" )  );
 title(main="H2O2 and chronological aging elevate LOH at different modes" )
 dev.off(); 
  
#################
 ### merge tb tb2m tb3m 
 tb$Cv.vs.Cb   =  tb2m$Cv.vs.Cb[match(tb$strain, tb2m$strains2)]
 tb$Cv.vs.CbGnls = tb3m$Cv.vs.Cb[match(tb$strain, tb3m$strains3)]
 tb$CvGnls = tb3m$Cv[match(tb$strain, tb3m$strains3)]
 #tb$Cb.vs.Cv = tb2m$Cb.vs.Cv[match(tb$strain, tb2m$strains2)]
 tb$Cb.vs.Cv = 1 / tb$Cv.vs.Cb; 
 tb$Cb.vs.CvGnls = 1 / tb$Cv.vs.CbGnls; 
 #This is for consistency, because 
 # 1/ratio and ratio before and after averaging, give slightly different values. 
 
 tb$Cb = tb2m$Cb[match(tb$strain, tb2m$strains2)]
 tb$Cv = tb2m$Cv[match(tb$strain, tb2m$strains2)] 
 tb$Cb.vs.Cv.mean = tb$Cb / tb$Cv

 plot( tb$Cb.vs.Cv ~ tb$Cb.vs.Cv.mean, xlim=c(0,2),ylim=c(0,2))

###### merge nat, tb2m, tb3m
 nat$Cv.vs.Cb = tb2m$Cv.vs.Cb[match(nat$strain, tb2m$strains2)]
 nat$Cv.vs.CbGnls = tb3m$Cv.vs.Cb[match(nat$strain, tb3m$strains3)]
 nat$CvGnls = tb3m$Cv[match(nat$strain, tb3m$strains3)]
 nat$CbGnls = tb3m$Cb[match(nat$strain, tb3m$strains3)]
 
 #tb$Cb.vs.Cv = tb2m$Cb.vs.Cv[match(tb$strain, tb2m$strains2)]
 nat$Cb.vs.Cv = 1 / nat$Cv.vs.Cb;  
 nat$Cb.vs.CvGnls = 1 / nat$Cv.vs.CbGnls;  
 nat$Cv = tb2m$Cv[match(nat$strain, tb2m$strains2)]
 nat$Cb = tb2m$Cb[match(nat$strain, tb2m$strains2)]
 
  
### regression analysis 
pTb = 1: length(tb[1,])
names(pTb) = names(tb)
 for( j in c(2:15,17:34) ) {
   m = lm( tb[, j] ~ tb$Cv.vs.Cb)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 #ARLS          Lmax        L0.all     CLS.vs.Tc      Cv.vs.Cb      Cb.vs.Cv 
 #3.546097e-02  3.877701e-02  3.137392e-02  4.748623e-02  0.000000e+00  7.306331e-05 
 #Cb.vs.Cv.mean 
 #5.729385e-04
 pTb[pTb>0.05 & pTb<0.15]

 tiff("ARLS-CvCb-20130618.tif",width=480,heigh=480)
 par(font=2) 
 plot( tb$ARLS ~ tb$Cv.vs.Cb , pch=19, col="red", main="H2O2-LOH ~ ARLS, 20130618", ylim=c(22,38),xlim=c(0.5,4.2)
       , ylab='ARLS',xlab='Cv/Cb Tolerance to H2O2-induced genomic instability')
 text( tb$Cv.vs.Cb+0.08, tb$ARLS+0.5, tb$strain)
 m = lm(tb$ARLS ~ tb$Cv.vs.Cb  )
 abline( m, col="blue")
 summary(m)
 text(3.5, 35,  "R2=0.37 p=0.035")
 dev.off()
 
 plot( tb$Cv.vs.Cb ~ tb$L0.all , pch=19, col="red", main="H2O2-LOH ~ mitotic asymmetry, 20130531" )
 text( tb$L0.all, tb$Cv.vs.Cb, tb$strain)
 m = lm(tb$Cv.vs.Cb ~ tb$L0.all )
 abline( m, col="blue")
 summary(m)
 text(0.2, 3.5, "R2=0.42 p=0.031")
 
 
 
 #now, how does Cb/Cv? 
 pTb2 = 1: length(tb[1,])
 names(pTb2) = names(tb)
 for( j in c(2:15,17:34) ) {
   m = lm( tb[, j] ~ tb$Cb.vs.Cv)
   sm = summary(m)
   pTb2[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb2[pTb2<0.05]
 
 ### Cv/Cb or Cb/Cv ~ robustness? I need a positive proxy
 summary(lm( tb$Cv.vs.Cb ~ tb$G ) )  #positive, p=0.37, 
 summary(lm( tb$Cb.vs.Cv ~ tb$G ) )  #negative, p=0.13
 
 ### 
 summary(lm( tb$Tg.vs.Tc ~ tb$G ) )  #negative, p=0.82
 summary(lm( tb$Tg.vs.Tc ~ tb$ARLS ) )  #negative, p=0.82
 summary(lm( 1/ tb$Tg.vs.Tc ~ tb$G ) )  #negative, p=0.82
 summary(lm( 1/ tb$Tg.vs.Tc ~ tb$ARLS ) )  #negative, p=0.82
 summary(lm( 1/ tb$Tg.vs.Tc ~ tb$ARLS ) )  #negative, p=0.82
 
 
 
 
 pTb = 1: length(tb[1,])
 names(pTb) = names(tb)
 for( j in c(2:15,17:34) ) {
   m = lm( tb[, j] ~ tb$Cv.vs.CbGnls)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 #pTb[pTb<0.05]
 #b.max        b.min     Tg.vs.Tc Cv.vs.CbGnls Cb.vs.CvGnls           Cb      frac.gc 
 #1.347212e-02 6.211504e-03 2.000139e-02 0.000000e+00 4.460174e-09 4.059540e-02 2.000139e-02 
 
 pTb[pTb>0.05 & pTb<0.15]
 
 summary(lm(tb$CLS ~ tb$Cv.vs.Cb)) #p=0.21
 
 pTb = 1: length(tb[1,])
 names(pTb) = names(tb)
 for( j in c(2:15,17:34) ) {
   m = lm( tb[, j] ~ tb$Cv)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 pTb[pTb>0.05 & pTb<0.15]

 
 pTb = 1: length(tb[1,])
 names(pTb) = names(tb)
 for( j in c(2:15,17:34) ) {
   m = lm( tb[, j] ~ tb$Cb)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 pTb[pTb>0.05 & pTb<0.15]
 
 pTb = 1: length(tb[1,])
 names(pTb) = names(tb)
 for( j in c(2:15,17:31) ) {
   m = lm( tb[, j] ~ tb$Cv.vs.CbGnls)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 
 
 summary( lm( nat$Cb.vs.Cv ~ nat$rls + nat$a + nat$b + nat$cls ) )
 summary( lm( nat$Cv ~ nat$rls + nat$a + nat$b + nat$cls ) )
 summary( lm( nat$Cb.vs.Cv ~ nat$rls + nat$cls ) ) #p=0.06
 summary( lm( nat$Cv ~ nat$rls + nat$cls ) ) 
 summary( lm( nat$Cv ~ nat$cls ) ) 
 summary( lm( nat$CvGnls ~ nat$cls ) ) 
 
 summary( lm( log(nat$Cv) ~ nat$cls ) ) 
 summary( lm( 1/ log(nat$Cv) ~ nat$cls ) ) 
 summary( lm( 1/ log(nat$Cb) ~ nat$cls ) ) 
 summary( lm( log(nat$Cb) ~ log(nat$cls )) ) 
 summary( lm( nat$Cb ~ nat$cls ) ) 
 plot( nat$Cb, nat$cls)
 points( nat$CbGnls, nat$cls, pch=19, col='red')
 
 plot( nat$Cv, nat$cls)
 plot( nat$CvGnls, nat$cls)
 points( nat$Cv, nat$cls, pch=19)
 
 summary( lm( nat$Cv.vs.CbGnls ~ nat$rls + nat$a + nat$b + nat$cls ) ) #p 0.71
 summary( lm( nat$Cv.vs.CbGnls ~ nat$rls + nat$cls ) ) #p=0.86
 summary( lm( nat$Cb.vs.CvGnls ~ nat$rls + nat$cls ) ) #p=0.81

 summary( lm( nat$Cb.vs.Cv ~ nat$rls ) ) #p=0.066 , negative correlation
 
 summary( lm( nat$Cb.vs.Cv ~ nat$cls ) ) #p=0.06 , negative correlation
 summary( lm( nat$Cv.vs.Cb ~ nat$cls ) ) #p=0.13, positive correlation

 plot( 1/ nat$Cv.vs.Cb ~ nat$cls, pch=19, col="red" )
 text( nat$cls, 1/nat$Cv.vs.Cb, nat$strain) 
 text( 10, 1.5, "p=0.06")
 
 summary( lm( nat$Cb.vs.CvGnls ~ nat$cls ) ) #p=0.6 , negative correlation
 summary( lm( nat$Cv.vs.CbGnls ~ nat$cls ) ) #p=0.6, positive correlation
 
 summary( lm( nat$Cb.vs.CvGnls ~ nat$rls ) ) #p=0.56
 summary( lm( 1 / nat$Cb.vs.CvGnls ~ nat$rls ) ) #p=0.68
 summary( lm( nat$Cv.vs.Cb ~ nat$rls ) ) #p=0.057 negative
 plot( nat$Cv.vs.Cb ~ nat$rls )
 points( nat$rls, nat$Cv.vs.CbGnls, col="red", pch=19 )
 
 summary( lm( nat$Cb.vs.Cv ~ nat$b ) ) #p=0.27
 summary( lm( nat$Cb.vs.Cv ~ nat$a ) ) #p=0.66
 summary( lm( nat$Cv.vs.Cb ~ nat$a ) ) #p=0.66
 
 summary( lm( tb$Cb.vs.Cv ~ tb$ARLS + tb$Tg + tb$Tc + tb$Tg.vs.Tc + tb$CLS ) )
 summary( lm( tb$CvGnls ~ tb$CLS + tb$ARLS + tb$Tg + tb$Tg.vs.Tc + tb$R0 + tb$G) ) #p=0.85

 summary( lm( tb$Cb.vs.Cv ~ tb$ARLS ) )  #p = 0.134
 summary( lm( tb$Cb.vs.Cv.mean ~ tb$ARLS ) )  #p = 0.4
 summary( lm( tb$Cb.vs.Cv ~ tb$R0 + tb$G ) ) # p =0.18
 summary( lm( tb$Cb.vs.Cv ~ tb$R0 * tb$G ) ) # p =0.38
 summary( lm( tb$Cb.vs.Cv ~ sqrt(tb$R0 * tb$G) ) ) # p =0.77
 
 summary( lm( tb$Cb.vs.Cv ~ tb$Tg.vs.Tc) ) #p = 0.4
 summary( lm( tb$Cb ~ tb$Tg.vs.Tc + tb$ARLS + tb$CLS ) ) #p=0.30
 
 summary( lm( tb$Cb ~ tb$CLS ) ) #p=0.23
 summary( lm( tb$Cv ~ tb$CLS ) ) #p=0.85
 summary( lm( tb$CvGnls ~ tb$CLS ) ) #p=0.85
 
 summary( lm( tb2m$Cb ~ tb2m$Cv ) ) #p0.09
 m = lm( tb2m$Cb ~ tb2m$Cv )
 plot( tb2m$Cb ~ tb2m$Cv )
 abline( m , col='red')
 
 summary( lm( tb$Cb.vs.Cv ~ tb$Tc ) ) #p=0.15
 summary( lm( tb$Cb.vs.Cv.mean ~ tb$Tc ) ) #p=0.15
 
 summary( lm( tb$CLS ~ tb$Cb.vs.Cv ) ) #p=0.1
 plot( tb$CLS ~ tb$Cb.vs.Cv, col="red" )
 text( jitter(tb$Cb.vs.Cv), jitter(tb$CLS), tb$strain)
 text( jitter(tb$Cb.vs.Cv), jitter(tb$CLS), seq(1:11))
 
 #summary( lm( tb$CLS ~ tb$Cb.vs.Cv ) ) #using averaged values, p=0.0705
 
 summary( lm( tb$CLS ~ tb$Cv.vs.Cb ) ) #p=0.2
 plot( tb$CLS ~ tb$Cv.vs.Cb )
 text( tb$Cv.vs.Cb, tb$CLS, tb$strain)
 
 
 pdf("Figure_CLS-CbCv.pdf", width=7, height=7)
 summary( lm( tb$CLS ~ tb$Cb.vs.Cv ) ) #p=0.024  !!!!! negative !!!why
 m = lm( tb$CLS ~ tb$Cb.vs.Cv )
 plot( tb$CLS ~ tb$Cb.vs.Cv, pch=19, xlim=c(0, 1.6),xlab='Cb/Cv',ylab='CLS' )
 abline( m, col='red')
 text( (tb$Cb.vs.Cv + 0.03), (tb$CLS + 0.3), tb$strain)
 text(1.0, 14, "R^2=0.54, p=0.024")
 dev.off()
   
 plot( tb$CLS ~ tb$Cb )
 plot( tb$CLS ~ tb$Cv )
 
 summary( lm( tb$L0.all ~ tb$Cb.vs.Cv ) ) #p=0.054 !!!!! positive 
 #this suggest H2O2 effect ~ asymetry --> asymmtry is caused by H2O2 distribution? 
 m = lm( tb$L0.all ~ tb$Cb.vs.Cv )
 pdf("Figure_L0-CbCv.pdf", width=7, height=7)
 plot( tb$L0.all ~ tb$Cb.vs.Cv, pch=19, xlim=c(0,1.6),ylim=c(0.28,0.01), xlab='Cb/Cv',ylab='L0' )  
 abline( m, col='red')
 text( (tb$Cb.vs.Cv + 0.02), (tb$L0.all + 0.008), tb$strain)
 text(0.2, 0.25, 'R^2=0.43, p=0.055')
 dev.off()
 
 
  
#quit("yes")

####
### END
####
 
