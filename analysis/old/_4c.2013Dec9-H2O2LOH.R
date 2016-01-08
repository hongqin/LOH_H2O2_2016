# 2013 Dec 9, using results from _2c2.20131209_update_fitting

# 2013Nov25, use summarized results from _3.*R, which contains standard deviations.  
# 2008 CLS data uses M14, 2013H2O2 uses SGU57, so only 10 strains overlaps

# Use 2013 fitting results and gnls reults
# 2013Nov14, use correlated regression
# Reverse L0 axis in plot, 2012 July 30

# 2012Feb25, Tg.vs.Tc ~ ln(R0) + G
# Partial correlations are all negative, what does this mean? Of couse, they should be negative correlations. 
# Tg/Tc is a measure of ability to maintian recombiation rate during aging

 rm( list = ls() );
 #setwd("~/projects/LOH-oxidants2012.osX/analysis")
 setwd("~/github/LOH_H2O2_2012-master/analysis")

###############
 #Load 2013 H2O2-LOH results
 list.files(pattern="csv")
 #tb2 = read.csv("output/LOHH2O2_averaged20131206.csv")
 #rownames(tb2) = as.character(tb2$strains)
 #tb2= tb2[-c(3,14), ]
 tb2 = read.csv("_merged.tb.20131209_9,15am.csv") 
 tb2$Cv.vs.Cb = tb2$Cv / tb2$Cb  
 
 #### analyze H2O2-LOH in tb2
 hist(tb2$Cv.vs.Cb, br=10)
 hist(log2(tb2$Cv.vs.Cb), br=10)
 summary(tb2) 
 hist(log2(1/tb2$Cv.vs.Cb), br=10)
 
##### #Load previous LOH-CLS results, Qin Plos One 2008
 tb = read.table("summary.new.by.strain.csv", header=T, sep="\t");
 tb.old = tb;
 labels = names( tb.old );
 tb = tb.old[, c("strain","ARLS","R0","G","CLS","Tc", "Tg","Tmmax","Tbmax", "Td", "Tdmax","TLmax","Lmax", 
 "b.max", "b.min", "strains", "L0.all", "L0.small" , "Pbt0","Pb0.5t0", "Pbt0.b") ];
 tb$CLS.vs.Tc = tb$CLS / tb$Tc; 
 tb$Tg.vs.Tc = tb$Tg / tb$Tc;
 tb$strain = as.character(tb$strain)
 tb.old = tb;
 tb = tb.old[1:13,] #remove rad52DD
 tb.test = tb[1:11,]   
 
 summary(lm(Tg.vs.Tc ~ ARLS, data=tb.test)) #re-run the old results, just to double-check
 plot( tb.test$Tg.vs.Tc ~ tb.test$ARLS )
 text(tb.test$ARLS, tb.test$Tg.vs.Tc, tb.test$strains)
 
 summary(lm( 1/Tg.vs.Tc ~ ARLS, data=tb.test)) #good, negative, p=0.012 
 
 ###############
 #load exg06 data
 nat = read.table("062705.rls.cls.tab", sep="\t", header=T, colClasses=c("character", rep(NA,4)) );
 
###### merge tb and tb2
 tb2$strain = tb2$Strain
 tb3a = merge(tb, tb2, by='strain', all=T)

summary(lm(tb3a$Cv.vs.Cb ~ tb3a$ARLS))
 summary(lm(1/tb3a$Cv.vs.Cb ~ tb3a$ARLS))

 summary(lm(tb3a$Cv.vs.Cb.Manu ~ tb3a$ARLS)) #p=0.002
 summary(lm(tb3a$Cv.vs.Cb.gnls ~ tb3a$ARLS)) #p=0.6997

plot( tb3a$Cv.vs.Cb.Manu ~ tb3a$ARLS, col='white' )
 text( tb3a$ARLS, tb3a$Cv.vs.Cb.Manu, tb3a$strain, offset=-5 )
# text( tb3a$ARLS, tb3a$Cv.vs.Cb.Manu, tb3a$Date, offset=10 )
tb3a[grep("M32",tb3a$file), ] 


 plot( tb3a$Cv.vs.Cb.gnls ~ tb3a$ARLS, col='white' )
 text( tb3a$ARLS, tb3a$Cv.vs.Cb.gnls, tb3a$strain )
 # text( tb3a$ARLS, tb3a$Cv.vs.Cb.Manu, tb3a$Date, offset=10 )
 tb3a[grep("M32",tb3a$file), ] 
tb3a$Cv.vs.Cb.gnls[ tb3a$Cv.vs.Cb.gnls > 5 ] = NA
  summary(lm(tb3a$Cv.vs.Cb.gnls ~ tb3a$ARLS)) #p=0.39

 summary(lm(tb3a$Cv.vs.Cb.gnls ~ tb3a$Cv.vs.Cb.Manu)) #p=0.53

 plot(tb3a$Cv.vs.Cb.gnls ~ tb3a$Cv.vs.Cb.Manu)
 
 plot(tb3a$Cv ~ tb3a$CvManu); text( tb3a$CvManu, tb3a$Cv, tb3a$strain)
 tb3a = tb3a[! is.na(tb3a$file), ] 
 tb3a[tb3a$Cv>0.4, c('file','Cv')] 
 tb3a$Cv[tb3a$Cv>0.4] = tb3a$CvManu
 
  plot(tb3a$Cb ~ tb3a$CbManu)
 
 summary(lm( tb3a$CvManu / tb3a$Cb ~ tb3a$ARLS))

 summary(lm( tb3a$Cv / tb3a$CbManu ~ tb3a$ARLS))
 
# summary(lm(1/tb3a$Cb.vs.Cv.byMeans ~ tb3a$ARLS))
#summary(lm(tb3a$Cv.vs.Cb.byMeans ~ tb3a$ARLS))
 
 #nat$strains = nat$strain
 #tb3b = merge(tb2, nat, by='strains', all=T)
 
 # compare Cb/Cv and Tg/Tc 
 hist(tb$Tg.vs.Tc, br=10)  
 t.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00035
 wilcox.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00042
 #Tg/Tc > 1

 ### side by side bar-plots of Tg/Tc Cb/Cv
 mystep=0.2
 my.breaks = seq( 0.1,  round(max( c( 1/tb2$Cv.vs.Cb, tb$Tg.vs.Tc ) + 0.1, 1)) ,by= mystep ); 
 my.breaks3 = seq( 0.1,  round(max( c( 1/tb2$Cv.vs.Cb, tb$Tg.vs.Tc, 1/tb2$Cv.vs.Cb.byMeans ) + 0.1, 1)) ,by= mystep ); 
 h.H2O2  <- hist( 1/tb2$Cv.vs.Cb, br= my.breaks, xlab = "Cb/Cv", ylab = "relative density", freq=F ) ; 
 h.H2O2Mean  <- hist( 1/tb2$Cv.vs.Cb.byMeans, br= my.breaks, xlab = "Cv/CbByMeans", ylab = "relative density", freq=F ) ; 
 
 h.aging <- hist(tb$Tg.vs.Tc, br= my.breaks, xlab = "Tg/Tc",  ylab = "relative density", freq=F ) ;
 
 #generate the comparison table
 bins <-  data.frame( rbind(h.H2O2$density,h.aging$density) )  ;
 bins3 <-  data.frame( rbind(h.H2O2$density,h.aging$density, h.H2O2Mean$density) )  ;
 my.mids = my.breaks[-length(my.breaks)] + mystep/2
 #my.mids
 names( bins ) <- my.mids
 row.names(bins) <- c( "H2O2", "Chronological Aging" )
 row.names(bins3) <- c( "H2O2", "Chronological Aging", "H2O2 Mean" )
 bins
 bins3
 
 pdf("plots/Figure_sideBYside20131125.pdf", width=8, height=5)
 barplot( as.matrix(bins), beside=T, col=c("black","gray"), ylab="Relative Frequency", xlab="Ratios",
          legend= c( "Cb/Cv H2O2", "Tg/Tc Aging" )  );
 title(main="H2O2 and chronological aging elevate LOH at different modes" )
 dev.off(); 
  
# barplot( as.matrix(bins3), beside=T, col=c("black","gray",'red'), ylab="Relative Frequency", xlab="Ratios",
#         legend= c( "Cb/Cv H2O2", "Tg/Tc Aging",  "Cb/Cv H2O2 Mean" )  );
 
 
### regression analysis 
pTb = 1: length(tb3a[1,])
names(pTb) = names(tb3a)
 for( j in c(3:36) ) {
   m = lm( tb3a[, j] ~ tb3a$ARLS, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 #ARLS          Lmax        L0.all     CLS.vs.Tc      Cv.vs.Cb      Cb.vs.Cv 
 #3.546097e-02  3.877701e-02  3.137392e-02  4.748623e-02  0.000000e+00  7.306331e-05 
 #Cb.vs.Cv.mean 
 #5.729385e-04
 pTb[pTb>0.05 & pTb<0.15]

 ### regression analysis 
 pTb = 1: length(tb3b[1,])
 names(pTb) = names(tb3b)
 for( j in c(3:14,16:19) ) {
   m = lm( tb3b[, j] ~ tb3b$Cv.vs.Cb, na.rm=T)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 #ARLS          Lmax        L0.all     CLS.vs.Tc      Cv.vs.Cb      Cb.vs.Cv 
 #3.546097e-02  3.877701e-02  3.137392e-02  4.748623e-02  0.000000e+00  7.306331e-05 
 #Cb.vs.Cv.mean 
 #5.729385e-04
 pTb[pTb>0.05 & pTb<0.15]
 
 
 
 tiff("ARLS-CvCb-20130618.tif",width=480,height=480)
 par(font=2) 
 plot( tb$ARLS ~ tb$Cv.vs.Cb , pch=19, col="red", main="H2O2-LOH ~ ARLS, 20130618", ylim=c(22,38),xlim=c(0.5,4.2)
       , ylab='ARLS',xlab='Cv/Cb Tolerance to H2O2-induced genomic instability')
 text( tb$Cv.vs.Cb+0.08, tb$ARLS+0.5, tb$strain)
 m = lm(tb$ARLS ~ tb$Cv.vs.Cb  )
 abline( m, col="blue")
 summary(m)
 text(3.5, 35,  "R2=0.37 p=0.035")
 dev.off()
 
 
 tiff("L0-CvCb-20130819.tif",width=480,height=480)
 par(font=2)
 plot( tb$L0.all ~ tb$Cv.vs.Cb , pch=19, col="red", main="H2O2-LOH ~ mitotic asymmetry, 20130618"
       ,xlim=c(0.5,4.2),ylim=c(0.27, 0.045)
       , ylab='L0.all mitotic asymmetry',xlab='Cv/Cb Tolerance to H2O2-induced genomic instability')
 text( tb$Cv.vs.Cb+0.09, tb$L0.all+0.008,tb$strain)
 m = lm( tb$L0.all ~ tb$Cv.vs.Cb  )
 abline( m, col="blue")
 summary(m)
 text(3.5, 0.25, "R2=0.42 p=0.031")
 dev.off()
 
 
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
 
###########2014 Nov 14, corAR1()
require(nlme)
 
 smalltb = tb[, c("R0","G","Cv.vs.Cb")]
 smalltb = smalltb[! is.na(smalltb$Cv.vs.Cb), ]
 head(smalltb)
 # Adding autocorrelation
 #mod4 = gls(youth ~ cadult*minwage, correlation = corAR1(form=~1), data = un)
 mod4 = gls(Cv.vs.Cb ~ log10(R0)*G, correlation = corAR1(), data = smalltb)
 
 summary(mod4)
 
              
#quit("yes")

####
### END
####
 
