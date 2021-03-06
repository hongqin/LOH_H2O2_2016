# 2013 Dec 9, using results from _3b.2013Dec9.summarizeFitting.H2O2LOH.R

# 2013Nov25, use summarized results from _3.*R, which contains standard deviations.  
# 2008 CLS data uses M14, 2013H2O2 uses SGU57, so only 10 strains overlaps

# Use 2013 fitting results and gnls reults
# 2013Nov14, use correlated regression
# Reverse L0 axis in plot, 2012 July 30

# 2012Feb25, Tg.vs.Tc ~ ln(R0) + G
# Partial correlations are all negative, what does this mean? Of couse, they should be negative correlations. 
# Tg/Tc is a measure of ability to maintian recombiation rate during aging

 rm( list = ls() );
 setwd("~/github/LOH_H2O2_2016/analysis")
 
###############
 #Load 2013 H2O2-LOH results
 list.files(pattern="csv", path='output')
 tb2 = read.csv("output/LOHH2O2_averaged20131210_v1.csv") 
 
##### #Load previous LOH-CLS results, Qin Plos One 2008
# tb = read.table("summary.new.by.strain.csv", header=T, sep="\t");
 tb = read.table("021307.summary.by.strain.csv", header=T, sep="\t");
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
 tb2$strain = tb2$strains
 #tb3a = merge(tb, tb2, by='strain', all=T)
 tb3 = merge(tb, tb2, by='strain')
 tb3 = tb3[ , -(grep('strains', names(tb3))) ]
 
 #remove BY4743, will increae p-value
 #tb3 = tb3[ tb3$strain != "BY4743", ]
 
 summary(lm(tb3$Cv.vs.Cb ~ tb3$ARLS)) #p 0.078
 summary(lm(tb3$Cb.vs.Cv ~ tb3$ARLS)) #p 0.039

 #summary(lm(  tb3$ARLS ~ tb3$CvManu * tb3$CbManu )) 
 
 ### regression analysis 
 pTb = 1: length(tb3[1,])
 names(pTb) = names(tb3)
 for( j in c(2:41) ) {
   #m = lm( tb3a[, j] ~ tb3a$ARLS, na.rm=T)
   m = lm( tb3[, j] ~ tb3$Cb.vs.Cv, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.05]
 
 #ARLS       L0.all    CLS.vs.Tc        CvSTD     bminH2O2     bmaxH2O2     Cv.vs.Cb Cv.vs.Cb.STD     Cb.vs.Cv 
 #0.0392109400 0.0078946548 0.0354310798 0.0757766416 0.0139759187 0.0047376191 0.0001815387 0.0362846069 0.0000000000 
 
 pTb = 1: length(tb3[1,])
 names(pTb) = names(tb3)
 for( j in c(2:41) ) {
   #m = lm( tb3a[, j] ~ tb3a$ARLS, na.rm=T)
   m = lm( tb3[, j] ~ tb3$Cb0.5.vs.Cv, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb[pTb<0.1]
 
 
 pTb2 = 1: length(tb3[1,])
 names(pTb2) = names(tb3)
 for( j in c(2:38) ) {
   m = lm( tb3[, j] ~ ( tb3$Cb0.5Mean), na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb2[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb2[pTb2<0.05]
 # CvManu    CvManuSTD       CvMean        CvSTD       CbMean        CbSTD     Cv.vs.Cb Cv.vs.Cb.STD    Cb0.5Mean 
 # 0.0050616007 0.0005836548 0.0074151552 0.0007418509 0.0002412994 0.0035395976 0.0057037439 0.0002793327 0.0000000000 
 
 pTb2 = 1: length(tb3[1,])
 names(pTb2) = names(tb3)
 for( j in c(2:38) ) {
   m = lm( tb3[, j] ~ ( tb3$CbManu), na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb2[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb2[pTb2<0.075] 
 # L0.small     CvMean     CbManu     CbMean 
 # 0.01231221 0.03964427 0.00000000 0.03873183 
 
 
 pTb2 = 1: length(tb3[1,])
 names(pTb2) = names(tb3)
 for( j in c(2:41) ) {
   m = lm( tb3[, j] ~ ( tb3$CvMean), na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CvMean, na.rm=T)
   #m = lm( tb3a[, j] ~ tb3a$CbMean, na.rm=T)
   sm = summary(m)
   pTb2[j] = 1 - pf(sm$fsta[1], sm$fsta[2], sm$fsta[3])
 }
 pTb2[pTb2<0.075] 

 

 ### side by side bar-plots of Tg/Tc Cb/Cv
 mystep=0.2
 my.breaks = seq( 0.2,  round(max( c( tb3$Cb.vs.Cv, tb3$Tg.vs.Tc )) + 0.2, 1) ,by= mystep ); 
# my.breaks3 = seq( 0.1,  round(max( c( tb3$Cb.vs.Cv, tb$Tg.vs.Tc ) + 0.1, 1)) ,by= mystep ); 
 h.H2O2  <- hist( tb3$Cb.vs.Cv, br= my.breaks, xlab = "Cb/Cv", ylab = "relative density", freq=F ) ; 
# h.H2O2  <- hist( tb3$Cb.vs.Cv, br= 10, xlab = "Cb/Cv", ylab = "relative density", freq=F ) ; 
 
 h.aging <- hist(tb$Tg.vs.Tc, br= my.breaks, xlab = "Tg/Tc",  ylab = "relative density", freq=F ) ;
 
 #generate the comparison table
 bins <-  data.frame( rbind(h.H2O2$density,h.aging$density) )  ;
# bins3 <-  data.frame( rbind(h.H2O2$density,h.aging$density, h.H2O2Mean$density) )  ;
 my.mids = my.breaks[-length(my.breaks)] + mystep/2
 #my.mids
 names( bins ) <- my.mids
 row.names(bins) <- c( "H2O2", "Chronological Aging" )
# row.names(bins3) <- c( "H2O2", "Chronological Aging", "H2O2 Mean" )
 bins
# bins3
 
 #pdf("plots/Figure_sideBYside20131209.pdf", width=8, height=5)
 tiff("plots/Figure_sideBYside20131209.tif", width=480, height=480)
 barplot( as.matrix(bins), beside=T, col=c("black","gray"), ylab="Relative Frequency", xlab="Ratios",
          legend= c( "Cb/Cv H2O2", "Tg/Tc CLS" )  );
 title(main="H2O2 and CLS trigger LOH at different modes" )
 dev.off(); 

 ks.test(tb3$Cb.vs.Cv, tb3$Tg.vs.Tc) #p=0.031
 
  
 tiff("plots/ARLS-CbCv-20131209.tif",width=480,height=480)
 par(font=2) 
 plot( tb3$ARLS ~ tb3$Cb.vs.Cv , pch=19, col="red", main="H2O2-LOH ~ ARLS, 20131209", ylim=c(22,38), xlim=c(0.1, 2.5)
       , ylab='ARLS',xlab='Cb/Cv Tolerance to H2O2-induced genomic instability')
 text( tb3$Cb.vs.Cv+0.08, tb3$ARLS+0.5, tb3$strain)
 m = lm(tb3$ARLS ~ tb3$Cb.vs.Cv  )
 abline( m, col="blue")
 summary(m)
 text(1.75, 28,  "R2=0.36 p=0.039")
 dev.off()
 
 
 tiff("plots/L0-CbCv-20131209.tif",width=480,height=480)
 par(font=2)
 plot( tb3$L0.all ~ tb3$Cb.vs.Cv , pch=19, col="red", main="H2O2-LOH ~ mitotic asymmetry, 20131209"
       ,xlim=c(0.1,2.5),ylim=c(0.28, 0.05)
       , ylab='L0.all mitotic asymmetry',xlab='Cb/Cv Tolerance to H2O2-induced genomic instability')
 text( tb3$Cb.vs.Cv+0.09, tb3$L0.all+0.008,tb3$strain)
 m = lm( tb3$L0.all ~ tb3$Cb.vs.Cv  )
 abline( m, col="blue")
 summary(m)
 text(0.8, 0.25, "R2=0.56 p=0.008")
 dev.off()
 
 
 ### Cv/Cb or Cb/Cv ~ robustness? I need a positive proxy
 #summary(lm( tb3$Cv.vs.Cb ~ tb3$G ) )  #positive, p=0.20, 
 #summary(lm( tb3$Cb.vs.Cv ~ tb3$G ) )  #negative, p=0.25
 
 
 
 
 
              
#quit("yes")

####
### END
####
 
