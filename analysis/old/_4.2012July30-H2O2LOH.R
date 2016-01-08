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
 tb = tb.old[c(1:11), c("strain","ARLS","R0","G","CLS","Tc", "Tg","Tmmax","Tbmax", "Td", "Tdmax","TLmax","Lmax", 
 "b.max", "b.min", "strains", "L0.all", "L0.small" , "Pbt0","Pb0.5t0", "Pbt0.b") ];
 tb$CLS.vs.Tc = tb$CLS / tb$Tc; 
 tb$Tg.vs.Tc = tb$Tg / tb$Tc;
 tb$strain = as.character(tb$strain)

 #load exg06 data
 nat = read.table("062705.rls.cls.tab", sep="\t", header=T, colClasses=c("character", rep(NA,4)) );
 
 #Load H2O2-LOH results
 tb2 = read.csv("H2O2_Log_Plot_Summarized_data,2012Jan24.csv")
 tb2$Strain = as.character(tb2$Strain)
 names(tb2) = c("Date", "Strain", "Cv", "Cb", "OD600nm","Notes","repeat")
 tb2.old = tb2;
 tb2 = tb2[, c(1,2,3,4)]
 tb2$Cv.vs.Cb = tb2$Cv / tb2$Cb
 tb2$Cb.vs.Cv = tb2$Cb / tb2$Cv
 head(tb2)
 
 #check strains names, do they match? 
 strains2 = unique(tb2$Strain)
 intersect( strains2, tb$strain)
 intersect( nat$strain, strains2)
 
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
 t.test( 1 / tb2m$Cv.vs.Cb, mu=1, alternative="less") #p=0.051
 t.test( 1 / tb2m$Cv.vs.CbByMean, mu=1, alternative="less") #p=0.32
 
 t.test( log2(1 / tb2m$Cv.vs.Cb), mu=0, alternative="less") #p=0.022
 t.test( log2(1 / tb2m$Cv.vs.CbByMean), mu=0, alternative="less") #p=0.10
 
 wilcox.test( 1/ tb2m$Cv.vs.Cb, mu=1, alternative="less") #p=0.053
 wilcox.test( 1/ tb2m$Cv.vs.CbByMean, mu=1, alternative="less") #p=0.38
 # Cb/Cv < 1
 
 hist(tb$Tg.vs.Tc, br=10)  
 t.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00072
 wilcox.test( tb$Tg.vs.Tc, mu=1, alternative="greater") #p=0.00098
 #Tg/Tc > 1

 ### side by side bar-plots of Tg/Tc Cb/Cv
 mystep=0.2
 my.breaks = seq( 0.1,  round(max( c( 1/tb2m$Cv.vs.Cb, tb$Tg.vs.Tc ) + 0.1, 1)) ,by= mystep ); 
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
  
### merge tb tb2m
 tb$Cv.vs.Cb = tb2m$Cv.vs.Cb[match(tb$strain, tb2m$strains2)]
 #tb$Cb.vs.Cv = tb2m$Cb.vs.Cv[match(tb$strain, tb2m$strains2)]
 tb$Cb.vs.Cv = 1 / tb$Cv.vs.Cb; 
 #This is for consistency, because 
 # 1/ratio and ratio before and after averaging, give slightly different values. 
 
 tb$Cb = tb2m$Cb[match(tb$strain, tb2m$strains2)]
 tb$Cv = tb2m$Cv[match(tb$strain, tb2m$strains2)] 
 tb$Cb.vs.Cv.mean = tb$Cb / tb$Cv

 plot( tb$Cb.vs.Cv ~ tb$Cb.vs.Cv.mean, xlim=c(0,2),ylim=c(0,2))
 
### regression analysis 
 summary( lm( tb$Cb.vs.Cv ~ tb$ARLS + tb$Tg + tb$Tc + tb$Tg.vs.Tc + tb$CLS ) )
 summary( lm( tb$Cb.vs.Cv ~ tb$ARLS ) )  #p = 0.134
 summary( lm( tb$Cb.vs.Cv.mean ~ tb$ARLS ) )  #p = 0.4
 summary( lm( tb$Cb.vs.Cv ~ tb$R0 + tb$G ) ) # p =0.18
 summary( lm( tb$Cb.vs.Cv ~ tb$R0 * tb$G ) ) # p =0.38
 summary( lm( tb$Cb.vs.Cv ~ sqrt(tb$R0 * tb$G) ) ) # p =0.77
 
 summary( lm( tb$Cb.vs.Cv ~ tb$Tg.vs.Tc) ) #p = 0.4
 summary( lm( tb$Cb ~ tb$Tg.vs.Tc + tb$ARLS + tb$CLS ) ) #p=0.30
 
 summary( lm( tb$Cb ~ tb$CLS ) ) #p=0.23
 summary( lm( tb$Cv ~ tb$CLS ) ) #p=0.85
 
 summary( lm( tb2m$Cb ~ tb2m$Cv ) ) #p0.09
 m = lm( tb2m$Cb ~ tb2m$Cv )
 plot( tb2m$Cb ~ tb2m$Cv )
 abline( m , col='red')
 
 summary( lm( tb$Cb.vs.Cv ~ tb$Tc ) ) #p=0.15
 summary( lm( tb$Cb.vs.Cv.mean ~ tb$Tc ) ) #p=0.15
 
 summary( lm( tb$CLS ~ tb$Cb.vs.Cv ) ) #p=0.024  !!!!! negative !!!why
 #summary( lm( tb$CLS ~ tb$Cb.vs.Cv ) ) #using averaged values, p=0.0705
 
 summary( lm( tb$CLS ~ tb$Cv.vs.Cb ) ) #p=0.013  !!!!! positive, so it is consistent

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
 
 
  
quit("yes")

####
### END
####
 
###2011Feb25: 
 summary( lm( tb$Tg ~ tb$ARLS) )
 summary( lm( tb$Tg ~ tb$G + log(tb$R0)) )
 summary( lm( tb$Tc ~ tb$G + log(tb$R0)) )
 summary( lm( tb$ARLS ~ tb$G + log(tb$R0)) )

 summary( lm( tb$G ~ log10(tb$R0) + tb$Pbt0 + tb$L0.all + tb$Tg.vs.Tc  ) ) 
 

 # robustness to hyper recombination 1/b.max 1/L0
 #summary( lm( tb$b.max ~ log10(tb$R0) + tb$G  ) ) 
 summary( lm( 1/tb$b.max ~ log10(tb$R0) + tb$G  ) ) 
 summary( lm( -tb$b.max ~ log10(tb$R0) + tb$G  ) )  
 summary(lm( - tb$L0.all ~ tb$R0 + tb$G ))
 summary(lm( 1/ tb$L0.all ~ tb$R0 + tb$G ))

 summary(lm( 1/ tb$Pbt0 ~ tb$R0 + tb$G ))
 
 
 summary( lm( tb$Tg.vs.Tc ~ tb$ARLS ) )
 #p=0.008, R2=0.56 

 summary(lm( (tb$ARLS ~ log(tb$R0) + tb$G )))
 #both negaive correlation
 
 summary( lm( tb$Tg.vs.Tc ~ log10(tb$R0)   ) ) #p=0.11
 plot( tb$Tg.vs.Tc ~ log10(tb$R0) )
 
 summary( lm( tb$Tg.vs.Tc ~            tb$G  ) ) #p=0.68, negative 
 plot( tb$Tg.vs.Tc ~ tb$G )
 # negative is opposite to GG realibility model 
 #tb$Tg.vs.Tc[4] = NA; #remove outliers
 summary( lm( tb$Tg.vs.Tc ~            tb$G  ) ) #p=0.79, remove one outlier, it is positive, However, partial correlation is still negative. 
 plot( tb$Tg.vs.Tc ~ tb$G )
  
 summary( lm( tb$Tg.vs.Tc ~ tb$R0 + tb$G   ) )
 #p=0.059, which is similar to ARLS 
 
 summary( lm( tb$Tg.vs.Tc ~ log10(tb$R0) + tb$G  ) )
 #p=0.013, R2=0.66 good p-value  !!!!!
 

 #p=0.22, negative correlations
 
 q("no") 
 
 
 
 
 ######### L0 ~ CLS
#summary( lm( tb$L0.all ~ tb$CLS) ) #p 0.011  R2=0.53 <--

#summary( lm( tb$L0.all ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G ) ) # none
# only CLS shows connection

#summary( lm( tb$L0.small ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G ) ) #none

#summary( lm( tb$L0.all ~ tb$ARLS) ) #p 0.49

#summary( lm( tb$L0.all ~ tb$ARLS + tb$CLS) ) #p 0.05, L0~CLS partial is p=0.02

#summary( lm( tb$L0.all ~ tb$Tc) ) #p 0.33

#summary( lm( tb$L0.all ~ tb$Tg) ) #p 55           


####### b.max ~ ARLS
summary( lm( tb$b.max ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all ) ) #p 0.50

summary( lm( tb$b.max ~ tb$ARLS ))
#Residual standard error: 0.06868 on 9 degrees of freedom
#Multiple R-Squared: 0.4122,     Adjusted R-squared: 0.3469
#F-statistic: 6.311 on 1 and 9 DF,  p-value: 0.0332


######## Lmax
summary( lm( tb$Lmax ~ tb$ARLS ))  #p 0.98


summary( lm( tb$TLmax ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) ) #none

summary( lm( tb$Lmax ~ tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) )
# 0.05 ???  but one degree of freedom?? So, this is questionable.


##### Tg/Tc ~ ARLS
summary( lm( tb$Tg.vs.Tc ~ tb$ARLS) );
# Residual standard error: 0.1861 on 9 degrees of freedom
# Multiple R-Squared: 0.5551,     Adjusted R-squared: 0.5057
# F-statistic: 11.23 on 1 and 9 DF,  p-value: 0.00851



summary( lm( tb$CLS.vs.Tc ~ tb$Tg.vs.Tc ))  #p 0.35

summary( lm( tb$Tg.vs.Tc ~ tb$Lmax + tb$ARLS + tb$Tc + tb$CLS + tb$R0 + tb$G + tb$L0.all + tb$b.max + tb$b.min + tb$Tdmax ) )
## 0.05 ???



postscript("110706.L0-CLS.ps", width=6,height=6, horizontal=F)
m.L0 = lm( tb$L0.all ~ tb$CLS)
plot(tb$L0.all ~ tb$CLS, xlab="CLS",ylab="L0",main="L0 ~ CLS", pch=16, xlim=c(2,18), ylim=c(0.05,0.26) )
abline( m.L0, col="red");

x = numeric( length(tb$CLS) ); names(x) = as.character( tb$strains ); 
y =x;
y[c("M34", "M5",  "M32", "M8")] = c(0.002, -0.01, -0.01, -0.01 )
strains = as.character(tb$strain); names(strains) = tb$strain
#pos = x; pos=NA;
strains["M5"] = "     M5";
text( tb$CLS, tb$L0.all+y+0.004, strains);

dev.off();

#########110606Mon regression:
 summary( lm( ARLS ~ Tg, data=tb) ); #p = 0.36
 summary( lm( ARLS ~ Tc, data=tb) ); # p 0.48
 summary( lm( ARLS ~ Tmmax, data=tb) ); #p 0.59
 summary( lm( ARLS ~ Tbmax, data=tb) ); #p 0.36
 summary( lm( ARLS ~ Td, data=tb) ); #p 0.93

### begin ####ARLS ~ b.max 
 summary( lm( ARLS ~ b.max, data=tb) ); #p 0.06455
 
 long.arls  = tb$b.max[ tb$ARLS > 33  ]
 short.arls = tb$b.max[ tb$ARLS < 33  ]
 t.test( long.arls, short.arls)   # 0.02863
 
 tb.tmp = tb[, c("ARLS", "b.max")] 
 tb.tmp$by.arls = ifelse( tb$ARLS  >33, 1, 0)
 tb.tmp$by.b    = ifelse( tb$b.max > 0.16, 1, 0)

 # wrong:: fisher.test( tb.tmp[,c("by.arls","by.b")] )
 tmp = table( tb.tmp[,c("by.arls","by.b")] )
 fisher.test( tmp ); #p 0.06

 pdf("3c.031908.arls.bmax.pdf", width=6, height=6, )
 m1 = lm( b.max ~ ARLS, data=tb);
 plot( b.max ~ ARLS, data=tb, pch=16, col="black", xlim=c(22,40), ylim=c(0.02, 0.3) );
 abline( m1 , col= "black", lty=2);
 y = tb$b.max +0.01 ;
 names(y) = as.character( tb$strain );
 y[c("M32" )] = y[c("M32" )] - 0.02 
 y[c("YPS163" )] = y[c("YPS163" )] - 0.02 
 x = tb$ARLS +0.5 ;
 x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain);
 text( 25, 0.27, "R2=0.34 p=0.059")
 dev.off();

###end #######ARLS ~ b.max 

### begin ALRS ~ (Tg-Tc)/Tc
 tb$frac.gc = (tb$Tg - tb$Tc) / tb$Tc;
 tb$frac.bm = (tb$Tbmax - tb$Tmmax) / tb$Tmmax; 
 
 summary( lm( ARLS ~ Tg:Tc, data=tb) );
 summary( lm( ARLS ~ frac.gc, data=tb)) # p 0.007
 summary( lm( ARLS ~ frac.bm, data=tb)) # p 0.01

 pdf("3a.031908.Tg-Tc.arls.pdf",width=6,height=6);
 plot( tb$Tg.vs.Tc ~ tb$ARLS, pch=16, col="black", xlim=c(22,40), ylim=c(0.8,1.9),xlab="ARLS",ylab="Tg/Tc" );
 m2 = lm( tb$Tg.vs.Tc ~ tb$ARLS );
 abline( m2 , col= "black", lty=2);
 y = tb$Tg.vs.Tc + 0.01
 names(y) = as.character( tb$strain );
 #y[c("M32" )] = y[c("M32" )] - 0.02 
 x = tb$ARLS ;
 #x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain, pos=4);
 text( 27, 1.7, "R2=0.56, p=0.008")
 dev.off();


 pdf("3b.031908.Trmax-Tmmax.arls.pdf",width=6,height=6);
 m3   = lm( frac.bm ~ ARLS, data=tb);
 tb$Trmax.vs.Tmmax = tb$Tbmax / tb$Tmmax
 m3.2 = lm( Trmax.vs.Tmmax ~ ARLS, data=tb);
 plot( Trmax.vs.Tmmax ~ ARLS, data=tb, pch=16, col="black", xlab="ARLS", ylab="Trmax/Tmmax",
         # main="(Tmax.g-Tmax.m/Tmax.m) ~ ARLS", 
         xlim=c(22,40) );
 abline( m3.2 , col= "black", lty=2);
 y = tb$Trmax.vs.Tmmax + 0.01
 names(y) = as.character( tb$strain );
 x = tb$ARLS ;
 text( x, y, tb$strain, pos=4);
 text( 27,1.4, "R2=0.52, p=0.012");
 dev.off();
### end 

############072307 added changes
 
 #### Tg - Tc
 pdf("072307.Tg.Tc.pdf",width=6,height=6);
 plot( tb$Tg ~ tb$Tc, pch=16,xlab="Tc",ylab="Tg",xlim=c(3,10),ylim=c(4,14) );
 m.gc = lm( tb$Tg ~ tb$Tc );
 abline( m.gc , col= "black", lty=2);
 y = tb$Tg + 0.01
 names(y) = as.character( tb$strain );
 #y[c("M32" )] = y[c("M32" )] - 0.02 
 x = tb$Tc ;
 #x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain, pos=4);
 text( 5, 13, "R2=0.59, p=0.006")
 dev.off();

 ##### Trmax - Tmmax
 pdf("072307.Trmax.Tmmax.pdf",width=6,height=6);
 plot( tb$Tbmax ~ tb$Tmmax, pch=16,xlab="Tmmax",ylab="Trmax",xlim=c(4,12),ylim=c(4,14) );
 m.bm = lm( tb$Tbmax ~ tb$Tmmax );
 abline( m.bm , col= "black", lty=2);
 y = tb$Tbmax + 0.01
 names(y) = as.character( tb$strain );
 #y[c("M32" )] = y[c("M32" )] - 0.02 
 x = tb$Tmmax ;
 #x[c("M1-2")] = x[c("M1-2")] - 2
 text( x, y, tb$strain, pos=4);
 text( 5, 13, "R2=0.60, p=0.005")
 dev.off();

##### Tdmax - Tc
 pdf("072307.Tdmax.Tc.pdf",width=6,height=6);
 plot( tb$Tdmax ~ tb$Tc, pch=16,xlab="Tc",ylab="Tdmax",xlim=c(3.5,10),ylim=c(4,12) );
 m.dc = lm( tb$Tdmax ~ tb$Tc );
 abline( m.dc , col= "black", lty=2);
 y = tb$Tdmax + 0.1
 names(y) = as.character( tb$strain );
 y[c("M14" )] = y[c("M14" )] - 0.4 
 y[c("M2-8" )] = y[c("M2-8" )] + 0.2 
 x = tb$Tc ;
 names(x) = as.character( tb$strain );
 x[c("M14")] = x[c("M14")] - 0.5
 text( x, y, tb$strain, pos=4);
 text( 5, 11, "R2=0.71, p=0.001")
 dev.off();



#quit("yes");
