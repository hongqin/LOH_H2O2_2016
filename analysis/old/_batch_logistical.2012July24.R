rm(list=ls());
library(nlme)

#setwd("~/projects/LOH-oxidants2012.osX/analysis")
setwd("~/github/LOH_H2O2_2012-master/")

######## functions
logistical.viability <- function( T, w, t ) { ret <- 1 /( 1 + ( t / T )^ w );  }
logistical.mortality <- function( T, w, t ) {  ( w * t^(w-1) ) / (( 1 + (t/T)^w ) * T^w ); }
logistical.black     <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
derivative.black <- function(b.max, b.min, T, w, t) {(b.max - b.min) * w * (t^(w -1)) / ((1 + (t / T )^w)^2 * T^w );}
genome.integrity <- function(b.max, b.min, T, w, t) {  2 * (b.max - b.min) / (1 + (t/T)^w) + (1 - 2 * b.max); }
############ end of functions

files = list.files(path="data.H2O2-LOH")
tb = data.frame(files)

#ctl.file = "../ctrl_file_H2O2LOH2012July24.csv";
#ctl.tb = read.csv(ctl.file)
#ctl.tb = ctl.tb[ctl.tb$Good==1, ]

#ii=2
for( ii in 1: length(tb$files) ) {    
  infile = tb$files[ii]
  fullFileName = paste('data.H2O2-LOH/',infile, sep='');
  mylabel = infile
  tb = read.csv(fullFileName, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
  names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
  
  tb$H2O2 = tb$H2O2stock/2
  tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
  tb.ori = tb; 
  tb = tb[ ! is.na(tb$White), ]
  
  tb$Dilution = tb$Dilution / tb$Dilution[1]
    
  ######## normalize all data
  mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
  for ( j in mycolumns) {
    tb[,j] = tb[,j] * tb$Dilution
  }
  
  ####### find out means
  H2O2 = sort(unique( tb$H2O2))
  tbm = data.frame(cbind(H2O2))
  for ( i in 1:length(H2O2)) {
    c = H2O2[i]
    tmp = tb[ tb$H2O2==c, ]  
    tbm$tot[i] = mean(tmp$tot, na.rm=T)
    tbm$White[i] = mean(tmp$White, na.rm=T)
    tbm$Black[i] = mean(tmp$Black, na.rm=T) 
    tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
    tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
    tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
    tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
  }
  
  ###### calculate fractions

  s= tbf$s
  
  model1 = function(v,w) {  1/( 1 + ( t / v )^ w ) }  
  fm.s <- gnls( s ~ model1(v,w) , start = list( v = ctl.tb$Cv[ii], w=ctl.tb$CvW[ii])  );
  ctl.tb$Cv[ii]= fm.s$coef[1];
  ctl.tb$CvW[ii] = fm.s$coef[2]
  
  #logistical.black <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
  b = tbf$Black
 
  #estimate starting values
  b.tmp = sort( b[ b>0 ] );
  b.min = b.tmp[1]/2 + b.tmp[2]/2;
  b.tmp = rev(b.tmp); 
  b.max = b.tmp[ 2 ] /2 + b.tmp[ 2] /2
  
  modelBlack = function(v, w) { b.max - (b.max-b.min)/(1 + (t/v)^w) }
  fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = ctl.tb$Cv[ii]/2, w=ctl.tb$CvW[ii])  );
  ctl.tb$Cb[ii]= fm.b$coef[1];
  ctl.tb$CbW[ii] = fm.b$coef[2]
  ctl.tb$b.min[ii] = b.min
  ctl.tb$b.max[ii] = b.max
  
}

write.csv(ctl.tb, '_ctl.tb_out.20130528.csv', quote=T, row.names=F)







#########################################
########### old old old old #############



######## calclulate T_g and w_g using Pb
 # Pb = b.max - (b.max-b.min) / (1 + (t/T.g)^w)

 fb = function( t, T.g, w ) { b.max - (b.max - b.min) / (1 + (t/T.g)^w) }
 fb2 = function( t, T.g)  { b.max - (b.max - b.min) / (1 + (t/T.g) ^ww[ss] ); }

  fm.b = gnls( b ~ fb2( t, T.g), start=list( T.g= t2) );

 fit.dPb = derivative.black(b.max, b.min, out2$Tg[ss], out2$wg[ss], tt );  #first point is Inf
 fit.g   = genome.integrity(b.max, b.min, out2$Tg[ss], out2$wg[ss], tt );
 fit.Rb  = fit.dPb / fit.g;     

 out2$Tbmax[ss] = tt[fit.Rb == max(fit.Rb)]

 ######## calculate Td, for half blacks
 out$g.e = genome.integrity(b.max, b.min, out2$Tg[ss], out2$wg[ss], out$t );
 out$R0.5   = out$R0.5.raw / out$g.e;                   ##

 out$dPb = derivative.black(b.max, b.min, out2$Tg[ss], out2$wg[ss], out$t );  #first point is Inf
 out$Rb     = out$dPb / out$g.e;   

 out$L     = out$R0.5 / out$Rb;   

 half.tmp = sort( out$R0.5[ out$R0.5>0 ] );
 half.min = half.tmp[1]/2 + half.tmp[2]/2;  #the first may be an outlier
 half.tmp = rev(half.tmp); 
 half.max = half.tmp[ 1 ] /2 + half.tmp[ 2 ] /2

 half = out$R0.5[ ! is.na(out$R0.5) ]
 t    = out$t[  ! is.na(out$R0.5) ]  ;
 fd = function( t, T.d, w ) { half.max - (half.max - half.min) / (1 + (t/T.d)^w) }
 fd2 = function( t, T.d)  { half.max - (half.max - half.min) / (1 + (t/T.d) ^ww[ss] ); }


