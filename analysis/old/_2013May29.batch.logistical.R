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

ctl.tb = read.csv("_ctl.tb_initial.values.csv")
ctl.tb.input = ctl.tb
files = list.files(path="data.H2O2-LOH")
#ctl.tb2 = data.frame(files)

ctl.tb$update = match(ctl.tb$files, files)

#ii=2
for( file in files ) {   
  rm (fm.s, fm.b )
  #if ( grep("fm.s", ls())  ) { rm (fm.s) }
  #if ( grep("fm.b", ls()) ) { rm (fm.b) }
  ii = match( file, ctl.tb$files)
  infile = ctl.tb$files[ii]
  fullFileName = paste('data.H2O2-LOH/',infile, sep='');
  mylabel = infile
  tb = read.csv(fullFileName, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
  names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")

  ####replace NA with zeros
  mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other"); 
  for ( i in 1:length(tb[,1])) {
    for ( j in mycolumns) {
      if( is.na(tb[i,j]) ) { tb[i,j]= 0; }
    }
  }
    
  
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
  
  tbm = tbm[tbm$tot>1, ] #remove plates with zero colonies
  
  
  ###### calculate fractions
  tbf = tbm; 
  tbf$s = tbf$tot / max(tbf$tot)
  for ( j in 3:8) {
    tbf[, j] = tbf[,j] / tbf$tot
  }  
  tbf$Black[tbf$Black<0]=NA;  #remove weird experimental data, such as low-lead concentration effect
  
  ######## calclulate T0.5 w_s
  t= tbf$H2O2
  s= tbf$s
  
  model1 = function(v,w) {  1/( 1 + ( t / v )^ w ) }  
  fm.s <- gnls( s ~ model1(v,w) , start = list( v = ctl.tb$Cv[ii], w=ctl.tb$CvW[ii])  );
  #    fm.s <- gnls( s ~ model1(v,w) , start = list( v = 0.01, w=1.5)  );
  ctl.tb$Cv[ii]= fm.s$coef[1];
  ctl.tb$CvW[ii] = fm.s$coef[2]
  
  #logistical.black <- function(b.max, b.min, T, w, t ) { ret <- b.max - (b.max - b.min) /( 1 + ( t / T )^ w );  }
  b = tbf$Black
 
  #fill the missing values with imputing
  for( bi in 2: (length(b)-1 ) ) { 
    if( is.na(b[bi]) )  { b[bi] = b[bi-1]/2 + b[bi+1]/2 }
  }
  if( is.na(b[length(b)]) )  { b[length(b)] = b[length(b)-1] }
                              
  #estimate starting values
  b.tmp = sort( b[ b>0 ] );
  b.min = b.tmp[1]/2 + b.tmp[2]/2;
  b.tmp = rev(b.tmp); 
  b.max = b.tmp[ 1 ] /2 + b.tmp[ 2] /2;  #20130529 correct a typo here
  
  modelBlack = function(v, w) { b.max - (b.max-b.min)/(1 + (t/v)^w) }
  #fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = ctl.tb$Cv[ii], w=ctl.tb$CvW[ii])  ); #20130530
  fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = ctl.tb$Cb[ii], w= ctl.tb$CbW[ii])  );

  #fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = fm.s$coef[1], w= fm.s$coef[2]   ) );
  #fm.b <- gnls( b ~ modelBlack(v,w) , start = list( v = 0.1, w= 2 ) );
  
  ctl.tb$Cb[ii]= fm.b$coef[1];
  ctl.tb$CbW[ii] = fm.b$coef[2]
  ctl.tb$b.min[ii] = b.min
  ctl.tb$b.max[ii] = b.max
  
}

write.csv(ctl.tb, '_ctl.tb_out.20130530a.csv', quote=T, row.names=F)
