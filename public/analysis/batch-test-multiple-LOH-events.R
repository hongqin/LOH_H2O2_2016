

# for some datasets, the 3/4 counts are too low. 

rm=(list=ls())
setwd("~/github/LOH_H2O2_2016/analysis")
debug = 0;

FileList = list.files( path="../data.H2O2-LOH/");  FileList; 
output = data.frame(FileList)
output$p_ttest = NA; 

if( debug > 5) {FileList = FileList[1:3]}

for( ii in 1:length(FileList)) {  
  infile = FileList[[ii]] #for list
  print( paste("ii= ", ii ))
  
  fullFileName = paste('../data.H2O2-LOH/',infile, sep='');
  print( fullFileName )
  mylabel = infile
  
  tb = read.csv(fullFileName, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
  names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
  
  ######## set zeros
  mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other"); 
  for( i in 1:length(tb[,1])) {
    for ( j in mycolumns) {
      if( is.na(tb[i,j]) ) { tb[i,j]= 0 }
    }
  }
  
  # Remove ZERO counts for 3/4 black analysis, 2016March22
  tb$Black[tb$Black<=0] = NA
  tb$halfBlack[tb$halfBlack<=0] = NA
  tb$quarterBlack[tb$quarterBlack<=0] = NA
  tb$ThreeQBlack[tb$ThreeQBlack<=0] = NA
  
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
  #s = H2O2
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
  
  tryCatch(
    { 
      H0TQB<-tbf$halfBlack*tbf$quarterBlack
      tt = t.test( tbf$ThreeQBlack, H0TQB, pairwise=T, alternative = "greater")
      print(tt)
      output$p_ttest[ii] = tt$p.value
    }, error = function(e) {e}
  )
  
}

head(output)
summary(output)
output$flag = '';
output$flag = ifelse( output$p_ttest <=0.05, 'signicant', output$flag )
write.csv(output, "__batch_ttest-threeQBlack.csv")
library(xlsx)
write.xlsx(output,"__batch_ttest-threeQBlack.xlsx" )
