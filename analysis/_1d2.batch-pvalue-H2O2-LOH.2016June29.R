#July 7, fisher exact test, bonferroni correction. 
# for fisher exact test, I use the original counts (before normalization)

#June 29, 2016. Export P(1/2), P(1/4), P(3/4) for Weibiao Wu. 

#batch plot, 2016 Feb 16, testing p-value for 
#if We can estimate p-value of observed ¾ blacks from ¼ and ½ black distributions. 
# assuming multinomial distribution 

# for some datasets, the 3/4 counts are too low. 

rm=(list=ls())
setwd("~/github/LOH_H2O2_2016/analysis")
debug = 0;

FileList = list.files( path="../data.H2O2-LOH/");  FileList; 
output = data.frame(FileList)
output$p_fishertest = NA; 

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
  
  tb$H2O2 = tb$H2O2stock/2
  tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
  tb.ori = tb; 
  tb = tb[ ! is.na(tb$White), ]
  
  tb$Dilution = tb$Dilution / tb$Dilution[1]
  
  ####### merge results by [H2O2] for 3/4LOH test
  H2O2 = sort(unique( tb$H2O2))
  tbMerge = data.frame(cbind(H2O2))
  for ( i in 1:length(H2O2)) {
    c = H2O2[i]
    tmp = tb[ tb$H2O2==c, ]	
    tbMerge$tot[i] = sum(tmp$tot, na.rm=T)
    tbMerge$White[i] = sum(tmp$White, na.rm=T)
    tbMerge$Black[i] = sum(tmp$Black, na.rm=T) 
    tbMerge$halfBlack[i] = sum(tmp$halfBlack, na.rm=T)  
    tbMerge$quarterBlack[i] = sum(tmp$quarterBlack, na.rm=T)
    tbMerge$ThreeQBlack[i] = sum(tmp$ThreeQBlack, na.rm=T)  
    # tbMerge$QQBlack[i] = sum(tmp$QQBlack, na.rm=T);  
  } 
  
  tbMerge$p = NA; 
  
  for( i in 1:length(tbMerge[,1])) {
   if ( (tbMerge$halfBlack[i] >=3) & (tbMerge$quarterBlack[i]>=3) & (tbMerge$ThreeQBlack[i]>=3)) {
    # do fisher exact test
    my.tot = tbMerge$tot[i]; 
    my.P.half= tbMerge$quarterBlack[i];
    my.P.quarter = tbMerge$quarterBlack[i];
    my.P.threeQ = tbMerge$ThreeQBlack[i];
    my.Other = my.tot - (my.P.half + my.P.quarter + my.P.threeQ)
    mytab = matrix(, nrow=2, ncol=2)
    mytab[1,1] = my.P.threeQ;     mytab[1,2] = my.P.half;  #Half True
    mytab[2,1] = my.P.quarter;    mytab[2,2] = my.Other    #Half False
    #1/4 True                       1/4 False
    tryCatch(
      {
        tmp = fisher.test(mytab);
        tbMerge$p[i] = tmp$p.value
      }, error = function(e) {e}
    )
   }
  }
  
  tbMerge = tbMerge[, c("H2O2", "tot", "halfBlack", "quarterBlack", "ThreeQBlack", "p")]
  
  # Bonferronie correction
  alpha = 0.05 #the significant level
  tmp = tbMerge$p; tmp = tmp[!is.na(tmp)]; num_of_trials = length( tmp )
  tbMerge$call = ifelse ( tbMerge$p > alpha/num_of_trials, 'nonsignicant', 'significant')
  
  outfile = paste("output.fisher.test/_fisher_results_", infile, sep='') 
  print (paste("write outfile: ", outfile))
  write.csv(tbMerge, outfile)

  # Now, summaryize the Fisher test results
  # I will describe the range of the p-values, followed by Bonferroni corrections. 
  
  
}#for loop





