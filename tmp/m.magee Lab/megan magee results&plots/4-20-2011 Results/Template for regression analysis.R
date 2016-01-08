tb = read.csv("mutation.csv") #read in the data

tb$MuR = tb$White / (tb$White + tb$Red) #calculate mutation rate
tb$Density = (tb$White + tb$Red) * 10 * tb$Dilution #calculate cell density

groups = sort(unique(tb$Group)); #how many groups
myReport = data.frame(groups); #generate a class report

#the follow codes generate report for the class
for( i in 1:6) {
 print(i)
 sub = tb[tb$Group==groups[i], ];
 myReport$MeanMu[i] = mean(sub$MuR, na.rm=T);
 myReport$StdMu[i] = sd(sub$MuR, na.rm=T);

 myReport$MeanDensity[i] = mean(sub$Density, na.rm=T)/1E6;
 myReport$StdDensity[i] = sd(sub$Density, na.rm=T)/1E6;
 
 myReport$OD[i] = sub$OD[1]
}

myReport
pairs(myReport)

#Do the density -OD plot. Q: Are the two cell density measures consistent? 
plot( myReport$MeanDensity ~ myReport$OD, col='white', main='Density ~ OD', xlab="OD600nm", ylab="Density xE6/ml" );
text(myReport$OD, myReport$MeanDensity, myReport$groups);
m = lm(myReport$MeanDensity ~ myReport$OD);
abline(m, col="red");
summary(m)
text( 1,1000, "R2=0.07, p=0.611");


# Do the mutation rate ~ OD plot
plot( myReport$MeanMu ~ myReport$OD, col='white', main='Mutation ~ OD', xlab="OD600nm", ylab="mutation rate" );
text(myReport$OD, myReport$MeanMu, myReport$groups);
m2 = lm(myReport$MeanMu ~ myReport$OD);
abline(m2, col="red");
summary(m2)
text( 1,0.25, "R2=0.35, p=0.21");







