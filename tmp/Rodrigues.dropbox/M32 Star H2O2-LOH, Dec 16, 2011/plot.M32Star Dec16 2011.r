### need to add standard deviation

rm=(list=ls())

infile = "M32Star Dec16 2011.csv"
tb = read.csv(infile);
tb = read.csv(infile, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")

tb$H2O2 = tb$H2O2stock/2
tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
tb.ori = tb; 

tb$Dilution = tb$Dilution / tb$Dilution[1]

###### some manual curations here
#tb$Dilution[c(19,20,21)] = c(10,10,10)  ##this dilution should be 10 times more
#tb$Black[25] = NA #This number is not right!!!!

######## normalize all data
mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
for ( j in mycolumns) {
 tb[,j] = tb[,j] * tb$Dilution
}

####### find out means
H2O2 = unique( tb$H2O2)
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

###### some manual curations here
tbm$halfBlack[tbm$halfBlack==0 & tbm$H2O2>0] = NA;

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
for ( j in 3:8) {
  tbf[, j] = tbf[,j] / tbf$tot
}

mylabel = 'M32star Dec162011'

pdf("M32star.12162011.black.pdf", width=5, height=5)
plot( tbf$s ~ tbf$H2O2, col="blue", axes=F, xlab='', ylab='', ylim=c(0, 1.2), type='p');
lines( tbf$s ~ tbf$H2O2, col="blue",lty=2)
axis( 4, at=pretty(c(0, 1.2)))
legend ( max(H2O2)*0.7, 0.5, c("viability","black"), col=c("blue","black"), lty=c(2,1), pch=c(1,16) )
par(new=T)
plot( tbf$Black ~ tbf$H2O2, pch=16, xlab='H2O2',ylab="black")
lines(tbf$Black ~ tbf$H2O2)
title(mylabel)
dev.off()

tbf$H2O2[tbf$H2O2==0] = min(H2O2[-1])/10;

with( tbf, plot( tot ~ log10(H2O2), col="blue"));
par(new=T)
with( tbf, plot( Black ~ log10(H2O2), pch=16))

pdf("M32star.12162011.Black-log.pdf", width=5, height=5)
tbf$H2O2[tbf$H2O2==0] = min(H2O2[-1])/10;

with( tbf, plot( s ~ log10(H2O2), col="blue", axes=F, xlab="H2O2", ylab='viability'), ylim=c(-0.5, 1.1) );
xlabels = log10(c(0.001,0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15))
axis(1, at = xlabels, labels= 10^xlabels)
ylabels = c(0, 0.25, 0.5, 0.75, 1.0)
axis(2, at=ylabels, labels=ylabels);
par(new=T)
with( tbf, plot( Black ~ log10(H2O2), pch=16, axes=F, xlab='', ylab=''))
axis(4, pretty(range(tbf$Black)))
title(mylabel)
dev.off()
	

	
#quit("yes")	
   
   