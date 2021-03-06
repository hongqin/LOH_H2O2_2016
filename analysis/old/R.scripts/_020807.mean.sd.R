#for loh.new.csv, including TLmax and Lmax

 rm( list = ls() );
 tb = read.table( "loh.new.csv", header=T, sep="\t");

 tb$L0 = tb$half.0.raw / tb$b0;

 tb.old = tb;

 tb[as.character(tb$expt) == "M34.100506.37C.tab",    c(1,3:length(tb[1,])) ] = NA;
 tb[as.character(tb$expt) == "YPS128.121205.37C.tab", c(1,3:length(tb[1,])) ] = NA;
 # tb[as.character(tb$expt) == "YPS128,053006",         c(1,3:length(tb[1,])) ] = NA;

 tb = tb[! is.na(tb$strain), ]
 strains = sort( unique( tb$strain) )

parameters =c("strain", "Nloh", "Tc", "Tc.sd", "Tg", "Tg.sd", "Tmmax","Tmmax.sd", "Tbmax","Tbmax.sd", 
"Td", "Td.sd", "Tdmax", "Tdmax.sd", "b.max", "b.max.sd", "b.min", "b.min.sd", "TLmax", "TLmax.sd", "Lmax", "Lmax.sd" );
out= data.frame( matrix(  nrow=length(strains), ncol=length(parameters) ) );
names(out) = parameters;
row.names(out) = strains;

out$strain = strains;

for( ss in 1:length(strains) ) {
  my.Tc = tb$Tc[tb$strain == strains[ss]];
  out$Tc[ss]     = mean( my.Tc );
  out$Tc.sd[ss]  =   sd( my.Tc );
  out$Nloh[ss] = length( my.Tc); 

  my.Tg = tb$Tg[tb$strain == strains[ss]];
  out$Tg[ss]     = mean( my.Tg );
  out$Tg.sd[ss]  =   sd( my.Tg );

  my.Tmmax = tb$Tmmax[tb$strain == strains[ss]];
  out$Tmmax[ss]    =   mean( my.Tmmax);
  out$Tmmax.sd[ss] =     sd( my.Tmmax);

  my.b.max = tb$b.max[tb$strain == strains[ss]];
  out$b.max[ss]    =   mean( my.b.max);
  out$b.max.sd[ss] =     sd( my.b.max);

  my.b.min = tb$b.min[tb$strain == strains[ss]];
  out$b.min[ss]    =   mean( my.b.min);
  out$b.min.sd[ss] =     sd( my.b.min);

  my.Tbmax = tb$Tbmax[tb$strain == strains[ss]];
  out$Tbmax[ss]    =   mean( my.Tbmax);
  out$Tbmax.sd[ss] =     sd( my.Tbmax);

  flag="out";
  # if( ( as.character(strains[ss]) != "RAD52DD" ) && ( as.character(strains[ss]) != "SGU57") ) {
  if( ( strains[ss] != "RAD52DD" ) && ( strains[ss] != "SGU57") ) {
  flag="in";
  my.Td = tb$Td[tb$strain == strains[ss]];
  out$Td[ss]     = mean( my.Td );
  out$Td.sd[ss]  =   sd( my.Td );

  my.Tdmax = tb$Td[tb$strain == strains[ss]];
  out$Tdmax[ss]    =   mean( my.Tdmax);
  out$Tdmax.sd[ss] =     sd( my.Tdmax);

  my.TLmax = tb$TLmax[tb$strain == strains[ss]];
  out$TLmax[ss]     = mean( my.TLmax );
  out$TLmax.sd[ss]  =   sd( my.TLmax );

  my.Lmax = tb$Lmax[tb$strain == strains[ss]];
  out$Lmax[ss]     = mean( my.Lmax );
  out$Lmax.sd[ss]  =   sd( my.Lmax );

  }
 write.table( out, "021307.loh.mean.sd.csv",  quote=F, sep="\t", row.names=F);
}


write.table( out, "021307.loh.mean.sd.csv",  quote=F, sep="\t", row.names=F);


#quit("yes");
