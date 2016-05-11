
#ref: 
# http://stackoverflow.com/questions/10508474/step-halving-issue-in-gnlsnlme

require(nlme)
x1 = rnorm(20)
y1 = x1 + rnorm(20)/10

#weird ones
x2 = rnorm(5)+5
y2 = rnorm(5)

y=c(y1,y2)
x=c(x1,x2)
mydata = data.frame(cbind(y,x))

foo = function(a,b) { y = x*a + b }
#gnls( foo, mydata, start=c(a=1,b=0)  )
model1 = gnls( y ~foo(a,b), start=list(a=1,b=0))

model2 = gnls( y ~foo(a,b), data=mydata, start=list(a=1,b=0), weights=varPower(form = ~x))

summary(model1)
summary(model2)

model1
model2
