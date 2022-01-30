restart
needsPackage "NCAlgebra"
makeMonic = g -> (leadCoefficient g)^(-1)*g
overlap = (a,b,c) -> (a*(b*c%Igb)-(a*b%Igb)*c)%Igb

kk = frac(QQ[c1,c2,c3,c4,c5,c6,c7,c8,c9])
A = kk{x6,x5,x4,x2,x7,x3,x1}

r1 = c1*x1^3+c4*x2*x4*x5+c7*x2*x4*x6
r3 = c2*x3^3+c5*x4*x5*x2+c8*x4*x6*x2
r7 = c3*x7^3+c6*x5*x2*x4+c9*x6*x2*x4
r2 = c4*x4*x5*x1+c7*x4*x6*x1+c5*x3*x4*x5+c8*x3*x4*x6+c6*x4*x7*x5+c9*x4*x7*x6
r4 = c4*x5*x1*x2+c7*x6*x1*x2+c5*x5*x2*x3+c8*x6*x2*x3+c6*x7*x5*x2+c9*x7*x6*x2
r5 = c4*x1*x2*x4+c5*x2*x3*x4+c6*x2*x4*x7
r6 = c7*x1*x2*x4+c8*x2*x3*x4+c9*x2*x4*x7

s6 = r6-(c7/c4)*r5

I = ncIdeal{r1,r2,r3,r4,r5,s6,r7}
Igb = ncGroebnerBasis(I, InstallGB => true)

--- Checking the 1 degree 3 overlap
f1 = overlap(x2,x3*x4,x5)

I = ncIdeal{r1,r2,r3,r4,r5,s6,r7,f1}
Igb = ncGroebnerBasis(I, InstallGB => true)

g1 = overlap(x1^2,x1,x2*x4)
g2 = overlap(x3^2,x3,x4*x5)
g3 = overlap(x7^2,x7,x5*x2)
g4 = overlap(x7*x5,x2,x3*x4)
g5 = overlap(x1,x2*x4,x7*x5)
g6 = overlap(x2*x4,x7*x5,x2)

I = ncIdeal{r1,r2,r3,r4,r5,s6,r7,f1,g1,g2,g3,g4,g5,g6}
Igb = ncGroebnerBasis(I, InstallGB => true)

overlap(x5*x2*x3,x3*x4,x5)
overlap(x7*x5,x2,x4*x7*x5)