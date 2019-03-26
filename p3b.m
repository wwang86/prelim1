kex=0.054;
R_T=1.6;
tau=0.8;
K_T=57;
w1=0.26;
w2=300;
K=0.30;
n=1.5;
I=0.0001:0.0001:10;
f_I=I.^n./(K^n+I.^n);
u=(w1+w2.*f_I)./(1+w1+w2.*f_I);
rx=5.88*10^(-5).*u;
mRNA=rx.*3600./8.35;
rl=kex.*R_T.*(mRNA./(tau.*K_T+(tau+1).*mRNA));
protein_c=rl./9.9*10^(-3);
semilogx(I,protein_c)
xlabel('inducer concentration(mM)')
ylabel('protein concentration(uM)')