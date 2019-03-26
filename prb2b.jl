include("parameters.jl")
include("problem2a.jl")

#Set initial values for x and I
m1 = 0.0
m2 = 0.0
m3 = 0.0
p1 = 0.0
p2 = 0.0
p3 = 0.0
I = 0.0
#setup sensitivity matrix
sens_phase1 = zeros(6,23);
sens_phase2early = zeros(6,23);
sens_phase2late = zeros(6,23);
#setup changes of the parameters

Parameters=[KL;
Kx;
vx;
vl;
kd_m;
kd_p
W_I1;
W_12;
W_13;
W_23;
W_11;
k_I1;
n_I1;
k_12;
n_12;
k_13;
n_13;
k_23;
n_23;
Lx_1;
Lx_2;
Lx_3
]

#h step size
h=0.5*[KL;
Kx;
vx;
vl;
kd_m;
kd_p
W_I1;
W_12;
W_13;
W_23;
W_11;
k_I1;
n_I1;
k_12;
n_12;
k_13;
n_13;
k_23;
n_23;
Lx_1;
Lx_2;
Lx_3
]

#Things need to be changed
KL=0.5*KL
Kx=0.5*Kx
W_I1=0.5*W_I1
W_12=0.5*W_12
W_13=0.5*W_13
W_23=0.5*W_23
W_11=0.5*W_11
k_I1=0.5*k_I1
n_I1=0.5*n_I1
k_12=0.5*k_12
n_23=0.5*n_23
k_23=0.5*k_23
Lx_1=0.5*Lx_1
Lx_2=0.5*Lx_2
Lx_3=0.5*Lx_3

#set up time range
trange=(1:1:560)

#output array
x_i=zeros(length(trange),6);

#A matrix
A0=[-kd_m-miu 0 0 0 0 0;0 -kd_m-miu 0 0 0 0; 0 0 -kd_m-miu 0 0 0; 0 0 0 -kd_p-miu 0 0; 0 0 0 0 -kd_p-miu 0; 0 0 0 0 0 -kd_p-miu]

#S matrix
S0=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]

#f values
f_I1=(I^n_I1)/(k_I1^n_I1+I^n_I1);
f_12=(p1^n_12)/(k_12^n_12+p1^n_12);
f_13=(p1^n_13)/(k_13^n_13+p1^n_13);
f_23=(p2^n_23)/(k_23^n_23+p2^n_23);

#find control functions
u_m1 = (W_11 + W_I1*f_I1)/(1+W_11 + W_I1*f_I1);
u_m2 = (W_22 + W_12*f_12)/(1+W_22 + W_12*f_12);
u_m3 = (W_33 + W_13*f_13)/(1 + W_33 + W_23*f_23 + W_13*f_13);
#control function is assmed 1 for proteins
u_p1=u_p2=u_p3=1;

#r values
r_x1=kE_m1*u_m1*RNAP_gDW*gene_gDW/(Kx*tau_m1+(tau_m1+1)*gene_gDW)
r_x2=kE_m2*u_m2*RNAP_gDW*gene_gDW/(Kx*tau_m2+(tau_m2+1)*gene_gDW)
r_x3=kE_m3*u_m3*RNAP_gDW*gene_gDW/(Kx*tau_m3+(tau_m3+1)*gene_gDW)
r_l1=kE_p1*rib_gDW*m1/(KL*tau_p1+(tau_p1+1)*m1)
r_l2=kE_p2*rib_gDW*m2/(KL*tau_p2+(tau_p2+1)*m2)
r_l3=kE_p3*rib_gDW*m3/(KL*tau_p3+(tau_p3+1)*m3)

#creat a matrix for r
r=[r_x1; r_x2; r_x3; r_l1; r_l2; r_l3]

#I think I was able to setup the initial conditions. From here,
