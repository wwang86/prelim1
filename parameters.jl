#Problem 2 parameters

DT = 0.66*60; #doubling time (min)
copies_per_cell = 200; #copy number of plasmids (copies/cell)
Lx_ch=1000; #characteristic length (nt)
LT_ch=333; #characteristic length (AA)
Rib = 45000; #copy numbers of ribosome (copies/cell)
Lx_1 = 1200; #nt
Lx_2 = 2400; #nt
Lx_3 = 600; #nt
Ll_1 = 400; #aa
Ll_2 = 800; #aa
Ll_3 = 200; #aa
M = 2.8e-13; #cell weighs (g/cell)
water_fraction = 0.7;
Dry_w = M*(1-water_fraction); #dry weight of cell (g/cell)
halflife_x = 2.1; #half_life of mRNA (min)
halflife_l = 24*60; #half_life of protein (min)
vx = 60*60; #RNAP elongation rate (nt/min)
vl = 16.5*60; #ribosome elongation rate (aa/min)
RNAP = 1150; #copy numbers of RNAP (copies/cell)
RNAPn = RNAP/(6.02e23)*10.0^9; #convert to nmol
RNAP_gDW = RNAPn/Dry_w
gene_gDW = copies_per_cell/(6.05e23)*10.0^9/Dry_w
rib_gDW = Rib/(6.05e23)*10.0^9/Dry_w
kd_m = log(2)/halflife_x
kd_p = log(2)/halflife_p
miu = log(2)/DT
#calculate kEj for each step
ke_mchar = vx/Lx_ch
ke_pchar = vl/LT_ch
Lx_factor1 = Lx_ch/Lx_1
Lx_factor2 = Lx_ch/Lx_2
Lx_factor3 = Lx_ch/Lx_3
Lp_factor1 = LT_ch/Ll_1
Lp_factor2 = LT_ch/Ll_2
Lp_factor3 = LT_ch/Ll_3
kE_m1 = ke_mchar*Lx_factor1
kE_m2 = ke_mchar*Lx_factor2
kE_m3 = ke_mchar*Lx_factor3
kE_p1 = ke_pchar*Lp_factor1
kE_p2 = ke_pchar*Lp_factor2
kE_p3 = ke_pchar*Lp_factor3



#Saturation Constants
Kx = 0.24; #nmol/gDW
KL = 454.64; #nmol/gDW

#tau for mRNA and protein
tau_m1 = 2.7
tau_m2 = 2.7
tau_m3 = 2.7
tau_p1 = 0.8
tau_p2 = 0.8
tau_p3 = 0.8


#Promoter Weights
W_I1 = 100;
W_11 = 0.000001;
W_12 = 10.0;
W_13 = 5.0;
W_22 = 0.000001;
W_23 = 50.0;
W_33 = 0.000001;

#Promoter binding
k_I1 = 0.30;
n_I1 = 1.5;
k_12 = 1000.0;
n_12 = 1.5;
k_13 = 1000.0;
n_13 = 1.5;
k_23 = 1000.0;
n_23 = 10.0;
