#flux.jl is provided by Dr. Varner
include("Flux.jl");
using PyPlot
Stoichiometric_matrix= [-1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
-1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 -1.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 1.0 -1.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 2.0 0.0 0.0 2.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0;
0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 -2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0;
0.0 0.0 0.0 0.0 1.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0;
0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 0.0 -1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0;
0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0];

#given parameters
Ll = 308; #protein length (aa)
Lx = 924; #gene length (nt)
Gc = 5*10e-3; #plasmid concentratioin (uM)
CV = 15; #volume of cell(uL)
Rx = 0.15; #RNA polymerase concentration (uM)
RL = 1.6; #ribosome concentration (uM)
vx = 60; #transcription elongation rate (nt/s)
vl = 16.5; #translation elongation rate (aa/s)
Kx = 0.3; #saturation constant transcription (uM)
Kl = 57; #saturation constant translation (uM)
taux = 2.7; #time constant transcription
taul = 0.8; #time constant translation
kdx = 8.35/3600; #degradation constant mRNA (/s)
kdl = 9.9*10e-3/3600; #degradation constant protein (/s)
Lx = 1000; #characteristic gene length (nt)
Lt = 330; #characteristic gene length (aa)
w1 = 0.26;
w2 = 300.0;
k = 0.3*10e3;
n = 1.5

I = collect(0.1:10:10000);
r_l = zeros(length(I));


#to do optimization
for i in 1:length(I)
  fI = (I[i]^n)/((I[i]^n)+k^n);
  u = (w1+w2*fI)/(1+w1+w2*fI);
  v2 = u*(vx/Lx)*Rx*((Gc)/(Kx*taux+(1+taux)*Gc))
  mRNA = v2/kdx
  v5 = (vl/Ll)*RL*(mRNA/(taul*Kl+(1+taul)*mRNA));#Bounds given as -100000.0≤b≤100000.0
  Default_bounds_array = [0	100000/3600;
  v2	v2;
  0	100000/3600	;
  0	100000/3600	;
  0	v5;
  0	100000/3600 ;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600;
  -100000/3600 100000/3600];
  species_bounds_array = zeros(17,2)
  Objective_coefficient_array = [0.0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
  optimize = calculate_optimal_flux_distribution(Stoichiometric_matrix, Default_bounds_array, Species_bounds_array, Objective_coefficient_array)
  r=optimize[2];
  r_l[i] = r[5];
end

proteinC = r_l./kdl;
I = I./10e3;

figure()
semilogx(I,proteinC)
xlable("Inducer[mM]")
ylable("Protein[uM]")
savefig("figure2.png")
