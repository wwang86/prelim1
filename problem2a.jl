#This code is adapted from the one given by Dr. Varner for hw3

# Sample script to run a model using GRNSimKit -
using GRNSimKit
using PyPlot

# time step -
time_step_size = 1.0*(1/60)

path_to_model_file = "/Users/12618/Desktop/ICT1.json"

# Build a data dictionary from a model file -
ddd = build_discrete_dynamic_data_dictionary_from_model_file(time_step_size, path_to_model_file)

# Run the model to steady-state, before we do anything -
steady_state = GRNSteadyStateSolve(ddd)

# Run the model 100 time units *before* we add inducer -
ddd[:initial_condition_array] = steady_state
(T0, X0) = GRNDiscreteDynamicSolve((0.0,1.0,time_step_size), ddd)

# Add inducer T = 1 mM for 100 time units -
ddd[:initial_condition_array] = X0[end,:]
ddd[:initial_condition_array][7] = 10.0
tstart_1 = T0[end]
tstop_1 = tstart_1 + 5.0
(T1, X1) = GRNDiscreteDynamicSolve((tstart_1,tstop_1, time_step_size), ddd)

# # Wash inducer out -no inducer wash out step in this problem
#ddd[:initial_condition_array] = X1[end,:]
#ddd[:initial_condition_array][7] = 0.0
#tstart_2 = T1[end]
#tstop_2 = tstart_2 + 10.0
#(T2, X2) = GRNDiscreteDynamicSolve((tstart_2,tstop_2, time_step_size), ddd)

# Package -
T = [T0 ; T1]
X = [X0 ; X1]

# make a plot -
plot(T,X[:,4]*1e-3,linewidth=2, linestyle="--")
plot(T,X[:,5]*1e-3,linewidth=2, linestyle="--")
plot(T,X[:,6]*1e-3,linewidth=2, linestyle="--")

# axis -
xlabel("Time (hr)", fontsize=16)
ylabel("Protein (nmol/gDW)", fontsize=16)
savefig("cheme.png")
