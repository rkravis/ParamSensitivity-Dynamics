cd(@__DIR__)
using CSV
using DataFrames
using PowerSimulationsDynamics
using PowerSystems
using PowerFlows
using CSV
using LaTeXStrings
using Plots
using LinearAlgebra
using Dates 
using JSON
using JLD2 

const PSY = PowerSystems;
const PSID = PowerSimulationsDynamics;
include("dynamic_data.jl")
include("helper_functions.jl")
cd(@__DIR__)


####################
###### 9 BUS ######
####################


params = generate_9bus_parameter_dataframe()
# Original configuration - SM at bus 1, GFM at bus 2, GFL at bus 3
df, tsim, filepath = run_experiment("data_files/WSCC_9bus.raw", params, ["Bus 2"], ["Bus 3"], ["Bus1"], create_vsm_gfm, create_gfl, dyn_gen_marconato_simple, dyn_gen_marconato, apply_param_update, false, true, "_og");
generate_A_v_D_eig_full_color(filepath, true, (-0.5,0.5), (-0.5,0.5))
generate_parameter_plots(filepath, true)
generate_parameter_plots_damping(filepath,true)
make_9bus_eig_and_pf_plot(filepath, true)


####################
###### 39 BUS ######
####################

params = generate_39bus_parameter_dataframe()

# original configuration 
df, tsim, filepath = run_experiment("data_files/IEEE 39 bus.RAW", params, ["30", "32","33"], ["34", "35","36"], ["37", "38", "39", "31"], create_vsm_gfm, create_gfl, dyn_gen_roundrotor30, dyn_gen_roundrotor30, apply_param_update, false, true,"_og");
#filepath = "IEEE 39 bus_2025-12-08T16:35:28.847_og"
generate_A_v_D_eig_full_color(filepath, true)
generate_parameter_plots(filepath, true)
generate_parameter_plots_damping(filepath,true)
make_39bus_eig_and_pf_plot(filepath,true)
