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
include("parameter_generation.jl")
cd(@__DIR__)

update_fcn, param_df = gfm_darco_gfl_pscad()

df, tsim, filepath = run_experiment("data_files/WSCC_9bus.raw", param_df, ["Bus 2"], ["Bus 3"], ["Bus1"], create_vsm_gfm, create_gfl, dyn_gen_marconato_simple, dyn_gen_marconato, update_fcn, false, true, "update_nov_25");

generate_A_v_D_eig_full_color(filepath, true,(-0.5,0.5), (-0.5,0.5))
generate_parameter_plots(filepath, true)
generate_parameter_plots_damping(filepath, true)

