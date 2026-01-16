cd(@__DIR__)
using PowerSimulationsDynamics
using PowerSystems

const PSY = PowerSystems
const PSID = PowerSimulationsDynamics

## Simulation functions 

function add_generator_directly(sys, number, bus_name, Pset, Qset, device_function)
    bus = get_component(ACBus, sys, string(bus_name))
    new_g = create_generator("generator-"*string(bus_name)*"-"*string(number), bus)
    new_g.active_power = Pset/new_g.base_power
    new_g.reactive_power = Qset/new_g.base_power
    add_component!(sys, new_g)
    add_component!(sys, device_function(new_g),new_g)
end

function create_base_system(filename, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function, fast_dynamics_toggle)
    sys = System(filename) 
    for l in get_components(StandardLoad, sys)
        transform_load_to_constant_impedance(l)
    end

    for g in get_components(Generator,sys)
        if g.bus.name in gfm_buses 
            add_component!(sys, gfm_function(g), g)
        elseif g.bus.name in gfl_buses 
            add_component!(sys, gfl_function(g), g)
            #g.dynamic_injector.outer_control.active_power_control.ωz = 0.132 * 2 * pi * 50
            #g.dynamic_injector.outer_control.reactive_power_control.ωf = 0.132 * 2 * pi * 50
            #g.dynamic_injector.freq_estimator.ω_lp = 0.132 * 2 * pi * 50
            #g.dynamic_injector.filter.lg = 0.02
            #g.dynamic_injector.filter.rg = 0.003
            
        elseif g.bus.name in sm_buses 
            add_component!(sys, sm_function(g), g)
        else
            throw("Generator "*g.name*" not assigned as either SM, GFM, GFL")
        end
    end

    if fast_dynamics_toggle
        for l in get_components(Line, sys)
            dyn_branch = DynamicBranch(l)
            add_component!(sys, dyn_branch)
        end
        for d in get_components(DynamicInverter, sys)
            d.filter.ext["is_filter_differential"] = true 
        end
    else 
        for d in get_components(DynamicInverter, sys)
            d.filter.ext["is_filter_differential"] = false 
        end
    end
    return sys 
end

function apply_param_update(sys,df_row,load_dict, line_dict, line_update_fcn, gfm_buses, gfl_buses, sm_buses)
    # Edits the sys object in place 
    # Retrieve parameter values 
    kq = df_row.gfm_kq
    Ta = df_row.gfm_Ta
    kd = df_row.gfm_kd
    cf = df_row.cf
    lf = df_row.lf
    kpv = df_row.gfm_kpv
    kiv = df_row.gfm_kiv
    kpc = df_row.gfm_kpc
    kic = df_row.gfm_kic
    kp = df_row.gfl_kpp
    ki = df_row.gfl_kip
    wz = df_row.wz
    gfm_share = df_row.gfm_share
    inv_share = df_row.inv_share
    load_scale = df_row.load_s
    line_scale = df_row.line_s

    # Scale line using function 
    line_update_fcn(sys, line_dict, line_scale)

    # Scale loads 
    total_p = 0
    total_q = 0
    for l in get_components(PSY.StandardLoad, sys)
        # scale load 
        l.impedance_active_power = load_dict[l.name][1] * load_scale;
        l.impedance_reactive_power = load_dict[l.name][2] * load_scale;
        total_p += l.impedance_active_power*l.base_power
        total_q += l.impedance_reactive_power*l.base_power
    end

    # update power injections 

    # define fractions of each generator type 
    gfm_eta = gfm_share*inv_share
    gfl_eta = (1-gfm_share)*inv_share
    sm_eta = 1-inv_share 

    #xP = find_injections(total_p, gP)
    #xQ = find_injections(total_q, gQ)
    num_GFM = length(gfm_buses)
    num_GFL = length(gfl_buses)
    num_SM = length(sm_buses)
    
    # Adjust injections for each generator 
    for g in get_components(Generator, sys)
        if g.bus.name in gfm_buses 
            g.active_power = gfm_eta*total_p/num_GFM/g.base_power # real power*scale / num of gen at that bus 
            g.reactive_power = gfm_eta*total_q/num_GFM/g.base_power;
            g.dynamic_injector.outer_control.reactive_power_control.kq = kq;
            g.dynamic_injector.inner_control.kpc = kpc;
            g.dynamic_injector.inner_control.kpv = kpv;
            g.dynamic_injector.inner_control.kic = kic;
            g.dynamic_injector.inner_control.kiv = kiv;
            g.dynamic_injector.filter.lf = lf 
            g.dynamic_injector.filter.cf = cf 
            g.dynamic_injector.outer_control.active_power_control.kd = kd;
            g.dynamic_injector.outer_control.active_power_control.Ta = Ta;
        elseif g.bus.name in gfl_buses 
            g.active_power = gfl_eta*total_p/num_GFL/g.base_power; 
            g.reactive_power = gfl_eta*total_q/num_GFL/g.base_power;
            g.dynamic_injector.outer_control.active_power_control.Ki_p = ki; 
            g.dynamic_injector.outer_control.active_power_control.Kp_p = kp;
            g.dynamic_injector.outer_control.reactive_power_control.Ki_q = ki; 
            g.dynamic_injector.outer_control.reactive_power_control.Kp_q = kp; 
            g.dynamic_injector.filter.lf = lf 
            g.dynamic_injector.filter.cf = cf 
            g.dynamic_injector.inner_control.kpc = kpc 
            g.dynamic_injector.inner_control.kic = kic 
            g.dynamic_injector.outer_control.active_power_control.ωz = wz
            g.dynamic_injector.outer_control.reactive_power_control.ωf = wz
        elseif g.bus.name in sm_buses
            g.active_power = sm_eta*total_p/num_SM/g.base_power;
            g.reactive_power = sm_eta*total_q/num_SM/g.base_power;
        end
    end
end

function apply_param_update_droop(sys, df_row, load_dict, line_dict, line_update_fcn, gfm_buses, gfl_buses, sm_buses)

    Rp = df_row.gfm_rp 
    kq = df_row.gfm_kq 
    cf = df_row.cf
    lf = df_row.lf
    kpv = df_row.gfm_kpv
    kiv = df_row.gfm_kiv
    kpc = df_row.gfm_kpc
    kic = df_row.gfm_kic
    kp = df_row.gfl_kpp
    ki = df_row.gfl_kip
    wz = df_row.wz
    gfm_share = df_row.gfm_share
    inv_share = df_row.inv_share
    load_scale = df_row.load_s
    line_scale = df_row.line_s

    # Scale line using function 
    line_update_fcn(sys, line_dict, line_scale)

    # Scale loads 
    total_p = 0
    total_q = 0
    for l in get_components(PSY.StandardLoad, sys)
        # scale load 
        l.impedance_active_power = load_dict[l.name][1] * load_scale;
        l.impedance_reactive_power = load_dict[l.name][2] * load_scale;
        total_p += l.impedance_active_power*l.base_power
        total_q += l.impedance_reactive_power*l.base_power
    end

    # update power injections 

    # define fractions of each generator type 
    gfm_eta = gfm_share*inv_share
    gfl_eta = (1-gfm_share)*inv_share
    sm_eta = 1-inv_share 

    #xP = find_injections(total_p, gP)
    #xQ = find_injections(total_q, gQ)
    num_GFM = length(gfm_buses)
    num_GFL = length(gfl_buses)
    num_SM = length(sm_buses)
    
    # Adjust injections for each generator 
    for g in get_components(Generator, sys)
        if g.bus.name in gfm_buses 
            g.active_power = gfm_eta*total_p/num_GFM/g.base_power # real power*scale / num of gen at that bus 
            g.reactive_power = gfm_eta*total_q/num_GFM/g.base_power;
            g.dynamic_injector.outer_control.reactive_power_control.kq = kq;
            g.dynamic_injector.inner_control.kpc = kpc;
            g.dynamic_injector.inner_control.kpv = kpv;
            g.dynamic_injector.inner_control.kic = kic;
            g.dynamic_injector.inner_control.kiv = kiv;
            g.dynamic_injector.filter.lf = lf 
            g.dynamic_injector.filter.cf = cf 
            g.dynamic_injector.outer_control.active_power_control.Rp = Rp;
        elseif g.bus.name in gfl_buses 
            g.active_power = gfl_eta*total_p/num_GFL/g.base_power; 
            g.reactive_power = gfl_eta*total_q/num_GFL/g.base_power;
            g.dynamic_injector.outer_control.active_power_control.Ki_p = ki; 
            g.dynamic_injector.outer_control.active_power_control.Kp_p = kp;
            g.dynamic_injector.outer_control.reactive_power_control.Ki_q = ki; 
            g.dynamic_injector.outer_control.reactive_power_control.Kp_q = kp; 
            g.dynamic_injector.filter.lf = lf 
            g.dynamic_injector.filter.cf = cf 
            g.dynamic_injector.inner_control.kpc = kpc 
            g.dynamic_injector.inner_control.kic = kic 
            g.dynamic_injector.outer_control.active_power_control.ωz = wz
            g.dynamic_injector.outer_control.reactive_power_control.ωf = wz
        elseif g.bus.name in sm_buses
            g.active_power = sm_eta*total_p/num_SM/g.base_power;
            g.reactive_power = sm_eta*total_q/num_SM/g.base_power;
        end
    end

end

function get_small_signal_results(filename, param_df, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function_a, sm_function_d, param_update_fcn, return_sim=false)

    start_time = time()

    # Make copy of parameter df to store results 
    df = deepcopy(param_df)
    # Create new df columns for results 
    df[!,:stable] .= [(NaN,NaN)]
    df[!,:gfm_eta] .= NaN 
    df[!,:gfl_eta] .= NaN 
    df[!,:sm_eta] .= NaN 
    df[!,:status] .= ""
    df[!,:damping] .= NaN # damping of least stable eigenvalue 
    df[!,:stable_d] .= [(NaN,NaN)]
    df[!,:status_d] .= ""
    df[!,:damping_d] .= NaN
    df[!,:gfm_eta] .= df.gfm_share.* df.inv_share
    df[!,:gfl_eta] .= (1 .-df.gfm_share).*df.inv_share
    df[!,:sm_eta] .= 1 .- df.gfm_eta .- df.gfl_eta
    df[!,:blue_corresponding_a] .= [(NaN,NaN)]
    df[!,:red_corresponding_d] .= [(NaN,NaN)]
    df[!,:diff] .= NaN 
    df[!,:least_damped_a] .= [(NaN,NaN)]
    df[!,:least_damped_d] .= [(NaN,NaN)]

    sys_a = create_base_system(filename, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function_a, false)
    # define perturbation - not important for small signal analysis 
    g_a = first(get_components(Generator,sys_a))
    perturbation_a = ControlReferenceChange(0.025, g_a.dynamic_injector, :P_ref, 0.9);

    sys_d = create_base_system(filename, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function_d, true)
    # define perturbation - not important for small signal analysis 
    g_d = first(get_components(Generator,sys_d))
    perturbation_d = ControlReferenceChange(0.025, g_d.dynamic_injector, :P_ref, 0.9);

    # Generate dicts to store original line and load data 
    load_dict_a = generate_load_dict(sys_a)
    line_dict_a, line_update_fcn_a = generate_line_dict(sys_a)

    load_dict_d = generate_load_dict(sys_d)
    line_dict_d, line_update_fcn_d = generate_line_dict(sys_d)

    # Iterate through parameters 
    last_sim_a = 0
    last_sim_d = 0
    
    for n = 1:size(df)[1];

        param_update_fcn(sys_a, df[n,:], load_dict_a, line_dict_a, line_update_fcn_a, gfm_buses, gfl_buses, sm_buses)
        param_update_fcn(sys_d, df[n,:], load_dict_d, line_dict_d, line_update_fcn_d, gfm_buses, gfl_buses, sm_buses)
        
        # Static model 
        sim_a = PSID.Simulation(
                MassMatrixModel, #Type of model used
                sys_a, #system
                pwd(), #folder to output results
                (0.0,1.0), #time span
                perturbation_a #, #Type of perturbation
                );
        
        # Dynamic model
        sim_d = PSID.Simulation(
                MassMatrixModel, #Type of model used
                sys_d, #system
                pwd(), #folder to output results
                (0.0,1.0), #time span
                perturbation_d #, #Type of perturbation
                );

        # Update df with results 
        # static model 
        ss_a = small_signal_analysis(sim_a);
        ls_eig_a, ls_eig_a_idx = return_least_stable_nonzero_eig(ss_a, sim_a.status)
        df[n,:stable] = (ls_eig_a.re, ls_eig_a.im)
        df[n,:status] = string(sim_a.status);
        df[n,:damping] = -100*real(ls_eig_a)/abs(ls_eig_a);

        # dynamic model 
        ss_d = small_signal_analysis(sim_d);
        ls_eig_d, ls_eig_d_idx = return_least_stable_nonzero_eig(ss_d, sim_d.status)
        df[n,:stable_d] = (ls_eig_d.re, ls_eig_d.im)
        df[n,:status_d] = string(sim_d.status);
        df[n,:damping_d] = -100*real(ls_eig_d)/abs(ls_eig_d);

        # construct diff metric 
        df[n,:diff] = real(ls_eig_a) - real(ls_eig_d)
        # save sim object 
        last_sim_a = sim_a 
        last_sim_d = sim_d

        # find least damped mode in static and dynamic model 
        df[n,:least_damped_a] = find_least_damped_eigenvalue(ss_a.eigenvalues)
        df[n,:least_damped_d] = find_least_damped_eigenvalue(ss_d.eigenvalues)
        
        # match eigenvalues for blue / red points 
        match_up_to = 30 # top X eigenvalues to match against 
        if df[n,:stable][1] <= 0 && df[n,:stable_d][1] > 0 # blue point 
            # find the match for the eigenvalue 
            closest_eig, idx = find_closest_eigenvalue(ss_a, ls_eig_d, match_up_to)
            # retrieve damping ratio for that eigenvalue 
            df[n,:blue_corresponding_a] = (closest_eig.re, closest_eig.im)
            
        end
        if df[n,:stable][1] > 0 && df[n,:stable_d][1] <= 0 # red point 
            # find the match for the eigenvalue 
            match_up_to = 30
            closest_eig, idx = find_closest_eigenvalue(ss_d, ls_eig_a, match_up_to)
            # retrieve damping ratio for that eigenvalue 
            df[n,:red_corresponding_d] = (closest_eig.re, closest_eig.im)
        end
    end
    end_time = time();
    return df, end_time - start_time, last_sim_a, last_sim_d
end

# Wrapper for running multiple simulations 
function run_experiment(filename, params, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function_a, sm_function_d, param_update_fcn=apply_param_update, return_sim=false, save_data=true, abbrv="")

    # Algebraic model first 
    results, tsim, sim_a, sim_d = get_small_signal_results(filename, params, gfm_buses, gfl_buses, sm_buses, gfm_function, gfl_function, sm_function_a, sm_function_d, param_update_fcn, return_sim);

    annotate_status(results)
    filepath = ""
    if save_data
        save_prefix = filename[collect(findlast("/",filename))[1]+1:end-4]*"_"*string(Dates.now())*abbrv
        mkdir(save_prefix)
        save_to_jld2(save_prefix*"/data.jld2",results)
        filepath = save_prefix
        to_json(sim_a.sys, save_prefix*"/alg_sys.json", force=true) # save a copy of the sys so we can verify generator locations if needed
        mdata = Dict(
            "gfm_buses" => gfm_buses, 
            "gfl_buses" => gfl_buses, 
            "sm_buses" => sm_buses, 
            "gfm_fcn" => string(gfm_function),
            "gfl_fcn" => string(gfl_function),
            "sm_a_fcn" => string(sm_function_a),
            "sm_d_fcn" => string(sm_function_d),
            "filename" => filename )
        

        open(filepath*"/output.json", "w") do f
            JSON.print(f, mdata, 4)
        end
    end
    if return_sim 
        return results, tsim, filepath, sim_a, sim_d
    else
        return results, tsim, filepath
    end
    
end

### Results analysis functions ## 

function damping_ratio(x::Vector)
    # Converts vector of tuples where first entry is real part second is imaginary to a vector of corresponding damping ratios
    y = stack(x,dims=1)
    y1 = y[:,1] + y[:,2].*im
    dr = -100 .*real.(y1)./abs.(y1);
    return dr 
end

function find_closest_eigenvalue(sa::PSID.SmallSignalOutput, target::ComplexF64, topX=5)
    # finds closest match to -target- eigenvalue in list of -top X- eigenvalues in -sa- object 
    min_dist = Inf 
    min_idx = NaN 
    max_i = min(topX, length(sa.eigenvalues))
    for i in range(0,max_i-1)
        eig = sa.eigenvalues[end-i]
        dist = norm(eig-target)
        if dist < min_dist 
            min_dist = dist 
            min_idx = i 
        end
    end
    return sa.eigenvalues[end-min_idx], min_idx
end

function find_least_damped_eigenvalue(eigs::Vector)::Tuple 
    # returns the value and index of the least damped eigenvalue in a list of eigenvalues 
    dr = -100 .*real.(eigs)./abs.(eigs);
    dr[ isnan.(dr)] .= Inf # set NaNs to Inf 
    (m, idx) = findmin(dr)
    return real(eigs[idx]), imag(eigs[idx])
end

function return_least_stable_nonzero_eig(ss::PSID.SmallSignalOutput, status::PowerSimulationsDynamics.BUILD_STATUS);
    # Finds the least stable eigenvalue that is nonzero, if the model built correctly
    eig_val = NaN 
    i = 0
    if status != PSID.BUILD_FAILED;
        i = 0
        while real(ss.eigenvalues[end-i]) == 0.0;
            # Choose next largest eig != 0 
            i += 1
        end
        return ss.eigenvalues[end-i], i
    else
        return eig_val, NaN 
    end
end

function verify_parameters(param_df::DataFrame)
    # Prints unique values for each parameter in param dataframe to quickly verify which parameters are being used 
    d = Dict() 
    for c in names(param_df)[1:16]
        d[c] = tuple(unique(param_df[:,c]))
    end
    return d
end

function save_to_jld2(loc, df)
    jldopen(loc, "w") do file
        file["df"] = df     
    end
end

function load_from_folder(filepath::String)
    output = load(filepath*"/data.jld2")
    return output["df"]
    # jldopen(filepath*"/data.jld2") do file
    #     output = file["df"]
    #     end
    # return output 
end

function load_from_folder(filepaths::Vector)
    df = load_from_folder(filepaths[1]*"/data.jld2")
    for i in range(2,length(filepaths))
        df = vcat(df, load_from_folder(filepaths[i]))
    end
    return df 
end     


### Plotting/visualization functions ###

function generate_A_v_D_eig_full_color(filepath::String, save_fig=false, xspan=(-1,1),yspan=(-1,1), xloc=0.65, yloc=0.45,xsize=0.3,ysize=0.45)
    x = load_from_folder(filepath)
    stb = get_stable(x)
    unstb = get_unstable(x)
    margin_blue = get_edge_blue(x)
    margin_red = get_edge_red(x)
    xmin = minimum([minimum(stack(x.stable,dims=1)[:,1]), minimum(stack(x.stable_d,dims=1)[:,1])])
    xmax = maximum([maximum(stack(x.stable,dims=1)[:,1]), maximum(stack(x.stable_d,dims=1)[:,1])])

    #xpts = range(xmin-0.5, xmax*1.05, 15)
    #p = plot(xpts, xpts, fillrange=ones(length(xpts))*maximum(x.stable_d)*1.08, fc=:grey, linealpha=0, alpha=0.3,primary=false,framestyle=:box,legend=:topleft)
    msw_val = 1.0
    p = plot()
    if !isempty(stb)
        scatter!(stack(stb.stable,dims=1)[:,1], stack(stb.stable_d,dims=1)[:,1], color="green", label=L"\lambda_{A,max} < 0, \lambda_{D,max} < 0", msw=msw_val,framestyle=:box)
    end
    if !isempty(unstb)
        Plots.scatter!(stack(unstb.stable,dims=1)[:,1], stack(unstb.stable_d,dims=1)[:,1], color="orange", label=L"\lambda_{A,max} > 0, \lambda_{D,max} > 0",msw=msw_val)
    end
    if !isempty(margin_blue)
        Plots.scatter!(stack(margin_blue.stable,dims=1)[:,1], stack(margin_blue.stable_d,dims=1)[:,1], color="blue", label=L"\lambda_{A,max} < 0, \lambda_{D,max} > 0", msw=msw_val)
    end
    if !isempty(margin_red)
        Plots.scatter!(stack(margin_red.stable,dims=1)[:,1], stack(margin_red.stable_d,dims=1)[:,1], color="red", label=L"\lambda_{A,max} > 0, \lambda_{D,max} < 0", msw=msw_val)
    end

    ylabel!(L"\lambda_{D,max}")
    xlabel!(L"\lambda_{A,max}")

    # Add inset 
    add_inset(xspan, yspan,xloc, yloc,xsize,ysize)
    display(p)
    if save_fig 
        savefig(filepath*"/"*string(Dates.now())*"_full_eig_plot.png")
    end
end

function add_inset(xspan::Tuple, yspan::Tuple,xloc=0.65, yloc=0.45,xsize=0.3,ysize=0.45)
    lens!([xspan[1], xspan[2]], [yspan[1], yspan[2]];
    inset=(1, bbox(xloc, yloc, xsize, ysize)), # bbox(x, y, width, height) where numbers are % of parent 
    subplot=2, ticks=true, framestyle=:box,
    lw=1, ls=:solid)
end

function show_top_participation_factors(sa::PSID.SmallSignalOutput, mode_idx::Int64,plot_on=false, cutoff=nothing)
    # Returns the top participation factor in mode at position -mode_idx- in eigenvalue vector. optionally plots the top -cutoff- participation factors in a bar chart  
    x = sort(summary_participation_factors(sa)[:,[1,end-mode_idx+1]], 2, rev=true)
    print("Eigenvalue = ", sa.eigenvalues[end-mode_idx])
    print(x)
    
    if cutoff==nothing
        cutoff = size(x)[1]
    end
    if plot_on
        bar(x[1:cutoff,1], x[1:cutoff,2], xrotation=90,color=:grey, legend=false,ylabel="participation factor",size=(600,600))
    end
end

function make_damping_plot(filepath::String, save_fig=false)
    y = load_from_folder(filepath)
    p = scatter(damping_ratio(y.least_damped_a), y.diff, color=y.quadrant,primary=false, framestyle=:box)
    vline!([3],color=:pink,width=1, label="3%")
    xlabel!(L"\zeta_{min, static}\ (\%)")
    ylabel!(L"\Delta \lambda")
    display(p)
    if save_fig
        savefig(filepath*"/"*string(Dates.now())*"_damping_ratio_plot.png")
    end
end

function average_parameters(x::DataFrame, cols_to_select::Vector)
    # Averages parameters per column (over rows) for dataframe over selected columns 
    y = x[:,cols_to_select]
    col_names = names(y)
    df = DataFrame([name => [] for name in col_names]) # create new df to store averages 
    means = [mean(y[!, col]) for col in names(y)]
    push!(df, means) 
    return df
end

function generate_parameter_plots(filepath::String, save_fig=false)
    # Plots error as a function of each parameter that is varied, generating one figure per parameter
    if save_fig
        mkpath(filepath*"/parameter_sensitivity_plots")
    end
    df = load_from_folder(filepath)
    x = df#vcat(get_edge_blue(df), get_edge_red(df), get_stable(df))
    all_param_names = names(x)[1:16]
    n_df = size(x)[1]
    print(n_df)
    for c in all_param_names
        # sort by all the other columns, plot points N at a time, moving through the df, where N is the number of unique values of the c parameter 
        n_p = length(unique(x[!,c])) # number of x values
        if n_p != 1 
            remaining = setdiff(all_param_names, [c])
            sort!(x, remaining) # sort in place 
            Plots.plot(framestyle=:box)
            # iterate through in increments of n_p 
            for i in 1:n_p:n_df-n_p-1
                if unique(x[i:i+n_p-1,:quadrant]) != ["orange"] # eliminate lines that are unstable to unstable 
                    Plots.plot!(x[i:i+n_p-1, c], x[i:i+n_p-1,:diff], lw=0.3, markershape=:cross, markercolor=x[i:i+n_p-1, :quadrant], color="gray", alpha=0.4,legend=false)
                end
            end

            xlabel!(label_dict[c])

            ylabel!(L"\Delta \lambda")
            labelfontsize=16;
            ticksize = 14;
            legendsize = 14;
            Plots.plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize, legendtitlesize=legendsize)
            display(Plots.plot!(legend=false))
            if save_fig
                savefig(filepath*"/parameter_sensitivity_plots/"*c*".png")
            end
        end
    end
end

function generate_parameter_plots_damping(filepath::String, save_fig=false)
    # Plots the least damped static mode damping ratio as a function of parameters, generating one figure for each parameter
    if save_fig
        mkpath(filepath*"/parameter_sensitivity_plots_damping")
    end
    df = load_from_folder(filepath)
    x = df#vcat(get_edge_blue(df), get_edge_red(df), get_stable(df))
    all_param_names = names(x)[1:16]
    n_df = size(x)[1]
    for c in all_param_names
        n_p = length(unique(x[!,c])) # number of x values
        if n_p != 1 
            remaining = setdiff(all_param_names, [c])
            sort!(x, remaining) # sort in place 
            Plots.plot(framestyle=:box)
            # iterate through in increments of n_p 
            for i in 1:n_p:n_df-n_p-1
                if unique(x[i:i+n_p-1,:quadrant]) != ["orange"] # eliminate lines that are unstable to unstable 
                    Plots.plot!(x[i:i+n_p-1, c], damping_ratio(x[i:i+n_p-1,:least_damped_a]), lw=0.3, markershape=:cross, markercolor=x[i:i+n_p-1, :quadrant], color="gray", alpha=0.4,primary=false)
                end
            end
            xlabel!(label_dict[c])
            ylabel!(L"\zeta_{min,static}\ (\%)")
            hline!([3],color=:pink,width=2,label="3%")
            labelfontsize=16;
            ticksize = 14;
            legendsize = 14;
            Plots.plot!(xguidefontsize=labelfontsize, yguidefontsize=labelfontsize, tickfontsize=ticksize,legendfontsize=legendsize, legendtitlesize=legendsize,legend=:top)
            display(Plots.plot!())
            if save_fig
                savefig(filepath*"/parameter_sensitivity_plots_damping/"*c*".png")
            end
        end
    end
end

function make_9bus_eig_and_pf_plot(filepath::String, save_fig=false)
    df = load_from_folder(filepath)
    x = get_edge_blue(df)
    df1, tsim, _, sim1a, sim1d = run_experiment("data_files/WSCC_9bus.raw", x[1:2,:], ["Bus 2"], ["Bus 3"], ["Bus1"], create_vsm_gfm, create_gfl, dyn_gen_marconato_simple, dyn_gen_marconato, apply_param_update, true, false);
    sa1a = small_signal_analysis(sim1a)
    sa1d = small_signal_analysis(sim1d)
    xspan = range(0,0.3,15)
    ymax = maximum(abs.(imag(sa1d.eigenvalues)))
    plot(xspan,ones(length(xspan)).*-ymax, fillrange=ones(length(xspan)).*ymax, fillalpha=0.3, fillcolor=:grey,primary=false, lc=:white, framestyle=:box)
    scatter!(real(sa1a.eigenvalues), imag(sa1a.eigenvalues), label=L"\mathrm{static}", alpha=0.8)
    scatter!(real(sa1d.eigenvalues), imag(sa1d.eigenvalues), label=L"\mathrm{dynamic}", alpha=0.8, shape=:utriangle)
    ylabel!(L"\Im \ (\lambda)")
    xlabel!(L"\Re \ (\lambda)")

    #add_inset((-0.1,0.1),(-7,7),0.15,0.3,0.4,0.5)
    add_inset((-1,0.5),(-20,20),0.15,0.3,0.4,0.5)
    if save_fig
        mkpath(filepath*"/eig_comparison_fig")
        savefig(filepath*"/eig_comparison_fig/eig_comparison_blue.png")
        save_to_jld2(filepath*"/eig_comparison_fig/eig_data.jld2", DataFrame(df1[end,:]))
    end
    show_top_participation_factors(sa1d,1,true, 10)
    savefig(filepath*"/pf_blue_point_mode_of_interest.png")
end


function make_39bus_eig_and_pf_plot(filepath::String, save_fig=false)
    df = load_from_folder(filepath)
    x = get_edge_blue(df)
    df1, tsim, _, sim1a, sim1d = run_experiment("data_files/IEEE 39 bus.RAW", x[1:2,:], ["30", "32","33"], ["34", "35","36"], ["37", "38", "39", "31"], create_vsm_gfm, create_gfl, dyn_gen_roundrotor30, dyn_gen_roundrotor30, apply_param_update, true, false);
    sa1a = small_signal_analysis(sim1a)
    sa1d = small_signal_analysis(sim1d)
    xspan = range(0,0.3,15)
    ymax = maximum(abs.(imag(sa1d.eigenvalues)))
    plot(xspan,ones(length(xspan)).*-ymax, fillrange=ones(length(xspan)).*ymax, fillalpha=0.3, fillcolor=:grey,primary=false, lc=:white, framestyle=:box)
    scatter!(real(sa1a.eigenvalues), imag(sa1a.eigenvalues), label=L"\mathrm{static}", alpha=0.8)
    scatter!(real(sa1d.eigenvalues), imag(sa1d.eigenvalues), label=L"\mathrm{dynamic}", alpha=0.8, shape=:utriangle)
    ylabel!(L"\Im \ (\lambda)")
    xlabel!(L"\Re \ (\lambda)")

    #add_inset((-0.1,0.1),(-7,7),0.15,0.3,0.4,0.5)
    #add_inset((-1,1),(-30,30),0.15,0.3,0.4,0.5)
    xlims!(-1,1)
    ylims!(-25,25)
    annotate_box(0,-20,0.8,8)
    annotate_box(0,20,0.8,8)
    
    if save_fig
        mkpath(filepath*"/eig_comparison_fig")
        savefig(filepath*"/eig_comparison_fig/eig_comparison_blue.png")
        save_to_jld2(filepath*"/eig_comparison_fig/eig_data.jld2", DataFrame(df1[end,:]))
    end
    show_top_participation_factors(sa1d,1,true, 10)
    savefig(filepath*"/pf_blue_point_mode_of_interest.png")
end

function annotate_box(x,y,width,height,color=:red)
    xlb = x-width/2
    xub = x+width/2
    ylb = y-height/2
    yub = y+height/2
    plot!([xlb,xub], [ylb, ylb], color=color,primary=false)
    plot!([xlb,xub], [yub, yub], color=color,primary=false)
    plot!([xlb,xlb], [ylb, yub], color=color, primary=false)
    plot!([xub,xub], [ylb, yub], color=color,primary=false)
end


## Parameter generation functions 
function generate_39bus_parameter_dataframe() 

    #GFM VSM outer loop 
    gfm_kq_range = [0.01, 0.4]
    gfm_Ta_range = [0.5, 10]
    gfm_kd_range = [10, 400]

    # GFM inner voltage loop 
    #kpv_range = collect([0.59, 1.0])
    kpv = 0.59 
    v_loop_range = collect([0.001, 0.1])
    #kiv_range = collect([400,736])

    # GFM/GFL inner current loop
    kpc = 0.74
    i_loop_range = [0.1, 0.01]

    # GFL parameters - outer loop 
    kp_range = [0.59, 2]
    ki_range = [7.36, 20]
    wz_range = [0.1*2*pi*60, 2*pi*60]

    # Inverter filter parameters 
    cf = [0.074]
    lf = [0.08]

    # System parameters 
    line_scale_range = [1.0 1.2]
    inverter_share_range = [0.4 0.8]
    gfm_share_range = [0.5 1.0]
    load_scale_range = [0.8 1.0]

    col_names = ["line_s", "load_s", "gfm_Ta", "gfm_kd", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "wz", "inv_share", "gfm_share"]

    df = DataFrame([name => [] for name in col_names])

    foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_Ta_range, gfm_kd_range, gfm_kq_range, kpv, kpv./v_loop_range, kpc, kpc./i_loop_range, kp_range, ki_range, cf, lf, wz_range,inverter_share_range, gfm_share_range))

    return df 
end 

function generate_9bus_parameter_dataframe() 
    
    #GFM VSM outer loop 
    gfm_kq_range = [0.01, 0.4]
    gfm_Ta_range = [0.5, 10]
    gfm_kd_range = [10, 400]

    # GFM inner voltage loop 
    #kpv_range = collect([0.59, 1.0])
    kpv = 0.59 
    v_loop_range = collect([0.001, 0.1])
    #kiv_range = collect([400,736])

    # GFM/GFL inner current loop
    kpc = 0.74
    i_loop_range = [0.1, 0.01]

    # GFL parameters - outer loop 
    kp_range = [0.59, 2]
    ki_range = [7.36, 20]
    wz_range = [0.1*2*pi*60, 2*pi*60]

    # Inverter filter parameters 
    cf = [0.5, 1.5]*0.074
    #cf = [0.074]
    lf = [0.5, 1.5]*0.08
    #lf = [0.08]

    # System parameters 
    line_scale_range = [0.8, 1.0]
    inverter_share_range = [0.5, 0.8]
    gfm_share_range = [0.5, 1.0]
    load_scale_range = [0.8, 1.0]

    col_names = ["line_s", "load_s", "gfm_Ta", "gfm_kd", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "wz", "inv_share", "gfm_share"]

    df = DataFrame([name => [] for name in col_names])

    foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_Ta_range, gfm_kd_range, gfm_kq_range, kpv, kpv./v_loop_range, kpc, kpc./i_loop_range, kp_range, ki_range, cf, lf, wz_range, inverter_share_range, gfm_share_range))

    return df 
end

function generate_9bus_droop_parameter_dataframe()

    #GFM droop outer loop 
    gfm_kq_range = [0.01, 0.4] # Q-V droop 
    gfm_rp_range = [0.02, 0.05] # P-f droop 

    # GFM inner voltage loop 
    #kpv_range = collect([0.59, 1.0])
    kpv = 0.59 
    v_loop_range = collect([0.001, 0.1])
    #kiv_range = collect([400,736])

    # GFM/GFL inner current loop
    kpc = 0.74
    i_loop_range = [0.1, 0.01]

    # GFL parameters - outer loop 
    kp_range = [0.59, 2]
    ki_range = [7.36, 20]
    wz_range = [0.1*2*pi*60, 2*pi*60]

    # Inverter filter parameters 
    #cf = [0.5, 1.5]*0.074
    cf = [0.074]
    #lf = [0.5, 1.5]*0.08
    lf = [0.08]

    # System parameters 
    line_scale_range = [0.8, 1.0]
    inverter_share_range = [0.5, 0.8]
    gfm_share_range = [0.5, 1.0]
    load_scale_range = [0.8, 1.0]

    col_names = ["line_s", "load_s", "gfm_rp", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "wz", "inv_share", "gfm_share"]

    df = DataFrame([name => [] for name in col_names])

    foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_rp_range, gfm_kq_range, kpv, kpv./v_loop_range, kpc, kpc./i_loop_range, kp_range, ki_range, cf, lf, wz_range, inverter_share_range, gfm_share_range))

    return df 
end

# Plotting 
label_dict = Dict(
    "gfm_kiv"=> L"k_{iv}",
    "gfl_kpp"=> L"k_p",
    "gfl_kip"=> L"k_i",
    "load_s" => L"\mathrm{load\ scale}",
    "line_s" => L"\mathrm{line\ scale}",
    "gfm_Ta" => L"T_a",
    "gfm_kd" => L"k_d",
    "gfm_kq" => L"k_q",
    "gfm_kpv" => L"k_{pv}",
    "gfm_kpc" => L"k_{pc}",
    "gfm_kic" => L"k_{ic}",
    "cf" => L"c_f",
    "lf" => L"l_f",
    "wz" => L"\omega_z",
    "inv_share" => L"IBR\ share",
    "gfm_share" => L"GFM\ share",
    "gfm_rp" => L"R_p"
)


pf_label_dict = Dict(

)

### Simulation helper functions 
function generate_line_dict(sys)
    line_dict = Dict()
    
    if isempty(get_components(DynamicBranch,sys)) 
        for l in get_components(PSY.Line, sys)
            line_dict[l.name] = (l.r, l.x, l.b.from, l.b.to)
        end
        update_fcn =  static_line_update
    else
        for l in get_components(PSY.DynamicBranch, sys)
            line_dict[l.branch.name] = (l.branch.r, l.branch.x, l.branch.b.from, l.branch.b.to)
        end
        update_fcn = dynamic_line_update
    end
    return line_dict, update_fcn 
end

# Static line update 
function static_line_update(sys, line_dict, line_scale)
    for line in get_components(PSY.Line, sys)
        og_p = line_dict[line.name]
        line.r = og_p[1]*line_scale
        line.x = og_p[2]*line_scale 
        line.b = (from=og_p[3]*line_scale, to=og_p[4]*line_scale)
    end
end

# Dynamic line update 
function dynamic_line_update(sys, line_dict, line_scale)
    for line in get_components(PSY.DynamicBranch, sys)
        og_p = line_dict[line.branch.name]
        line.branch.r = og_p[1]*line_scale
        line.branch.x = og_p[2]*line_scale 
        line.branch.b = (from=og_p[3]*line_scale, to=og_p[4]*line_scale)
    end
end

# Load updates 
function generate_load_dict(sys)
    load_dict = Dict()

    for l in get_components(StandardLoad, sys)
        load_dict[l.name] = (l.impedance_active_power, l.impedance_reactive_power)
    end
    return load_dict 
end

function load_data(filename)
    x = DataFrame(CSV.File(filename))
    return x 
end 

# Analysis helper functions 

function get_stable(x)
    output = x[x.quadrant.=="green",:]
    #output = filter([:status, :status_d, :stable, :stable_d] => filter_stable, x)
    return output 
end

function get_unstable(x)
    output = x[x.quadrant.=="orange", :]
    #output = filter([:status, :status_d, :stable, :stable_d] => filter_unstable, x)
    return output 
end

function get_edge_blue(x)
    # where dynamic is unstable but algebraic is stable 
    output = x[x.quadrant.=="blue",:]
    #output = filter([:status, :status_d, :stable, :stable_d] => filter_edge, x)
    return output 
end 

function get_edge_red(x)
    # where algebraic is unstable and dynamic is stable 
    output = x[x.quadrant.=="red",:]
    #output = filter([:status, :status_d, :stable, :stable_d] => filter_edge, x)
    return output 
end 

function annotate_status(x)
    x[!,"quadrant"] .= ""
    for i in range(1,size(x)[1])
        st1 = x[i,:].status 
        st2 = x[i,:].status_d 
        eig1 = x[i,:].stable[1] # real part only 
        eig2 = x[i,:].stable_d[1] # real part only 
        isbuilt1 = st1 == "BUILT"
        isbuilt2 = st2 == "BUILT"
        a = (eig1 <= 0) && (eig2 > 0)  # algebraic stable, dynamic unstable  
        b = (eig2 <= 0) && (eig1 > 0) # algebraic unstable, dynamic stable  
        edge_condition_a = (isbuilt1 && isbuilt2) && a #algebraic stable, dynamic unstable  
        edge_condition_b = (isbuilt1 && isbuilt2) && b # algebraic unstable, dynamic stable  
        
        a = (eig1 > 0) && (eig2 > 0) 
        unstable_condition = (isbuilt1 && isbuilt2) && a
        a = (eig1 <= 0) && (eig2 <= 0) 
        stable_condition = (isbuilt1 && isbuilt2) && a

        if edge_condition_a
            x[i,"quadrant"] = "blue"
        elseif edge_condition_b
            x[i,"quadrant"] = "red"
        elseif unstable_condition
            x[i,"quadrant"] = "orange"
        elseif stable_condition
            x[i,"quadrant"] = "green"
        else
            x[i,"quadrant"] = "error"
        end 
    end 
end 

