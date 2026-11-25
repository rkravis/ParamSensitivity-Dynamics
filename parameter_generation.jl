




# Parameter generation functions 
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

function vary_param(x)
    #return [0.5, 1.5] * x # +/- 50% 
    return [0.3, 1.7] * x # +/- 70%
end

function dont_vary_param(x)
    return [x]
end

function gfm_darco_gfl_pscad()
    """Generate parameter dataframe for GFM D'Arco and GFL PSCAD models. Also returns a function for updating the parameters of a model."""

    function make_parameter_dataframe()
        #GFM VSM outer loop - D'ARCO VSM PAPER 
        gfm_kq_range = vary_param(0.2)
        gfm_Ta_range = vary_param(2)
        gfm_kd_range = vary_param(400)
        gfm_kw = dont_vary_param(20) # frequency droop gain
        gfm_wf = dont_vary_param(1000) # filter frequency cutoff, reactive power loop (rad)

        # GFM inner voltage loop 
        gfm_kpv = vary_param(0.59)
        gfm_kiv = vary_param(736) 

        # GFM inner current loop
        gfm_kpc = vary_param(1.27)
        gfm_kic = vary_param(14.3)

        # GFM and GFL filter parameters 
        gfm_cf = dont_vary_param(0.074)
        gfm_lf = dont_vary_param(0.08)
        gfm_rf = dont_vary_param(0.003)

        # PLL parameters (GFM-D'Arco)
        # low pass filter = 500 rad/s 
        # kp = 0.084 
        # ki = 4.69 
        gfm_kp_pll = dont_vary_param(0.084)
        gfm_ki_pll = dont_vary_param(4.69)
        gfm_wl_pll = dont_vary_param(500)

        # GFL paper source - WALLACE PAPER ON OPEN SOURCE PSCAD 
        # PLL (paper uses reduced pll)
        gfl_kp_pll = dont_vary_param(0.084)
        gfl_ki_pll = dont_vary_param(4.69)
        gfl_wl_pll = dont_vary_param(1.32*2*pi*60)

        # Outer loop 
        gfl_kp_p_range = vary_param(2) 
        gfl_ki_p_range = vary_param(20) 

        gfl_kp_q_range = vary_param(2) 
        gfl_ki_q_range = vary_param(20)

        # GFL inner current loop 
        gfl_kpc = vary_param(0.37)
        gfl_kic = vary_param(0.7)

        # GFL filter parameters 
        gfl_cf = dont_vary_param(2.5)
        gfl_rf = dont_vary_param(0.016)
        gfl_lf = dont_vary_param(0.009)

        # GFL power measurement filter frequency
        # Geng & Hiskens paper uses 50 rad
        gfl_wz_range = dont_vary_param(0.132*2*pi*60)

        # System parameters - we choose 
        #line_scale_range = [0.8, 1.0]
        line_scale_range = [1.0]
        #inverter_share_range = [0.5, 0.8]
        inverter_share_range = [0.8]
        #gfm_share_range = [0.5, 1.0]
        gfm_share_range = [0.5]
        #load_scale_range = [0.8, 1.0]
        load_scale_range = [1.0]

        col_names = [
            "line_s", 
            "load_s", 
            "gfm_Ta", 
            "gfm_kd", 
            "gfm_kq",
            "gfm_kw",
            "gfm_wf", 
            "gfm_kpv", 
            "gfm_kiv", 
            "gfm_kpc", 
            "gfm_kic", 
            "gfl_kpp", 
            "gfl_kip",
            "gfl_kpq",
            "gfl_kiq",
            "gfl_cf", 
            "gfl_lf", 
            "gfl_rf",
            "gfm_cf",
            "gfm_lf",
            "gfm_rf",
            "inv_share",
            "gfm_share",
            "gfm_kp_pll",
            "gfm_ki_pll",
            "gfm_wl_pll",
            "gfl_kp_pll",
            "gfl_ki_pll",
            "gfl_wl_pll",
            "gfl_kpc",
            "gfl_kic",
            "gfl_wz"]

        df = DataFrame([name => [] for name in col_names])

        foreach(x -> push!(df, x), Iterators.product(
            line_scale_range, 
            load_scale_range, 
            gfm_Ta_range, 
            gfm_kd_range,
            gfm_kq_range,
            gfm_kw,
            gfm_wf,
            gfm_kpv,
            gfm_kiv,
            gfm_kpc,
            gfm_kic,
            gfl_kp_p_range,
            gfl_ki_p_range,
            gfl_kp_q_range,
            gfl_ki_q_range,
            gfl_cf,
            gfl_lf,
            gfl_rf, 
            gfm_cf,
            gfm_lf,
            gfm_rf,
            inverter_share_range,
            gfm_share_range,
            gfm_kp_pll,
            gfm_ki_pll,
            gfm_wl_pll,
            gfl_kp_pll,
            gfl_ki_pll,
            gfl_wl_pll,
            gfl_kpc,
            gfl_kic,
            gfl_wz_range))

        return df 
    end 

    function apply_param_update(
        sys,
        df_row,
        load_dict,
        line_dict,
        line_update_fcn,
        gfm_buses,
        gfl_buses,
        sm_buses)

        # Edits the sys object in place 
        # GFM parameters 
        gfm_kq = df_row.gfm_kq
        gfm_Ta = df_row.gfm_Ta
        gfm_kd = df_row.gfm_kd
        gfm_kw = df_row.gfm_kw
        gfm_wf = df_row.gfm_wf
        gfm_cf = df_row.gfm_cf
        gfm_lf = df_row.gfm_lf
        gfm_rf = df_row.gfm_rf
        gfm_kpv = df_row.gfm_kpv
        gfm_kiv = df_row.gfm_kiv
        gfm_kpc = df_row.gfm_kpc
        gfm_kic = df_row.gfm_kic
        gfm_kp_pll = df_row.gfm_kp_pll
        gfm_ki_pll = df_row.gfm_ki_pll
        gfm_wl_pll = df_row.gfm_wl_pll

        # GFL parameters 
        gfl_kpp = df_row.gfl_kpp
        gfl_kip = df_row.gfl_kip
        gfl_kpq = df_row.gfl_kpq
        gfl_kiq = df_row.gfl_kiq
        gfl_wz = df_row.gfl_wz
        gfl_kp_pll = df_row.gfl_kp_pll
        gfl_ki_pll = df_row.gfl_ki_pll
        gfl_wl_pll = df_row.gfl_wl_pll
        gfl_lf = df_row.gfl_lf
        gfl_cf = df_row.gfl_cf
        gfl_rf = df_row.gfl_rf
        gfl_kpc = df_row.gfl_kpc
        gfl_kic = df_row.gfl_kic

        # System parameters 
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

        num_GFM = length(gfm_buses)
        num_GFL = length(gfl_buses)
        num_SM = length(sm_buses)
        
        # Adjust injections for each generator 
        for g in get_components(Generator, sys)
            if g.bus.name in gfm_buses 
                g.active_power = gfm_eta*total_p/num_GFM/g.base_power # real power*scale / num of gen at that bus 
                g.reactive_power = gfm_eta*total_q/num_GFM/g.base_power;
                # Virtual inertia 
                g.dynamic_injector.outer_control.reactive_power_control.kq = gfm_kq;
                g.dynamic_injector.inner_control.kpc = gfm_kpc;
                g.dynamic_injector.inner_control.kpv = gfm_kpv;
                g.dynamic_injector.inner_control.kic = gfm_kic;
                g.dynamic_injector.inner_control.kiv = gfm_kiv;
                g.dynamic_injector.filter.lf = gfm_lf 
                g.dynamic_injector.filter.cf = gfm_cf
                g.dynamic_injector.filter.rf = gfm_rf
                g.dynamic_injector.outer_control.active_power_control.kd = gfm_kd;
                g.dynamic_injector.outer_control.active_power_control.Ta = gfm_Ta;
                g.dynamic_injector.outer_control.active_power_control.kω = gfm_kw;
                g.dynamic_injector.outer_control.reactive_power_control.ωf = gfm_wf;
                g.dynamic_injector.freq_estimator.kp_pll = gfm_kp_pll
                g.dynamic_injector.freq_estimator.ki_pll = gfm_ki_pll
                g.dynamic_injector.freq_estimator.ω_lp = gfm_wl_pll
            elseif g.bus.name in gfl_buses 
                g.active_power = gfl_eta*total_p/num_GFL/g.base_power; 
                g.reactive_power = gfl_eta*total_q/num_GFL/g.base_power;
                g.dynamic_injector.outer_control.active_power_control.Ki_p = gfl_kip; 
                g.dynamic_injector.outer_control.active_power_control.Kp_p = gfl_kpp;
                g.dynamic_injector.outer_control.reactive_power_control.Ki_q = gfl_kiq; 
                g.dynamic_injector.outer_control.reactive_power_control.Kp_q = gfl_kpq;
                g.dynamic_injector.filter.lf = gfl_lf 
                g.dynamic_injector.filter.cf = gfl_cf 
                g.dynamic_injector.filter.rf = gfl_rf
                g.dynamic_injector.inner_control.kpc = gfl_kpc 
                g.dynamic_injector.inner_control.kic = gfl_kic 
                g.dynamic_injector.outer_control.active_power_control.ωz = gfl_wz
                g.dynamic_injector.outer_control.reactive_power_control.ωf = gfl_wz
                g.dynamic_injector.freq_estimator.kp_pll = gfl_kp_pll 
                g.dynamic_injector.freq_estimator.ki_pll = gfl_ki_pll
                g.dynamic_injector.freq_estimator.ω_lp = gfl_wl_pll
            elseif g.bus.name in sm_buses
                g.active_power = sm_eta*total_p/num_SM/g.base_power;
                g.reactive_power = sm_eta*total_q/num_SM/g.base_power;
            end
        end
    end

    df = make_parameter_dataframe()
    return apply_param_update, df 
end

function gfm_markovic_gfl_pscad()
    """Generate parameter dataframe for GFM Markovic (Droop) and GFL PSCAD models. Also returns a function for updating the parameters of a model."""

    function make_parameter_dataframe()
        # GFM outer loop - Markovic (Droop) paper
        gfm_kq_range = vary_param(0.02) #2%
        gfm_kp_range = vary_param(0.001) #0.1%

        # GFM inner voltage loop 
        # voltage control has feedforward 
        gfm_kpv = vary_param(0.52)
        gfm_kiv = vary_param(1.16) 
        gfm_kff_v = dont_vary_param(1) #feedforward

        # GFM inner current loop
        gfm_kpc = vary_param(0.74)
        gfm_kic = vary_param(1.19)
        gfm_kff_i = dont_vary_param(0) # feedforward

        # GFM filter parameters 
        gfm_cf = dont_vary_param(0.074)
        gfm_lf = dont_vary_param(0.08)
        gfm_rf = dont_vary_param(0.03) # 0.03

        # # PLL parameters (GFM)
        gfm_kp_pll = dont_vary_param(0.084)
        gfm_ki_pll = dont_vary_param(4.69)

        # GFM virtual impedance parameters
        gfm_rv = dont_vary_param(0.001)
        gfm_xv = dont_vary_param(0.2)

        # GFL paper source - WALLACE PAPER ON OPEN SOURCE PSCAD 
        # PLL (paper uses reduced pll)
        gfl_kp_pll = dont_vary_param(0.084)
        gfl_ki_pll = dont_vary_param(4.69)
        gfl_wl_pll = dont_vary_param(1.32*2*pi*60)

        # Outer loop 
        gfl_kp_p_range = vary_param(2) 
        gfl_ki_p_range = vary_param(20) 

        gfl_kp_q_range = vary_param(2) 
        gfl_ki_q_range = vary_param(20)

        # GFL inner current loop 
        gfl_kpc = vary_param(0.37)
        gfl_kic = vary_param(0.7)

        # GFL filter parameters 
        gfl_cf = dont_vary_param(2.5)
        gfl_rf = dont_vary_param(0.016)
        gfl_lf = dont_vary_param(0.009)

        # GFL power measurement filter frequency
        # Geng & Hiskens paper uses 50 rad
        gfl_wz_range = dont_vary_param(0.132*2*pi*60)

        # System parameters - we choose 
        #line_scale_range = [0.8, 1.0]
        line_scale_range = [1.0]
        #inverter_share_range = [0.5, 0.8]
        inverter_share_range = [0.8]
        #gfm_share_range = [0.5, 1.0]
        gfm_share_range = [0.5]
        #load_scale_range = [0.8, 1.0]
        load_scale_range = [1.0]

        col_names = [
            "line_s", 
            "load_s",  
            "gfm_kp", 
            "gfm_kq",
            "gfm_kpv", 
            "gfm_kiv",
            "gfm_kff_v", 
            "gfm_kpc", 
            "gfm_kic",
            "gfm_kff_i",
            "gfm_rv",
            "gfm_xv",
            "gfl_kpp", 
            "gfl_kip",
            "gfl_kpq",
            "gfl_kiq",
            "gfl_cf", 
            "gfl_lf", 
            "gfl_rf",
            "gfm_cf",
            "gfm_lf",
            "gfm_rf",
            "inv_share",
            "gfm_share",
            "gfm_kp_pll",
            "gfm_ki_pll",
            "gfl_kp_pll",
            "gfl_ki_pll",
            "gfl_wl_pll",
            "gfl_kpc",
            "gfl_kic",
            "gfl_wz"]

        df = DataFrame([name => [] for name in col_names])

        foreach(x -> push!(df, x), Iterators.product(
            line_scale_range, 
            load_scale_range, 
            gfm_kp_range,
            gfm_kq_range,
            gfm_kpv,
            gfm_kiv,
            gfm_kff_v,
            gfm_kpc,
            gfm_kic,
            gfm_kff_i,
            gfm_rv,
            gfm_xv,
            gfl_kp_p_range,
            gfl_ki_p_range,
            gfl_kp_q_range,
            gfl_ki_q_range,
            gfl_cf,
            gfl_lf,
            gfl_rf, 
            gfm_cf,
            gfm_lf,
            gfm_rf,
            inverter_share_range,
            gfm_share_range,
            gfm_kp_pll,
            gfm_ki_pll,
            gfl_kp_pll,
            gfl_ki_pll,
            gfl_wl_pll,
            gfl_kpc,
            gfl_kic,
            gfl_wz_range))

        return df 
    end 

    function apply_param_update(
        sys,
        df_row,
        load_dict,
        line_dict,
        line_update_fcn,
        gfm_buses,
        gfl_buses,
        sm_buses)

        # Edits the sys object in place 
        # GFM parameters 
        gfm_kq = df_row.gfm_kq
        gfm_kp = df_row.gfm_kp
        gfm_cf = df_row.gfm_cf
        gfm_lf = df_row.gfm_lf
        gfm_rf = df_row.gfm_rf
        gfm_kpv = df_row.gfm_kpv
        gfm_kiv = df_row.gfm_kiv
        gfm_kffv = df_row.gfm_kff_v
        gfm_kpc = df_row.gfm_kpc
        gfm_kic = df_row.gfm_kic
        gfm_kffi = df_row.gfm_kff_i
        gfm_rv = df_row.gfm_rv
        gfm_xv = df_row.gfm_xv
        gfm_kp_pll = df_row.gfm_kp_pll
        gfm_ki_pll = df_row.gfm_ki_pll

        # GFL parameters 
        gfl_kpp = df_row.gfl_kpp
        gfl_kip = df_row.gfl_kip
        gfl_kpq = df_row.gfl_kpq
        gfl_kiq = df_row.gfl_kiq
        gfl_wz = df_row.gfl_wz
        gfl_kp_pll = df_row.gfl_kp_pll
        gfl_ki_pll = df_row.gfl_ki_pll
        gfl_wl_pll = df_row.gfl_wl_pll
        gfl_lf = df_row.gfl_lf
        gfl_cf = df_row.gfl_cf
        gfl_rf = df_row.gfl_rf
        gfl_kpc = df_row.gfl_kpc
        gfl_kic = df_row.gfl_kic

        # System parameters 
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

        num_GFM = length(gfm_buses)
        num_GFL = length(gfl_buses)
        num_SM = length(sm_buses)
        
        # Adjust injections for each generator 
        for g in get_components(Generator, sys)
            if g.bus.name in gfm_buses 
                g.active_power = gfm_eta*total_p/num_GFM/g.base_power # real power*scale / num of gen at that bus 
                g.reactive_power = gfm_eta*total_q/num_GFM/g.base_power;
                g.dynamic_injector.outer_control.reactive_power_control.kq = gfm_kq;
                g.dynamic_injector.outer_control.active_power_control.Rp = gfm_kp;
                g.dynamic_injector.inner_control.kpc = gfm_kpc;
                g.dynamic_injector.inner_control.kpv = gfm_kpv;
                g.dynamic_injector.inner_control.kic = gfm_kic;
                g.dynamic_injector.inner_control.kiv = gfm_kiv;
                g.dynamic_injector.filter.lf = gfm_lf 
                g.dynamic_injector.filter.cf = gfm_cf
                g.dynamic_injector.filter.rf = gfm_rf
                g.dynamic_injector.freq_estimator.kp_pll = gfm_kp_pll
                g.dynamic_injector.freq_estimator.ki_pll = gfm_ki_pll
                g.dynamic_injector.inner_control.lv = gfm_xv
                g.dynamic_injector.inner_control.rv = gfm_rv 
                g.dynamic_injector.inner_control.kffv = gfm_kffi
                g.dynamic_injector.inner_control.kffi = gfm_kffv 

            elseif g.bus.name in gfl_buses 
                g.active_power = gfl_eta*total_p/num_GFL/g.base_power; 
                g.reactive_power = gfl_eta*total_q/num_GFL/g.base_power;
                g.dynamic_injector.outer_control.active_power_control.Ki_p = gfl_kip; 
                g.dynamic_injector.outer_control.active_power_control.Kp_p = gfl_kpp;
                g.dynamic_injector.outer_control.reactive_power_control.Ki_q = gfl_kiq; 
                g.dynamic_injector.outer_control.reactive_power_control.Kp_q = gfl_kpq;
                g.dynamic_injector.filter.lf = gfl_lf 
                g.dynamic_injector.filter.cf = gfl_cf 
                g.dynamic_injector.filter.rf = gfl_rf
                g.dynamic_injector.inner_control.kpc = gfl_kpc 
                g.dynamic_injector.inner_control.kic = gfl_kic 
                g.dynamic_injector.outer_control.active_power_control.ωz = gfl_wz
                g.dynamic_injector.outer_control.reactive_power_control.ωf = gfl_wz
                g.dynamic_injector.freq_estimator.kp_pll = gfl_kp_pll 
                g.dynamic_injector.freq_estimator.ki_pll = gfl_ki_pll
                g.dynamic_injector.freq_estimator.ω_lp = gfl_wl_pll
            elseif g.bus.name in sm_buses
                g.active_power = sm_eta*total_p/num_SM/g.base_power;
                g.reactive_power = sm_eta*total_q/num_SM/g.base_power;
            end
        end
    end

    df = make_parameter_dataframe()
    return apply_param_update, df 
end

function gfm_darco_gfl_bri()
    """Generate parameter dataframe for GFM D'Arco and GFL Bri models. Also returns a function for updating the parameters of a model."""

    function make_parameter_dataframe()
        #GFM VSM outer loop - D'ARCO VSM PAPER 
        gfm_kq_range = vary_param(0.2)
        gfm_Ta_range = vary_param(2)
        gfm_kd_range = vary_param(400)
        gfm_kw = dont_vary_param(20) # frequency droop gain
        gfm_wf = dont_vary_param(1000) # filter frequency cutoff, reactive power loop (rad)

        # GFM inner voltage loop 
        gfm_kpv = vary_param(0.59)
        gfm_kiv = vary_param(736) 

        # GFM inner current loop
        gfm_kpc = vary_param(1.27)
        gfm_kic = vary_param(14.3)

        # GFM and GFL filter parameters 
        gfm_cf = dont_vary_param(0.074)
        gfm_lf = dont_vary_param(0.08)
        gfm_rf = dont_vary_param(0.003)

        # PLL parameters (GFM-D'Arco)
        # low pass filter = 500 rad/s 
        # kp = 0.084 
        # ki = 4.69 
        gfm_kp_pll = dont_vary_param(0.084)
        gfm_ki_pll = dont_vary_param(4.69)
        gfm_wl_pll = dont_vary_param(500)

        # GFL paper source - Bri/Amir paper 
        # TODO what kind of PLL is used in the Bri/Amir paper?
        # PLL
        gfl_kp_pll = dont_vary_param(52)
        gfl_ki_pll = dont_vary_param(415)
        gfl_wl_pll = dont_vary_param(500)

        # Outer loop power control 
        gfl_kp_p_range = vary_param(2.1) 
        gfl_ki_p_range = vary_param(21) 

        gfl_kp_q_range = vary_param(2.08) 
        gfl_ki_q_range = vary_param(20.8)

        # GFL inner current loop 
        gfl_kpc = vary_param(0.38)
        gfl_kic = vary_param(0.7)
        gfl_kff = dont_vary_param(1)

        # GFL filter parameters 
        gfl_cf = dont_vary_param(2.55)
        gfl_rf = dont_vary_param(0.016)
        gfl_lf = dont_vary_param(0.009)

        # coupling parameters 
        # l = 0.0017
        # r = 0.003 

        # GFL power measurement filter frequency
        # Geng & Hiskens paper uses 50 rad
        gfl_wz_range = dont_vary_param(50*2*pi*60)

        # System parameters - we choose 
        #line_scale_range = [0.8, 1.0]
        line_scale_range = [1.0]
        #inverter_share_range = [0.5, 0.8]
        inverter_share_range = [0.8]
        #gfm_share_range = [0.5, 1.0]
        gfm_share_range = [0.5]
        #load_scale_range = [0.8, 1.0]
        load_scale_range = [1.0]

        col_names = [
            "line_s", 
            "load_s", 
            "gfm_Ta", 
            "gfm_kd", 
            "gfm_kq",
            "gfm_kw",
            "gfm_wf", 
            "gfm_kpv", 
            "gfm_kiv", 
            "gfm_kpc", 
            "gfm_kic", 
            "gfl_kpp", 
            "gfl_kip",
            "gfl_kpq",
            "gfl_kiq",
            "gfl_cf", 
            "gfl_lf", 
            "gfl_rf",
            "gfm_cf",
            "gfm_lf",
            "gfm_rf",
            "inv_share",
            "gfm_share",
            "gfm_kp_pll",
            "gfm_ki_pll",
            "gfm_wl_pll",
            "gfl_kp_pll",
            "gfl_ki_pll",
            "gfl_wl_pll",
            "gfl_kpc",
            "gfl_kic",
            "gfl_wz"]

        df = DataFrame([name => [] for name in col_names])

        foreach(x -> push!(df, x), Iterators.product(
            line_scale_range, 
            load_scale_range, 
            gfm_Ta_range, 
            gfm_kd_range,
            gfm_kq_range,
            gfm_kw,
            gfm_wf,
            gfm_kpv,
            gfm_kiv,
            gfm_kpc,
            gfm_kic,
            gfl_kp_p_range,
            gfl_ki_p_range,
            gfl_kp_q_range,
            gfl_ki_q_range,
            gfl_cf,
            gfl_lf,
            gfl_rf, 
            gfm_cf,
            gfm_lf,
            gfm_rf,
            inverter_share_range,
            gfm_share_range,
            gfm_kp_pll,
            gfm_ki_pll,
            gfm_wl_pll,
            gfl_kp_pll,
            gfl_ki_pll,
            gfl_wl_pll,
            gfl_kpc,
            gfl_kic,
            gfl_wz_range))

        return df 
    end 

    function apply_param_update(
        sys,
        df_row,
        load_dict,
        line_dict,
        line_update_fcn,
        gfm_buses,
        gfl_buses,
        sm_buses)

        # Edits the sys object in place 
        # GFM parameters 
        gfm_kq = df_row.gfm_kq
        gfm_Ta = df_row.gfm_Ta
        gfm_kd = df_row.gfm_kd
        gfm_kw = df_row.gfm_kw
        gfm_wf = df_row.gfm_wf
        gfm_cf = df_row.gfm_cf
        gfm_lf = df_row.gfm_lf
        gfm_rf = df_row.gfm_rf
        gfm_kpv = df_row.gfm_kpv
        gfm_kiv = df_row.gfm_kiv
        gfm_kpc = df_row.gfm_kpc
        gfm_kic = df_row.gfm_kic
        gfm_kp_pll = df_row.gfm_kp_pll
        gfm_ki_pll = df_row.gfm_ki_pll
        gfm_wl_pll = df_row.gfm_wl_pll

        # GFL parameters 
        gfl_kpp = df_row.gfl_kpp
        gfl_kip = df_row.gfl_kip
        gfl_kpq = df_row.gfl_kpq
        gfl_kiq = df_row.gfl_kiq
        gfl_wz = df_row.gfl_wz
        gfl_kp_pll = df_row.gfl_kp_pll
        gfl_ki_pll = df_row.gfl_ki_pll
        gfl_wl_pll = df_row.gfl_wl_pll
        gfl_lf = df_row.gfl_lf
        gfl_cf = df_row.gfl_cf
        gfl_rf = df_row.gfl_rf
        gfl_kpc = df_row.gfl_kpc
        gfl_kic = df_row.gfl_kic

        # System parameters 
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

        num_GFM = length(gfm_buses)
        num_GFL = length(gfl_buses)
        num_SM = length(sm_buses)
        
        # Adjust injections for each generator 
        for g in get_components(Generator, sys)
            if g.bus.name in gfm_buses 
                g.active_power = gfm_eta*total_p/num_GFM/g.base_power # real power*scale / num of gen at that bus 
                g.reactive_power = gfm_eta*total_q/num_GFM/g.base_power;
                # Virtual inertia 
                g.dynamic_injector.outer_control.reactive_power_control.kq = gfm_kq;
                g.dynamic_injector.inner_control.kpc = gfm_kpc;
                g.dynamic_injector.inner_control.kpv = gfm_kpv;
                g.dynamic_injector.inner_control.kic = gfm_kic;
                g.dynamic_injector.inner_control.kiv = gfm_kiv;
                g.dynamic_injector.filter.lf = gfm_lf 
                g.dynamic_injector.filter.cf = gfm_cf
                g.dynamic_injector.filter.rf = gfm_rf
                g.dynamic_injector.outer_control.active_power_control.kd = gfm_kd;
                g.dynamic_injector.outer_control.active_power_control.Ta = gfm_Ta;
                g.dynamic_injector.outer_control.active_power_control.kω = gfm_kw;
                g.dynamic_injector.outer_control.reactive_power_control.ωf = gfm_wf;
                g.dynamic_injector.freq_estimator.kp_pll = gfm_kp_pll
                g.dynamic_injector.freq_estimator.ki_pll = gfm_ki_pll
                g.dynamic_injector.freq_estimator.ω_lp = gfm_wl_pll
            elseif g.bus.name in gfl_buses 
                g.active_power = gfl_eta*total_p/num_GFL/g.base_power; 
                g.reactive_power = gfl_eta*total_q/num_GFL/g.base_power;
                g.dynamic_injector.outer_control.active_power_control.Ki_p = gfl_kip; 
                g.dynamic_injector.outer_control.active_power_control.Kp_p = gfl_kpp;
                g.dynamic_injector.outer_control.reactive_power_control.Ki_q = gfl_kiq; 
                g.dynamic_injector.outer_control.reactive_power_control.Kp_q = gfl_kpq;
                g.dynamic_injector.filter.lf = gfl_lf 
                g.dynamic_injector.filter.cf = gfl_cf 
                g.dynamic_injector.filter.rf = gfl_rf
                g.dynamic_injector.inner_control.kpc = gfl_kpc 
                g.dynamic_injector.inner_control.kic = gfl_kic 
                g.dynamic_injector.outer_control.active_power_control.ωz = gfl_wz
                g.dynamic_injector.outer_control.reactive_power_control.ωf = gfl_wz
                g.dynamic_injector.freq_estimator.kp_pll = gfl_kp_pll 
                g.dynamic_injector.freq_estimator.ki_pll = gfl_ki_pll
                g.dynamic_injector.freq_estimator.ω_lp = gfl_wl_pll
            elseif g.bus.name in sm_buses
                g.active_power = sm_eta*total_p/num_SM/g.base_power;
                g.reactive_power = sm_eta*total_q/num_SM/g.base_power;
            end
        end
    end

    df = make_parameter_dataframe()
    return apply_param_update, df 
end





# function gfm_markovic_gfl_pscad()
#     " +/- 50% of the nominal values from papers listed"
#     min_s = 0.5
#     max_s = 1.5
#     #GFM Droop outer loop - Understanding small signal stability 

#     gfm_kp_range = [min_s,max_s]*0.2
#     gfm_kq_range = [min_s,max_s]*0.001 

#     # GFM inner voltage loop 
#     kpv = [min_s,max_s]*0.52
#     kiv = [min_s,max_s]*1.16

#     # GFM inner current loop
#     kpc = [min_s,max_s]*0.74
#     kic = [min_s,max_s]*1.19

#     # GFM and GFL filter parameters 
#     cf = [0.074]
#     lf = [0.08]

#     # GFL paper source - WALLACE PAPER ON OPEN SOURCE PSCAD 
#     # GFL outer loop 
#     kp_range = [min_s,max_s]*2 #[0.59, 2]
#     ki_range = [min_s,max_s]*20 #[7.36, 20]
#     # kp_range = [2]
#     # ki_range = [20]

#     # GFL PLL - using a reduced PLL 
#     kp_pll = [min_s,max_s]*2 
#     ki_pll = [min_s,max_s]*20 

#     # GFL inner current loop 
#     kpc_gfl = [min_s, max_s]*0.37
#     kic_gfl = [min_s, max_s]*0.7

#     # GFL power measurement filter frequency
#     # No reference for this parameter, Geng & Hiskens paper uses 50 rad
#     #wz_range = [min_s,max_s]*50 #[0.1*2*pi*60, 2*pi*60]
#     wz_range = [50]

#     # System parameters - we choose 
#     #line_scale_range = [0.8, 1.0]
#     line_scale_range = [1.0]
#     inverter_share_range = [0.5, 0.8]
#     gfm_share_range = [0.5, 1.0]
#     load_scale_range = [1.0]
#     #load_scale_range = [0.8, 1.0]

#     col_names = ["line_s", "load_s", "gfm_kp", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "wz", "inv_share", "gfm_share", "kp_pll", "ki_pll", "kpc_gfl", "kic_gfl"]
    

#     df = DataFrame([name => [] for name in col_names])

#     foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_kp_range, gfm_kq_range, kpv, kiv, kpc, kic, kp_range, ki_range, cf, lf, wz_range, inverter_share_range, gfm_share_range, kp_pll, ki_pll, kpc_gfl, kic_gfl))

#     return df  
# end

# function gfm_darco_gfl_green()

#     # GFM parameters D'Arco paper
#     gfm_kq_range = vary_param(0.2)
#     gfm_Ta_range = vary_param(2)
#     gfm_kd_range = vary_param(400)

#     # Inner voltage loop 
#     kpv = vary_param(0.59)
#     kiv = vary_param(736) #collect([0.001, 0.1])

#     # Inner current loop
#     kpc = vary_param(1.27)
#     kic = vary_param(14.3) #[0.1, 0.01]

#     # Filter parameters 
#     lf_gfm = [0.074]
#     cf_gfm = [0.08]

#     # GFL parameters: Tim Green MOR paper 
    
#     # Outer loop 
#     kp_range = vary_param(6.84)
#     ki_range = vary_param(1.37)

#     kp_q_range = [0.0] #vary_param(0)
#     ki_q_range = vary_param(0.4)

#     # GFL PLL - using a reduced PLL 
#     kp_pll = vary_param(0.2)
#     ki_pll = vary_param(0.01)

#     # GFL inner current loop 
#     kpc_gfl = vary_param(1.19)
#     kic_gfl = vary_param(7.17)

#     # Filter params 
#     lf = [0.0498]
#     cf = [0.12]

#     # GFL power measurement filter frequency
#     # No reference for this parameter, Geng & Hiskens paper uses 50 rad
#     #wz_range = [min_s,max_s]*50 #[0.1*2*pi*60, 2*pi*60]

#     wz_range = [50]

#     # System parameters - we choose these  
#     #line_scale_range = [0.8, 1.0]
#     line_scale_range = [1.0]
#     #inverter_share_range = [0.5, 0.8]
#     inverter_share_range = [0.5]
#     #gfm_share_range = [0.5, 1.0]
#     gfm_share_range = [0.5]
#     #load_scale_range = [0.8, 1.0]
#     load_scale_range = [1.0]

#     col_names = ["line_s", "load_s", "gfm_Ta", "gfm_kd", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "cf_gfm", "lf_gfm", "wz", "inv_share", "gfm_share", "kp_pll", "ki_pll", "kpc_gfl", "kic_gfl", "gfl_kpq", "gfl_kiq"]

#     df = DataFrame([name => [] for name in col_names])

#     foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_Ta_range, gfm_kd_range, gfm_kq_range, kpv, kiv, kpc, kic, kp_range, ki_range, cf, lf, cf_gfm, lf_gfm, wz_range, inverter_share_range, gfm_share_range, kp_pll, ki_pll, kpc_gfl, kic_gfl, kp_q_range, ki_q_range))

#     return df  
# end


# function gfm_markovic_gfl_green()
#     " +/- 50% of the nominal values from papers listed"
#     min_s = 0.5
#     max_s = 1.5
#     #GFM Droop outer loop - Understanding small signal stability 

#     gfm_kp_range = [min_s,max_s]*0.2
#     gfm_kq_range = [min_s,max_s]*0.001 

#     # GFM inner voltage loop 
#     kpv = [min_s,max_s]*0.52
#     kiv = [min_s,max_s]*1.16

#     # GFM inner current loop
#     kpc = [min_s,max_s]*0.74
#     kic = [min_s,max_s]*1.19

#     # GFM and GFL filter parameters 
#     cf = [0.074]
#     lf = [0.08]

#     # GFL paper source - Tim Green MOR paper 
#     # GFL outer loop 
#     kp_range = [min_s,max_s]*6.84 #[0.59, 2]
#     ki_range = [min_s,max_s]*1.37 #[7.36, 20]
#     # kp_range = [2]
#     # ki_range = [20]

#     # GFL PLL - using a reduced PLL 
#     kp_pll = [min_s,max_s]*0.2
#     ki_pll = [min_s,max_s]*0.01

#     # GFL inner current loop 
#     kpc_gfl = [min_s, max_s]*1.19
#     kic_gfl = [min_s, max_s]*7.17

#     # GFL power measurement filter frequency
#     # No reference for this parameter, Geng & Hiskens paper uses 50 rad
#     wz_range = [min_s,max_s]*50 #[0.1*2*pi*60, 2*pi*60]
#     wz_range = [50]

#     # System parameters - we choose 
#     #line_scale_range = [0.8, 1.0]
#     line_scale_range = [1.0]
#     inverter_share_range = [0.5, 0.8]
#     gfm_share_range = [0.5, 1.0]
#     load_scale_range = [1.0]
#     #load_scale_range = [0.8, 1.0]
#     col_names = ["line_s", "load_s", "gfm_kp", "gfm_kq", "gfm_kpv", "gfm_kiv", "gfm_kpc", "gfm_kic", "gfl_kpp", "gfl_kip", "cf", "lf", "wz", "inv_share", "gfm_share", "kp_pll", "ki_pll", "kpc_gfl", "kic_gfl"]
    

#     df = DataFrame([name => [] for name in col_names])

#     foreach(x -> push!(df, x), Iterators.product(line_scale_range, load_scale_range, gfm_kp_range, gfm_kq_range, kpv, kiv, kpc, kic, kp_range, ki_range, cf, lf, wz_range, inverter_share_range, gfm_share_range, kp_pll, ki_pll, kpc_gfl, kic_gfl))

#     return df  
# end

