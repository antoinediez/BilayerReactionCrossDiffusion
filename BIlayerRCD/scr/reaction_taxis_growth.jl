function reaction_B(u,v,c,param)
    #==== Linear KS (with inhibitor) ====#
    # r_u = param[:sc]*(c - param[:α]*v - u)
    # r_v = param[:sc]*(c - param[:β]*v)
    #====================================#

    #== Activator-Inhibitor (no cell) ==#
    # act_synth = max(0.0,min(param[:act_synth]+param[:act_autocat]*u - param[:inhib_act]*v,param[:max_act_synth]))
    # inhib_synth = max(0.0,min(param[:inhib_synth]+param[:inhib_autocat]*v + param[:inhib_act_crea]*u,param[:max_inhib_synth]))
    # r_u = param[:sc]*(act_synth - param[:act_deg]*u)
    # r_v = param[:sc]*(inhib_synth - param[:inhib_deg]*v)
    #====================================#

    #== Activator-Inhibitor (u,v->c) ==#
    # act_synth = max(0.0,min(param[:act_synth]+param[:act_autocat]*c - param[:inhib_act]*v,param[:max_act_synth]))
    # inhib_synth = max(0.0,min(param[:inhib_synth]+param[:inhib_autocat]*c + param[:inhib_act_crea]*c,param[:max_inhib_synth]))
    # r_u = param[:sc]*(act_synth - param[:act_deg]*u)
    # r_v = param[:sc]*(inhib_synth - param[:inhib_deg]*v)
    #====================================#

    #=========== Schnakenberg ===========#
    r_u = param[:sc] * (param[:a] - u + u^2 * v)
    r_v = param[:sc] * (param[:b] - u^2 * v)
    #====================================#

    #=========== No reaction ============#
    # r_u = 0.0
    # r_v = 0.0
    #====================================#
    return (r_u,r_v)
end

function reaction_S(u,v,param)
    # =========== Schnakenberg ===========#
    r_u = param[:sc] * (param[:a] - u + u^2 * v)
    r_v = param[:sc] * (param[:b] - u^2 * v)
    #====================================#

    #=== Linear Activator-Inhibitor =====#
    # act_synth = max(0.0,min(param[:act_synth]+param[:act_autocat]*u - param[:inhib_act]*v,param[:max_act_synth]))
    # inhib_synth = max(0.0,min(param[:inhib_synth]+param[:inhib_autocat]*v + param[:inhib_act_crea]*u,param[:max_inhib_synth]))
    # r_u = param[:sc]*(act_synth - param[:act_deg]*u)
    # r_v = param[:sc]*(inhib_synth - param[:inhib_deg]*v)
    #====================================#
    
    #=========== No reaction ============#
    # r_u = 0.0
    # r_v = 0.0
    #====================================#
    return (r_u,r_v)
end

function taxis(u,v,param)
    return u
    # return u/(1 + v)
    # return 4.0 * u^4/(1 + u^4)
    # return IfElse.ifelse(u<u0,0.0,(u-u0)^2/(u -u0 + 1.0))
    # return 0.0
end

function growth(c,param)
    return param[:sc] * c * (param[:c_eq] - c)
    # return 0.0
end

function taxis_h(c,param)
    return c
    # return c*exp(-param[:taxis_decay_rate]*c)
    # return 0.0
end