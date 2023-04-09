using NLsolve
using ProgressMeter

function equilibria(reac_fun,param;U0,kwargs...)
	
    function f!(F,U)
        F .= reac_fun(U,param;kwargs...)
	end

    equi = nlsolve(f!,U0)
    attempt = 1
    while ((~equi.f_converged) || sum(equi.zero .< 0)>0) && (attempt<1e4)
        equi = nlsolve(f!,10.0.*rand(length(U0)).*U0)
        attempt += 1
    end
    if equi.f_converged
        return equi.zero
    else
        println("Cannot find an equilibrium...")
        return nothing
    end
end


function equilibria_cpl(param;A,B,U0,kwargs...)
    n = size(A,1)
    m = length(U0) - n
    Δ = @MVector zeros(n)
    function f!(F,U;Δ=Δ)
        var_S = @view U[1:n]
        var_B = @view U[(n+1):(n+m)]
        r_S = reaction_S(var_S,param)
        r_B = reaction_B(var_B,param)
        
        for i in 1:n
            Δ[i] = U[n+i] - U[i]
        end
        cpl_S = A*Δ
        cpl_B = B*Δ
        for i in 1:n
            F[i] = r_S[i] + param[:η]*cpl_S[i]
            F[n+i] = r_B[i] - param[:η]*cpl_B[i]
        end
        for i in (n+1):m
            F[n+i] = r_B[i]
        end
	end

    equi = nlsolve(f!,U0)
    attempt = 1
    while ((~equi.f_converged) || sum(equi.zero .< 0)>0) && (attempt<1e4)
        equi = nlsolve(f!,10.0.*rand(length(U0)).*U0)
        attempt += 1
    end
    if equi.f_converged
        return equi.zero
    else
        println("Cannot find an equilibrium...")
        return nothing
    end
end


function critical(reac_fun,diff_fun,param,η_sb,δ_sb,var;η_val=0.0,U0=10*rand(2),kwargs...)

    @variables var_symb[1:length(var)]

    X,η,δ,coeffs,_ = compute_coeffs(
            reac_fun,diff_fun,
            var_symb,param,
            η_sb,δ_sb;
            kwargs...
    )

    ∂X = Differential(X)
    a0 = substitute(coeffs[1],Dict(η=>η_val,var_symb=>var))
    dXa0 = expand_derivatives(∂X(a0))

    a0_fun2 = eval(build_function(a0,X,δ;expression=Val{false}))
    dXa0_fun2 = eval(build_function(dXa0,X,δ;expression=Val{false}))
    
    function f!(F,U;X=X,a0=a0,dXa0=dXa0)
        # U = X,δ
        F[1] = a0_fun2(U[1],U[2])
        F[2] = dXa0_fun2(U[1],U[2])
    end

    crit = nlsolve(f!,U0)
    attempt = 1
    while ((~crit.f_converged) || (crit.zero[1]<0) || (crit.zero[2]<0)) && (attempt<1e4)
        crit = nlsolve(f!,10.0.*rand(length(U0)).*U0)
        attempt += 1
    end
    if crit.f_converged
        return crit.zero
    else
        println("Cannot find the critical values...")
        return nothing
    end
end


function critical_curve(η_all,reac_fun,diff_fun,param,η_sb,δ_sb,equi_fun;equi0,U0=10.0*rand(2),kwargs...)
    X_cr = Float64[]
    δ_cr = Float64[]
    η_cr = Float64[]

    @variables var[1:length(equi0)]

    X,η,δ,coeffs,_ = compute_coeffs(
            reac_fun,diff_fun,
            var,param,
            η_sb,δ_sb;
            kwargs...
    )

    ∂X = Differential(X)
    a0 = coeffs[1]
    dXa0 = expand_derivatives(∂X(a0))

    a0_fun2 = eval(build_function(a0,X,var,η,δ;expression=Val{false}))
    dXa0_fun2 = eval(build_function(dXa0,X,var,η,δ;expression=Val{false}))
    
    @showprogress for η in η_all
        param[η_sb] = η
        equi = equi_fun(param;U0=equi0,kwargs...)
        function f!(F,U;X=X,a0=a0,dXa0=dXa0)
            # U = X,δ
            F[1] = a0_fun2(U[1],equi,η,U[2])
            F[2] = dXa0_fun2(U[1],equi,η,U[2])
        end
        if length(X_cr)>0
            U0=[X_cr[end],δ_cr[end]]
        end
        crit = nlsolve(f!,U0)
        attempt = 1
        while ((~crit.f_converged) || (crit.zero[1]<0)) && (attempt<1e4)
            crit = nlsolve(f!,10.0.*rand(length(U0)).*U0)
            attempt += 1
        end
        if crit.f_converged
            append!(X_cr,crit.zero[1])
            append!(δ_cr,crit.zero[2])
            append!(η_cr,η)
        end
    end
    return η_cr,δ_cr,X_cr
end







