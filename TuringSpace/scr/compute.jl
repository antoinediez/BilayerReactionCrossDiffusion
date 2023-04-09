using Symbolics
using ProgressMeter

function compute_jac(reac_fun,var,param;kwargs...)
    return Symbolics.jacobian(collect(reac_fun(var,param;kwargs...)),var)
end

function compute_coeffs(reaction,diffusion,var_symb,param,η_sb,δ_sb;kwargs...)
    @variables λ, X
    η = parse_expr_to_symbolic(η_sb,Main)
    δ = parse_expr_to_symbolic(δ_sb,Main)
    param_symb = convert(Dict{Symbol,Real},param)
    param_symb[η_sb] = η
    param_symb[δ_sb] = δ
    J = compute_jac(reaction,var_symb,param_symb;kwargs...)
    D = diffusion(var_symb,param_symb;kwargs...)
    Δ = dispersion(λ,X,J,D)
    n = size(J,1)
    coeffs = λ_coeff(Δ,λ,n-1)
    H = η_coeff(Δ,λ,η,n-1)
    return X,η,δ,coeffs,H
end


function compute_slope(param,η_sb,δ_sb,equi0,X_cr,δ_cr;kwargs...)

    A = kwargs[:A]
    B = kwargs[:B]
    
    n = size(A,1)
    m = length(equi0) - n 
    equi_S = @view equi0[1:n]
    equi_B = @view equi0[(n+1):(n+m)]

    @variables uS[1:n], uB[1:m]

    X,_,δ,_,H = compute_coeffs(
            reaction_cpl,diffusion_cpl,
            [uS;uB],param,
            η_sb,δ_sb;
            kwargs...
    )

    JS = compute_jac(reaction_S,uS,param)
    JB = compute_jac(reaction_B,uB,param)

    JS_val = Symbolics.value.(substitute(JS,Dict(uS => equi_S)))
    JB_val = Symbolics.value.(substitute(JB,Dict(uB => equi_B)))

    uS_p = inv(JS_val) * A * (equi_S - equi_B[1:n])
    uB_p = inv(JB_val) * [B * (equi_B[1:n] - equi_S) ; zeros(m-n)]

    equi_p = [uS_p;uB_p]

    H_fun = eval(build_function(H,X,[uS;uB],δ;expression=Val{false}))[1]
    H1_val = H_fun(X_cr,[equi_S;equi_B],δ_cr)[2]

    ∂δ = Differential(δ)
    H0 = substitute(H[1],Dict(X => X_cr))
    dδ_H0 = expand_derivatives(∂δ(H0))
    dδ_H0_val = Symbolics.value(substitute(
        dδ_H0,
        Dict(
            uS => equi_S,
            uB => equi_B,
            δ => δ_cr)
    )
    )

    ∇u = [Differential(u) for u in [uS;uB]]
    ∇u_H0 = [expand_derivatives(∂u(H0)) for ∂u in ∇u]
    ∇u_H0_val = Symbolics.value.(substitute(
        ∇u_H0,
        Dict(
            uS => equi_S,
            uB => equi_B,
            δ => δ_cr)
    )
    )

    return - (H1_val + sum(equi_p .* ∇u_H0_val))/dδ_H0_val
end


function compute_TuringSpace!(
    zs,xs,ys,X_all,
    reac_fun,diff_fun,param,
    η_sb,δ_sb,
    n_var,
    equi_fun=nothing;
    U0=nothing,kwargs...)

    @variables var[1:n_var]

    X,η,δ,coeffs,_ = compute_coeffs(
            reac_fun,diff_fun,
            var,param,
            η_sb,δ_sb;
            kwargs...
    )
    coeffs_fun = eval(build_function(coeffs,X,var,η,δ;expression=Val{false}))[1]

    RH_fun = compute_RH_fun(n_var)

    if isnothing(U0)
        U0 = 10.0*rand(n_var)
    end
    
    if isnothing(equi_fun)
        equi0 = equilibria(reac_fun,param;U0=U0,kwargs...)
    else
        param[η_sb] = 0.0
        equi0 = equi_fun(param;U0=U0,kwargs...)
    end

    @showprogress for i in eachindex(xs)
        for j in eachindex(ys)
            param[δ_sb] = ys[j]
            param[η_sb] = xs[i]
            if isnothing(equi_fun)
                equi = equilibria(reac_fun,param;U0=equi0,kwargs...)
            else
                equi = equi_fun(param;U0=equi0,kwargs...)
            end
            zs[i,j] = unstable_modes(
                X_all,
                coeffs_fun,
                equi,xs[i],ys[j];
                RH_fun=RH_fun
            )
        end
    end

    return zs
end
