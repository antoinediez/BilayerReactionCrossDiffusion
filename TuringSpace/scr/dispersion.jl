using LinearAlgebra

function dispersion(λ,X,J,D)
    return det(λ*I - J + X.*D)
end

function λ_coeff(Δ,λ,K)
    ∂λ = Differential(λ)
    coeffs = Num[]
    for k in 0:K
        if k==0
            append!(coeffs,substitute(Δ,Dict(λ => 0.0)))
        else
            append!(coeffs,substitute(expand_derivatives((∂λ^k)(Δ)),Dict(λ => 0.0))/factorial(k))
        end
    end
    return coeffs
end

function η_coeff(Δ,λ,η,K)
    ∂η = Differential(η)
    coeffs = Num[]
    for k in 0:K
        if k==0
            append!(coeffs,substitute(Δ,Dict(λ => 0.0, η => 0.0)))
        else
            append!(coeffs,substitute(expand_derivatives((∂η^k)(Δ)),Dict(λ => 0.0, η => 0.0))/factorial(k))
        end
    end
    return coeffs
end

function dispersion_cpl(λ,X,η,JS,DS,JB,DB,A,B)
    n = size(JS)[1]
    m = size(JB)[1]

    J = [
        JS         zeros(n,m);
        zeros(m,n) JB
    ]

    D = [
        DS         zeros(n,m);
        zeros(m,n) DB
    ]

    E = [
        A           (-A)          zeros(n,m-n);
        (-B)         B            zeros(n,m-n);
        zeros(m-n,n) zeros(m-n,n) zeros(m-n,m-n)   
    ]

    return det(-J + X.*D + η.*E + λ*I)
end
