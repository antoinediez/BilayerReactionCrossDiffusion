using StaticArrays
using LinearAlgebra


#===================== TWO SCHNAKENBERG ========================================#

param = Dict{Symbol,Float64}(
    :aS => 0.2305,
    :bS => 0.7695,
    :scS => 2.0,
    :DuS => 1.0,
    :DvS => 15.0,
    :aB => 0.15,
    :bB => 0.2,
    :scB => 0.5,
    :DuB => 1.0,
    :DvB => 15.0
)

function reaction_S(var,param;kwargs...)
    # u,v = var
    f = param[:scS]*(param[:aS] - var[1] + var[1]^2*var[2])
    g = param[:scS]*(param[:bS] - var[1]^2*var[2]) 
return f,g
end

function diffusion_S(var,param;kwargs...)
    D = [
        param[:DuS] 0.0;
        0.0        param[:DvS]
    ]
    return D
end

function equilibria_S(param;kwargs...)
    u = param[:aS] + param[:bS]
    v = param[:bS]/u^2
    return u,v,param[:ceq]
end

#########################

function reaction_B(var,param;kwargs...)
    # u,v,c = var
    f = param[:scB]*(param[:aB] - var[1] + var[1]^2*var[2])
    g = param[:scB]*(param[:bB] - var[1]^2*var[2])
return f,g
end

function diffusion_B(var,param;kwargs...)
    D = [
        param[:DuB] 0.0;         
        0.0         param[:DvB]]
    return D
end

function equilibria_B(param;kwargs...)
    u = param[:aB] + param[:bB]
    v = param[:bB]/u^2
    return u,v
end

#===============================================================================#


#===================== TWO PSEUDO-LINEAR + CHEMO ===============================#

# param = Dict{Symbol,Float64}(
#     :a0S => 1.0,
#     :a1S => 2.0,
#     :a2S => 1.0,
#     :bS => 2.0,
#     :b0S => 0.25,
#     :c1S => 1.0,
#     :c2S => 2.0, 
#     :dS => 1.0,
#     :maxS => 10000.0,
#     :scS => 1.0,
#     :DuS => 1.0,
#     :DvS => 2.91,
#     :a0B => 3.0,
#     :a1B => 8.0,
#     :a2B => 7.0,
#     :bB => 2.0,
#     :b0B => 0.75,
#     :c1B => 0.5,
#     :c2B => 1.5, 
#     :dB => 1.0,
#     :maxB => 10000.0,
#     :scB => 1.0,
#     :DuB => 25.0,
#     :DvB => 35.0,
#     :DcB => 25.0,
#     :γ => 0.2,
#     :ceq => 0.25,
#     :χ => 120.0
# )

# # function smooth_ReLU(x;β=1.0)
# #     # return log(1.0 + exp(β*x)) / β
# #     return 0.5 * (x + sqrt(x^2 + 0.1))
# # end

# # function ϕ(x,M;β=1.0)
# #     return M - smooth_ReLU(-smooth_ReLU(x;β=β)+M;β=β)
# # end

# function ϕ(x,M)
#     return max(0.0,min(x,M))
# end


# function reaction_S(var,param;kwargs...)
#     u_synth = ϕ(param[:a0S] + param[:a1S]*var[1] - param[:bS]*var[2],param[:maxS])
#     f = param[:scS]*(u_synth - param[:a2S]*var[1])
#     v_synth = ϕ(param[:b0S] + param[:c1S]*var[2] + param[:dS]*var[1],param[:maxS])
#     g = param[:scS]*(v_synth - param[:c2S]*var[2])
# return f,g
# end

# function diffusion_S(var,param;kwargs...)
#     D = [
#         param[:DuS] 0.0;
#         0.0        param[:DvS]
#     ]
#     return D
# end

# #####################

# function reaction_B(var,param;kwargs...)
#     u_synth = ϕ(param[:a0B] + param[:a1B]*var[3] - param[:bB]*var[2],param[:maxB])
#     f = param[:scB]*(u_synth - param[:a2B]*var[1])
#     v_synth = ϕ(param[:b0B] + param[:c1B]*var[3] + param[:dB]*var[3],param[:maxB])
#     g = param[:scB]*(v_synth - param[:c2B]*var[2])
#     r = param[:γ] * var[3] * (param[:ceq] - var[3])
# return f,g,r 
# end

# function diffusion_B(var,param;kwargs...)
#     D = [
#         param[:DuB]              0.0         0.0;
#         0.0                      param[:DvB] 0.0;
#         (-param[:χ]*param[:ceq]) 0.0         param[:DcB]
#     ]
#     return D
# end



#===============================================================================#



#########################

A = @SArray [
    1.0 0.0;
    0.0 1.0
]

B = @SArray [
    1.0 0.0;
    0.0 1.0
]

function reaction_cpl(var,param;A=A,B=B,kwargs...)
    n = size(A,1)
    m = length(var) - n
    var_S = @view var[1:n]
    var_B = @view var[(n+1):(n+m)]
    r_S = reaction_S(var_S,param)
    r_B = reaction_B(var_B,param)
    Δ = SVector{n,typeof(var[1])}([var[size(A,1)+i] - var[i] for i in 1:size(A,1)])
    cpl_S = A*Δ
    cpl_B = B*Δ
    F = zeros(typeof(var[1]),n+m)
    for i in 1:n
        F[i] = r_S[i] + param[:η]*cpl_S[i]
        F[n+i] = r_B[i] - param[:η]*cpl_B[i]
    end
    for i in (n+1):m
        F[n+i] = r_B[i]
    end
    return F
end

function diffusion_cpl(var,param;n=size(A,1),kwargs...)
    m = length(var) - n
    var_S = @view var[1:n]
    var_B = @view var[(n+1):(n+m)]
    DS = diffusion_S(var_S,param)
    DB = diffusion_B(var_B,param)
    return [DS zeros(n,m);zeros(m,n) DB]
end

#############################################

function reaction_BS(var,param;A=A,B=B,kwargs...)
    n = size(A,1)
    m = length(var)
    invAp = [inv(A) zeros(n,m-n) ; zeros(m-n,n) zeros(m-n,m-n)]
    invBp = [inv(B) zeros(n,m-n) ; zeros(m-n,n) I]
    L = [(inv(A) + inv(B)) zeros(n,m-n); zeros(m-n,n) I]
    return inv(L)*(invAp * [collect(reaction_S(var[1:n],param));zeros(m-n)] .+ invBp * collect(reaction_B(var,param)))
end

function diffusion_BS(var,param;A=A,B=B,kwargs...)
    n = size(A,1)
    m = length(var)
    invAp = [inv(A) zeros(n,m-n) ; zeros(m-n,n) zeros(m-n,m-n)]
    invBp = [inv(B) zeros(n,m-n) ; zeros(m-n,n) I]
    L = [(inv(A) + inv(B)) zeros(n,m-n); zeros(m-n,n) I]
    DSp = [diffusion_S(var[1:n],param) zeros(n,m-n) ; zeros(m-n,n) zeros(m-n,m-n)]
    return inv(L)*(invAp * DSp  .+ invBp * diffusion_B(var,param))
end