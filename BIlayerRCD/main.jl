using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Symbolics, IfElse
using DifferentialEquations
using JLD2
using NLsolve


println("Initialize everything...")

#-------- Add the scripts------------------------------------------------------#
include("reaction_taxis_growth.jl")
include("scr/plot_save.jl")
include("scr/ODE_func.jl")
#------------------------------------------------------------------------------#

simu_name = init_directory(simu_name="simu")

#-------- Discretization parameters -------------------------------------------#
dim = 2    # Dimension of the bulk (1 or 2)
Nx = 128    # Number of cells in the x-direction (bulk and surface)
Ny = 128    # Number of cells in the y-direction (bulk)
Lx = 200.0     # Length of the domain in the x-direction (bulk and surface)
Ly = 200.0   # Length of the domain in the y-direction (bulk)
dx = Lx/Nx    # Cell length in the x-direction
dy = Ly/Ny    # Cell length in the y-direction

#------------------------------------------------------------------------------#


#-------- Coupling parameters ----------------------------------------------#

## Transfer coefficient at the S/B boundary
η = 0.01
# η = 0.0    # No coupling between the surface and the bulk

## Exchange coefficients
αS = 1.0
βS = 1.0
αB = 1.0
βB = 180.0

#------------------------------------------------------------------------------#


#-------- Time parameters -----------------------------------------------------#

T = 2000.0 # Simulation time 

#------------------------------------------------------------------------------#


#------------------- Model  ---------------------------------------------------#

#=================== SURFACE ==========================#

#********************************#
# reac_S_name = "Schnakenberg"
# a_S = 0.2305 
# b_S = 0.7695
# sc_S = 1.0
# param_reac_S = Dict(
#     :a => a_S,
#     :b => b_S,
#     :sc => sc_S
# )

# Diffusion parameters
# DuS = 1.0
# DvS = 20.0

## Critical diffusion
# fuS = param_reac_S[:sc]*(-1 + 2*uS_eq*vS_eq)
# fvS = param_reac_S[:sc]*uS_eq^2
# guS = param_reac_S[:sc]*(-2*uS_eq*vS_eq)
# gvS = param_reac_S[:sc]*(-uS_eq^2)
# DvS_cr = ((sqrt(fuS*gvS - fvS*guS) + sqrt(-fvS*guS))/fuS)^2

## Equilibria
# uS_eq = a_S + b_S 
# vS_eq = b_S/(a_S + b_S)^2
#********************************#


#********************************#
reac_S_name = "Linear Activator-Inhibitor"
synth_uS = 1.0
a1S = 2.0
bS = 2.0
a2S = 1.0
synth_vS = 0.25
c1S = 1.0
dS = 1.0
c2S = 2.0
param_reac_S = Dict(
    :act_synth => synth_uS,
    :act_autocat => a1S,
    :inhib_act =>  bS,
    :max_act_synth => 10.0,
    :act_deg => a2S,
    :inhib_synth => synth_vS,
    :inhib_autocat => c1S,
    :inhib_act_crea => dS,
    :max_inhib_synth => 10.0,
    :inhib_deg => c2S,
    :sc => 1.0
)


# Equilibria 
aS = a1S - a2S
cS = c2S - c1S
uS_eq = 1/(aS*cS - bS*dS) * (bS*synth_vS - synth_uS*cS)
vS_eq = 1/(aS*cS - bS*dS) * (aS*synth_vS - synth_uS*dS)

## Diffusion parameters
# DvS_cr = ((sqrt(bS*dS - aS*cS) + sqrt(bS*dS))^2 / aS^2)
DvS = 2.91
DuS = 1.0
#********************************#

#======================================================#


#===================== BULK ===========================#

#********************************#
# reac_B_name = "Schnakenberg"
# a_B = 0.2305 
# b_B = 0.7695
# sc_B = 1.0

# param_reac_B = Dict(
#     :a => a_B,
#     :b => b_B,
#     :sc => sc_B,
#     :c_eq => 1.0
# )

## Equilibria 
# uB_eq = a_B + b_B 
# vB_eq = b_B/(a_B + b_B)^2

## Diffusion parameters
# DuB = 1.0
# DvB = 15.0
# Dc = 1.0
# fuB = param_reac_B[:sc]*(-1 + 2*uB_eq*vB_eq)
# fvB = param_reac_B[:sc]*uB_eq^2
# guB = param_reac_B[:sc]*(-2*uB_eq*vB_eq)
# gvB = param_reac_B[:sc]*(-uB_eq^2)
# DvB_cr = ((sqrt(fuB*gvB - fvB*guB) + sqrt(-fvB*guB))/fuB)^2
#********************************#

#********************************#
reac_B_name = "Linear Activator-Inhibitor"
param_reac_B = Dict(
    :act_synth => 3.0,
    :act_autocat => 8.0,
    :inhib_act =>  2.0,
    :max_act_synth => 10.0,
    :act_deg => 7.0,
    :inhib_synth => 0.75,
    :inhib_autocat => 0.5,
    :inhib_act_crea => 1.0,
    :max_inhib_synth => 10.0,
    :inhib_deg => 1.5,
    :sc => 1.0,
    :c_eq => 0.25
)


DuB = 25.0
DvB = 35.0
Dc = 25.0
# aB = param_reac_B[:act_autocat] - param_reac_B[:act_deg]
# cB = param_reac_B[:inhib_deg] - param_reac_B[:inhib_autocat]
# a2B = param_reac_B[:act_deg]
# c2B = param_reac_B[:inhib_deg]
# DvB_cr = ((sqrt(bB*dB - aB*cB) + sqrt(bB*dB))^2 / aB^2) * DuB
#********************************#

## Chemotaxis
# χ = 0.0
χ = 120.0

## Growth rate
γ = 0.2   # Growth rate
# γ = 0.0

#=======================================================#

#------------------- Initial condition  ---------------------------------------#

## Compute the equilibrium
function f!(F,U;η=η,αS=αS,βS=βS,αB=αB,βB=βB)
    r_uS, r_vS = reaction_S(U[1],U[2],param_reac_S)
    r_uB, r_vB = reaction_B(U[3],U[4],param_reac_B[:c_eq],param_reac_B)
    F[1] = r_uS + η*αS*(U[3]-U[1])
    F[2] = r_vS + η*βS*(U[4]-U[2])
    F[3] = r_uB + η*αB*(U[1]-U[3])
    F[4] = r_vB + η*βB*(U[2]-U[4])
end
# equi = nlsolve(f!, [uB_eq; vB_eq; uS_eq; vS_eq])
global equi = nlsolve(f!,[rand(); rand(); rand(); rand()])
while ~equi.f_converged || equi.zero[1] < 0 || equi.zero[2] < 0 || equi.zero[3] < 0 || equi.zero[4] < 0
    global equi = nlsolve(f!, [500*rand(); 500*rand(); 500*rand(); 500*rand()])
end

uS_eq = equi.zero[1]
vS_eq = equi.zero[2]
uB_eq = equi.zero[3]
vB_eq = equi.zero[4]

###

if dim==2
    U_init = zeros(Nx,Ny+1,3)
    # The solution U is gathered in a big array of size (Nx,Ny+1,3). 
    # The individual solutions can be recovered with the following index table
    uS_index = :,1,1    
    vS_index = :,1,2    
    uB_index = :,2:(Ny+1),1 
    vB_index = :,2:(Ny+1),2
    c_index = :,2:(Ny+1),3

    U_init[uS_index...] = uS_eq * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
    U_init[vS_index...] = vS_eq * (ones(Nx) + 0.1 * (rand(Nx).-0.5))

    U_init[uB_index...] = uB_eq * (ones(Nx,Ny) + 0.1 * (rand(Nx,Ny).-0.5))
    U_init[vB_index...] = vB_eq * (ones(Nx,Ny) + 0.1 * (rand(Nx,Ny).-0.5))
    U_init[c_index...] = param_reac_B[:c_eq] * (ones(Nx,Ny) + 0.1 * (rand(Nx,Ny).-0.5))
else
    U_init = zeros(Nx,5)
    # The solution U is gathered in a big array of size (Nx,5). 
    # The individual solutions can be recovered with the following index table
    uS_index = :,1   
    vS_index = :,2    
    uB_index = :,3
    vB_index = :,4 
    c_index = :,5 

    U_init[uS_index...] = equi.zero[1] * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
    U_init[vS_index...] = equi.zero[2] * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
    U_init[uB_index...] = equi.zero[3] * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
    U_init[vB_index...] = equi.zero[4] * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
    U_init[c_index...] = param_reac_B[:c_eq] * (ones(Nx) + 0.1 * (rand(Nx).-0.5))
end


#------------------------------------------------------------------------------#


#--------------------- Create an ODE problem ----------------------------------#

if dim==2
    p = (Dc,DuS,DvS,DuB,DvB,χ,γ,η,dx,dy,Nx,Ny,param_reac_S,param_reac_B,αS,βS,αB,αB)    # Simulation parameters
    func! = func_2D!
else
    p = (DuS,DvS,DuB,DvB,Dc,χ,γ,η,dx,Nx,param_reac_S,param_reac_B,αS,βS,αB,βB)
    func! = func_1D!
end

println("Initialize ODE problem...")

dU_init = copy(U_init)
jac_sparsity = Symbolics.jacobian_sparsity((du,u)->func!(du,u,p,0.0),dU_init,U_init)    # Automatic sparsity detection

ode_func = ODEFunction(func!,jac_prototype=float.(jac_sparsity))
prob = ODEProblem(ode_func,U_init,(0.,T),p)

# prob = ODEProblem(func!,U_init,(0.,T),p)    # Without sparsity information

println("Ok let's solve the ODE...")

# algo = nothing    # Default choice
# algo = Euler() 
algo = ImplicitEuler()
# algo = RadauIIA3()

#------------------------------------------------------------------------------#


#----------- Finallly run the simulation and save data ------------------------#

sol = solve(prob,algo,alg_hints=[:stiff],saveat=1.0;progress=true,progress_steps=1)
# sol = solve(prob,algo,saveat=10.0;dt=0.2,progress=true,progress_steps=1)    # For explicit Euler, specify dt

println("Solving is done, now plotting...")

plot_save_sol(
    sol,U_init,c_index,uB_index,vB_index,uS_index,vS_index,dx,dy;
    range_uS=(0.498,0.506),
    range_vS=(0.749,0.755),
    range_uB=(0.4,0.9),
    range_vB=(0.69,0.86),
    range_c=(0.19,1.1),
    dir=simu_name,video_name=simu_name,
    dt=1.0)

println("Save data...")

data = Dict(
    "simu_name" => simu_name,
    "reac_B_name" => reac_B_name,
    "reac_S_name" => reac_S_name,
    "sol_final" => sol[end],
    "T" => T,
    "Dc" => Dc,
    "DuS" => DuS,
    "DvS" => DvS,
    "DuB" => DuB, 
    "DvB" => DvB,
    "χ" => χ,
    "γ" => γ,
    "η" => η,
    "dx" => dx,
    "dy" => dy,
    "Nx" => Nx,
    "Ny" => Ny,
    "param_reac_S" => param_reac_S,
    "param_reac_B" => param_reac_B,
    "αS" => αS,
    "βS" => βS,
    "αB" => αB,
    "βB" => βB,
)

save(simu_name*"/"*simu_name*".jld2",data)

println("All done.")

