using CairoMakie
using Symbolics
using LinearAlgebra

include("scr/dispersion.jl")
include("scr/RouthHurwitz.jl")
include("scr/solve.jl")
include("scr/compute.jl")

# 1. DEFINE A MODEL 

include("models.jl")

# 2. DEFINE THE (η,δ) Turing Space

n_var = 4    # Number of variables

################ THE SET OF MODES ###############################################

L = 1000.0    # Length of the space
p_all = collect(0:1200)    # Number of modes 
# p_all = [0.0]
X_all = (pi/L * p_all).^2

#################################################################################


################ BIFURCATION VARIABLES IN THE (η,δ) Turing space ################
### The Julia-symbols to access the variables in the dictionary of parameters ### 

η_sb = :η
δ_sb = :DvS

#################################################################################


################ THE SPATIAL DISCRETISATION OF THE (η,δ) Turing space ###########

xs = LinRange(0.0,5.0,200)    # η-axis
ys = LinRange(14.0,20.0,200)    # δ=axis

#################################################################################


################ THE REACTION, DIFFUSION, EQUILIRBIA FUNCTIONS AND KWARGS #######

reac_fun = reaction_cpl
diff_fun = diffusion_cpl
equi_fun = equilibria_cpl
kwargs = Dict(:A => A, :B => B, :n => size(A,1))

#################################################################################


# 3. COMPUTE AND PLOT STUFF.... 

#===============================================================================#

println("Compute Turing space....")

zs = zeros(length(xs),length(ys))

zs = compute_TuringSpace!(
    zs,xs,ys,X_all,
    reac_fun,diff_fun,param,
    η_sb,δ_sb,
    n_var,
    equi_fun;
    kwargs...
)

println("Plot Turing space....")

fig = Figure(resolution=(600,600),fontsize=28)
ax = Axis(
    fig[1,1],
    xgridcolor = :black,
    ygridcolor = :black,
    xgridwidth = 0.5,
    ygridwidth = 0.5,
    xminorgridcolor = :black,
    yminorgridcolor = :black,
    xminorgridvisible = true,
    yminorgridvisible = true,
    xminorgridwidth = 0.2,
    yminorgridwidth = 0.2,
    xminorticks = IntervalsBetween(5),
    yminorticks = IntervalsBetween(5),
    xlabel=String(η_sb),
    ylabel=String(δ_sb),
)
hm = heatmap!(ax,xs,ys,zs)
translate!(hm, 0, 0, -100)
Colorbar(fig[1,2],hm)
ylims!(ax,ys[1],ys[end])

#===============================================================================#


#===============================================================================#

println("Compute the critical values at η=0 ....")

param[η_sb] = 0.0
equi0 = equi_fun(param;U0=10.0*rand(n_var),kwargs...)
X_cr, δ_cr = critical(reac_fun,diff_fun,param,η_sb,δ_sb,equi0;kwargs...)
hlines!(ax,[δ_cr],color=:red,linestyle="--")    # add to current plot

#===============================================================================#

#===============================================================================#

println("Compute the slope at η=0 ....")

δ_p = compute_slope(param,η_sb,δ_sb,equi0,X_cr,δ_cr;kwargs...)
tgt = lines!(ax,xs,δ_cr.+δ_p.*xs,color=:red)    # add to current plot

#===============================================================================#

#===============================================================================#

# println("Compute the asymptotic critical values when η→+∞ ....")

U0 = rand(2)
equi = equilibria(reaction_BS,param;A=A,B=B,U0=U0)
X_BS_cr, δ_BS_cr = critical(reaction_BS,diffusion_BS,param,η_sb,δ_sb,equi;A=A,B=B,kwargs...)
hlines!(ax,[δ_BS_cr],color=:red)    # add to current plot

#===============================================================================#

#===============================================================================#

println("Compute the critical curve ....")

η_cc, δ_cc, X_cc = critical_curve(
    xs,
    reac_fun,diff_fun,param,η_sb,δ_sb,equi_fun;equi0,U0=[X_cr,δ_cr],kwargs...)
critic_curve = lines!(ax,η_cc,δ_cc,color=:red)     # add to current plot

#===============================================================================#

display(fig)
println("All done.")
