using BoundaryValueDiffEq
using Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Formatting
using Printf
using DifferentialEquations



include("params_def.jl")
pp =define_env(T=30,N_t=60)
include("Functions.jl")

solutionNK=solve_system(params =pp)



plot_IRF(solution=solutionNK,var =["x","\\pi","i"])

 plot_θ_cum(;var="x",θ_range=range(.1,2000,length=75),ϕ=ϕ,
                T_range=[0,1,5,10,T])


pp =define_env(T=300,N_t=150)

include("Function_quadratic.jl")
@time solution=solve_system_quad(;params=pp)
plot_IRF_quad(;solution=solution,var=["\\iota"])

@time plot_θ_cum_quad(;var="Y",theta_range=range(.1,500,length=2),κ_range=[10,10000],T_range=[0])

 plot_θ_cum_quad(;θ_range=range(.1,500,length=2),T_range=[0,T],κ_range=[30,300])


compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)

pp =define_env(T=40,N_t=70,ind_Taylor=0.0,S=0.0)

include("NK_FTPL_no_K.jl")



@time solutionNK_FTPL =solve_system(params=pp)

plot_IRF_FTPL(solution=solutionNK_FTPL)

plot_θ_cum(θ_range=range(.01,150,length=20),T_range=[0,30])

pp =define_env(T=70,N_t=30,ind_Taylor=0.0,S=0.0)

include("NK_FTPL_WITH_K.jl")


@time solutionNK_FTPL_quad =solve_system_quad_FTPL(params=pp)

plot_IRF_quad_FTPL(solution=solutionNK_FTPL_quad)