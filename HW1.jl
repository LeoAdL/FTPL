using BoundaryValueDiffEq
using Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Formatting
using Printf
using DifferentialEquations

include("params_def.jl")
pp =define_env(T=40,N_t=80)



include("Functions.jl")

@time solutionNK=solve_system(params =pp)



plot_IRF(solution=solutionNK,var =["x","\\pi","i"])

plot_θ_cum()


pp =define_env(T=200,N_t=30)

include("Function_quadratic.jl")
@time solution=solve_system_quad(;params=pp)
plot_IRF_quad(;solution=solution)

@time plot_θ_cum_quad()

 plot_θ_cum_quad(;θ_range=range(.1,500,length=2),T_range=[0,T],κ_range=[30,300])


compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)

pp =define_env(T=40,N_t=80,ind_Taylor=1.0,S=1.0,long_term=1.0)

include("NK_FTPL_no_K.jl")



@time solutionNK_FTPL =solve_system(params=pp)

plot_IRF_FTPL(solution=solutionNK_FTPL)

plot_θ_cum_FTPL()

pp =define_env(T=80,N_t=100,ind_Taylor=1.0,long_term=1.0)

include("NK_FTPL_WITH_K.jl")


@time solutionNK_FTPL_quad =solve_system_quad_FTPL(params=pp)

plot_IRF_quad_FTPL(solution=solutionNK_FTPL_quad)

plot_θ_cum_quad_FTPL(κ_range=[3,30,300])