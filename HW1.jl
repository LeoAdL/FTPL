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


pp =define_env(T=200,N_t=200)

include("Function_quadratic.jl")
 solution=solve_system_quad(;params=pp)
plot_IRF_quad(;solution=solution)

plot_θ_cum_quad()

 plot_θ_cum_quad(;θ_range=range(.1,1000,length=10),T_range=[0,T],κ_range=[3,30,300])


compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)

pp =define_env(T=45,N_t=90,ind_Taylor=0.0,S=0.0
                ,long_term=0.0)

include("NK_FTPL_no_K.jl")

plot_all_longterm()

plot_all_S()
solutionNK_FTPL =solve_system(params=pp)

p_2=plot_IRF_FTPL(solution=solutionNK_FTPL)
plot(p_1,p_2,layout = 8)
plot_θ_cum_FTPL()

pp =define_env(T=150,N_t=150,ind_Taylor=0.0,long_term=0.0,S=0.0)

include("NK_FTPL_WITH_K.jl")


 solutionNK_FTPL_quad =solve_system_quad_FTPL(params=pp)

plot_IRF_quad_FTPL(solution=solutionNK_FTPL_quad)

plot_θ_cum_quad_FTPL(κ_range=[3,30,300])

plot_all_longterm()

plot_all_S()
