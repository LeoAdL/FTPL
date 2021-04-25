using BoundaryValueDiffEq
using Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Formatting
using Printf
using DifferentialEquations



include("params_def.jl")
params =define_env(T=100,N_t=100)
include("Functions.jl")

solutionNK=solve_system()

plot_IRF(solution=solutionNK,var=["x","\\pi","i"])

params =define_env()
include("Function_quadratic.jl")

plot_IRF_quad(;var =["k","\\iota","Y"],
                solution=solve_system_quad(;params=define_env(κ=10)),
                T_end=200)

@time plot_θ_cum_quad(;var="Y",theta_range=range(.1,500,length=2),κ_range=[10,10000],T_range=[0])

 plot_θ_cum_quad(;θ_range=range(.1,500,length=2),T_range=[1,T],κ_range=[30,300])

 plot_θ_cum(;var="x",θ_range=range(10^(-3),500,length=2),ϕ=ϕ,
                T_range=[T])

compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)



