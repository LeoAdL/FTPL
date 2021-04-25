using BoundaryValueDiffEq
using Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Formatting
using Printf
using DifferentialEquations



include("params_def.jl")
include("Functions.jl")
include("Function_quadratic.jl")

@time solution=solve_system_quad(;params=define_env(κ=10))
plot_IRF_quad(;var =["k","\\iota","Y"],solution=solution,T_end=200)

@time plot_θ_cum_quad(;var="Y",theta_range=range(.1,500,length=2),κ_range=[10,10000],T_range=[0])

plot_θ_cum_quad(;theta_range=range(.1,500,length=10),T_range=[1])

compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)



