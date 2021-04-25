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

@time solution=solve_system_quad(;params=define_env())
plot_IRF_quad(;solution=solution)

plot_θ_cum_quad(;var="Y",theta_range=range(.1,100,length=3),κ_range=[.5,100])

plot_θ_cum_quad(;theta_range=range(.1,500,length=10),T_range=[1])

compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)



