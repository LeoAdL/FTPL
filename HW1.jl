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


solution=solve_system_quad(;params=define_env(;κ=.1))



plot_IRF_quad(;pos =[1,3,5,7,8],solution=solution)

plot_θ_cum_quad(theta_range=range(.1,500,length=5))

plot_θ_cum_quad(;theta_range=range(.1,500,length=10),T_range=[1])

compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)



