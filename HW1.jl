using BoundaryValueDiffEq
using Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Formatting
using Printf
using DifferentialEquations




include("params_def.jl")
params=define_env(;γ=1.0)
include("Functions.jl")
include("Function_quadratic.jl")


solution=solve_system_quad(;params=define_env())



plot_IRF_quad(;solution=solution)

plot_θ_cum()

compute_dev(;θ=10^(-5),T=10.0)

compute_dev(;θ=1000.0,T=50.0)


compute_dev(;θ=100,T=50.0)

compute_dev(;θ=0.001,T=50.0)



