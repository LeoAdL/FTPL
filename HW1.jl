using BoundaryValueDiffEq
using LinearAlgebra, Statistics, Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Distributions
using Formatting
using Printf

#Schaab + KMV 2018


include("Functions.jl")

solution=solve_system(;params=define_env(θ=1))



plot_IRF(;pos =[1],sol=solution.sol,SS=solution.SS)

plot_θ()

compute_dev(;θ=400.0,T=10.0)


