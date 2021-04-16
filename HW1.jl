using BoundaryValueDiffEq
using LinearAlgebra, Statistics, Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Distributions
using Formatting
using Printf
using DifferentialEquations



#Schaab + KMV 2018


include("Functions.jl")

solution=solve_system(;params=define_env())



plot_IRF(;pos =[1,2,3,4],solution=solution)

plot_θ()

compute_dev(;θ=1.0,T=50.0)

compute_dev(;θ=10.0,T=50.0)


compute_dev(;θ=0.1,T=50.0)



