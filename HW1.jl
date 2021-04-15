using DifferentialEquations
using LinearAlgebra, Statistics, Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Distributions
using Formatting
using Printf


function define_env(;ρ̄      =  .051,
                     θᵨ     =   ,
                     θᵢ      =   ,
                     ϵ      =   10,
                     ψ      =   1.0,
                     ϕ      =   1.25,
                     γ      =   1.0)
    σ   =   1/γ
    params  =   @with_kw (ρ̄      =  ρ̄,
                     θᵨ     =   θᵨ,
                     θᵢ      =   θᵢ,
                     ϵ      =   ϵ,
                     ψ      =   ψ,
                     ϕ      =   ϕ,
                     σ      =   σ,
                     γ      =   γ)
    return params
end

params  =   define_env()