using DifferentialEquations
using BoundaryValueDiffEq
using LinearAlgebra, Statistics, Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Distributions
using Formatting
using Printf

#Schaab + KMV 2018


function define_env(;ρ̄      =  .022,
                     θᵨ     =   .22,
                     θᵢ      =   1.0,
                     θ      =   100.0,
                     ϵ      =   10,
                     ψ      =   2.0,
                     ϕ      =   1.25,
                     γ      =   2.0)
    σ   =   1/γ
    params  =   @with_kw (ρ̄      =  ρ̄,
                     θᵨ     =   θᵨ,
                     θᵢ      =   θᵢ,
                     θ      =   θ,
                     ϵ      =   ϵ,
                     ψ      =   ψ,
                     ϕ      =   ϕ,
                     σ      =   σ,
                     γ      =   γ,
                     T      =   10.0,
                     dt     =   .01)
    return params
end

params  =   define_env()

function NK_Rote!(du,u,p,t)
    @unpack σ,ϵ,θ,σ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   params
    x    =   u[1]
    dx   =   du[1]

    π   =   u[2]
    dπ  =   du[2]

    i   =   u[3]
    di  =   du[3]

    ρ   =   u[4]
    dρ  =   du[4]

    dx  =   σ*(i-π-ρ)*x
    
    dπ  =   (ϵ-1)/θ*(x^(1\σ+ψ)-1)-π*((1-σ)*(i-π)+σ*ρ)

    di  =   -θᵢ*(i-ϕ*π)

    dρ  =   -θᵨ*(ρ-ρ̄)
end