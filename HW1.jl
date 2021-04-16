using BoundaryValueDiffEq
using LinearAlgebra, Statistics, Plots
using Parameters,LaTeXStrings
using NLsolve
using StatsBase, Random
using Distributions
using Formatting
using Printf

#Schaab + KMV 2018


function define_env(;ρ̄      =  .0022,
                     θᵨ     =   .22,
                     θᵢ      =   1.0,
                     θ      =   100.0,
                     ϵ      =   10,
                     ψ      =   2.0,
                     ϕ      =   1.25,
                     γ      =   2.0,
                     T      =   10.0,
                     σᵨ     =   0.003)
    init_ρ =   ρ̄+sqrt(σᵨ^2/(1-θᵨ^2))
    σ   =   1/γ
    dt  =   T/100.0
    params  =   @with_kw (ρ̄      =  ρ̄,
                     θᵨ     =   θᵨ,
                     θᵢ      =   θᵢ,
                     θ      =   θ,
                     ϵ      =   ϵ,
                     ψ      =   ψ,
                     ϕ      =   ϕ,
                     σ      =   σ,
                     γ      =   γ,
                     T      =   T,
                     dt     =   dt,
                     init_ρ=init_ρ)
    return params
end

params  =   define_env()

function NK_Rote!(du,u,p,t)
    @unpack σ,ϵ,θ,σ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   params
    x    =   u[1]

    π   =   u[2]

    i   =   u[3]

    ρ   =   u[4]

    du[1]  =   σ*(i-π-ρ)*x
    
    du[2]  =   (ϵ-1.0)/θ*(x^(1\σ+ψ)-1.0)-π*((1.0-σ)*(i-π)+σ*ρ)

    du[3]  =   -θᵢ*(i-ϕ*π)

    du[4]  =   -θᵨ*(ρ-ρ̄)
    
end

function SS()
    @unpack σ,ϵ,θ,σ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   params

    ρ_ss    =   ρ̄
    π_ss    =   ρ̄/(ϕ-1.0)
    i_ss    =   ρ̄*ϕ/(ϕ-1.0)
    x_ss    =   (1.0+θ/(ϵ-1.0)*π_ss*(σ*ρ̄+(1-σ)*(i_ss-π_ss)))^(1.0/(1.0/σ+ψ))

    return (ρ_ss=ρ_ss,π_ss=π_ss,i_ss=i_ss,x_ss=x_ss)
end

function bc1!(residual,u,p,t)
    @unpack ρ_ss,π_ss,i_ss,x_ss= SS()
    @unpack init_ρ  =   params
    residual[1] =   u[end][1]- x_ss
    residual[2] =   u[end][2]- π_ss
    residual[3] =   u[1][3]- i_ss
    residual[4] =   u[1][4]- init_ρ
end

@unpack ρ_ss,π_ss,i_ss,x_ss= SS()
@unpack T,dt,init_ρ = params

tspan   =   (0.0,T)


bvp1 = BVProblem(NK_Rote!, bc1!, [x_ss,π_ss,i_ss,init_ρ], tspan)

sol1 = solve(bvp1, GeneralMIRK4(), dt=dt)

SS_vec = [x_ss,π_ss,i_ss,ρ_ss]
dev =   (sol1[1:end].-SS_vec)./SS_vec
plot(sol1.t,dev')
