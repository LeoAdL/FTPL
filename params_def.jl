#Schaab + KMV 2018


function define_env(;ρ̄      =  2.2/100,
                     θᵨ     =   .22,
                     σᵨ     =   0.003,
                     θᵢ      =   1.0,
                     θ      =   100.0,
                     ϵ      =   11,
                     ψ      =   1/2.0,
                     ϕ      =   1.25,
                     γ      =   2.0,
                     T      =   50.0,
                     N_t    =   100.0,
                     κ      =   3.0,
                     δ      =   .1,
                     α      =   1/3)
    init_ρ =   ρ̄+sqrt(σᵨ^2/(2*θᵨ^2))
    σ   =   1/γ
    dt  =   T/N_t
    
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
                     κ      =   κ,
                     δ      =   δ,
                     α      =   α,
                     dt     =   dt,
                     init_ρ=init_ρ)
    return params
end
