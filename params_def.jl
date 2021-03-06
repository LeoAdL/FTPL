#Schaab + KMV 2018


function define_env(;ρ̄      =  2.2/100,
                     θᵨ         = .22,
                     σᵨ         = 0.003,
                     θᵢ         = 1.0,
                     θ          = 100.0,
                     ϵ          = 11,
                     ψ          = 1/2.0,
                     ϕ          = 1.25,
                     γ          = 2.0,
                     T          = 50,
                     N_t        = 80,
                     κ          = 3.0,
                     δ          = .1,
                     α          = 1/3,
                     ϕ_FTPL     = 1.0-.25,
                     s₀         = 3.0/100,
                     ind_Taylor = 0.0,
                     A          = .133,
                     S          = 0.1,
                     long_term  = 0.0)
              #S       = 2*ρ̄ +1
              i_target = 2.0/100 + ρ̄
              init_ρ   = ρ̄-sqrt(σᵨ^2/(2*θᵨ^2))/10
              σ        = 1.0/γ
              dt       = T/N_t
              χₙ       = A*(ϵ-1.0)/ϵ
           di(i,π)     = -θᵢ*(i-ϕ_FTPL*π)*ind_Taylor

       params   = @with_kw (ρ̄      =  ρ̄,
                     θᵨ         = θᵨ,
                     θᵢ         = θᵢ,
                     θ          = θ,
                     ϵ          = ϵ,
                     ψ          = ψ,
                     ϕ          = ϕ,
                     σ          = σ,
                     γ          = γ,
                     T          = T,
                     κ          = κ,
                     δ          = δ,
                     α          = α,
                     A          = A,
                     dt         = dt,
                     χₙ         = χₙ,
                     S          = S,
                     s₀         = s₀,
                     i_target   = i_target,
                     init_ρ     = init_ρ,
                     ind_Taylor = ind_Taylor,
                     di         = di,
                     ϕ_FTPL     = ϕ_FTPL,
                     long_term  = long_term)
    return params
end

