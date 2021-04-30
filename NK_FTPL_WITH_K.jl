function solve_system_quad_FTPL(;params)

    function static_funct(p)
    @unpack σ, ϵ, θ, ϕ, ψ, ρ̄, θᵨ, θᵢ, κ, δ, A, χₙ, α, γ,S, ind_Taylor, i_target, ϕ_FTPL, s₀ = p
        function ι(q::Float64)
            return((q-1.0)/κ)
        end

        function ℓ(C::Float64,k::Float64,q::Float64)
            return(k*((C/k+ι(q))/A)^(1.0/(1.0-α)))
        end

        function w(C::Float64,k::Float64,q::Float64)
            return (ℓ(C,k,q)^(ψ)*C^(γ))
        end


        function ν_k(C::Float64,k::Float64,q::Float64)
            return((1.0/q)*(α/(1.0-α))*w(C,k,q)*(ℓ(C,k,q)/k))
        end
            
        function χ(C::Float64,k::Float64,q::Float64)
            return((w(C,k,q)/(1.0-α))^(1.0-α)*(q*ν_k(C,k,q)/α)^(α))
        end

        function Y(C::Float64,k::Float64,q::Float64)
            return(A*(k)^(α)*(ℓ(C,k,q))^(1-α))
        end
        return (ι=ι,ℓ=ℓ,w=w,ν_k=ν_k,χ=χ,Y=Y)
    end



    function NK!(du,u,p,t)
         σ      ,ϵ, θ, ϕ, ψ,   ρ̄, θᵨ, θᵢ, κ, δ, A, χₙ,S, ind_Taylor, i_target, ϕ_FTPL, s₀ = p
         @unpack ι, ℓ, w, ν_k, χ,  Y  = static_funct(p)
                 q  = u[1]
                 C  = u[2]
                 k  = u[3]
                 π  = u[4]
                 ρ  = u[5]
                 i  = u[6]
                 v  = u[7]
                 vˡ  = u[8]


            du[1] = q*(i-π+(ι(q)+κ*(ι(q))^(2.0)/2.0)/q-ι(q)-ν_k(C,k,q))
            
            du[2] = σ*C*(i-π-ρ)
            
            du[3] = ι(q)*k

            du[4] = π*((1.0-σ)*(i-π)+σ*ρ)-((ϵ-1.0)/θ)*(χ(C,k,q)/χₙ-1.0)
            
            du[5] = -θᵨ*(ρ-ρ̄)

            du[6] = -θᵢ*(i-ϕ_FTPL*π)*ind_Taylor

            du[7] = v*(i-π) -(s₀+S*vˡ)

            du[8] = vˡ *(i-π) -(s₀+S*vˡ)

        
    end

    function SS(p)
        @unpack α, γ, σ, ϵ, θ, ϕ, ψ, ρ̄, θᵨ, θᵢ, κ, δ, A,S, ind_Taylor, i_target, ϕ_FTPL, s₀ = p

            k_c = (α/(ρ̄))*((ϵ-1.0)/ϵ)*(1.0+θ/(ϵ-1.0)*ρ̄^(2)/(ϕ_FTPL*ind_Taylor-1.0))
    
            k_l = (A*k_c)^(1.0/(1.0-α))

            q_ss = 1.0
            k_ss = (ρ̄*((1-α)/α)*(k_l)^(1+ψ)*(k_c)^(γ))^(1/(ψ+γ))
            C_ss = k_ss/k_c
            ρ_ss = ρ̄
            i_ss = ρ̄*ϕ_FTPL/(ϕ_FTPL-1.0)*ind_Taylor +(1-ind_Taylor)*i_target
            π_ss = i_ss - ρ_ss

            v_ss = s₀/(i_ss-π_ss-S)
            s_ss = s₀ + S*v_ss

       return(π_ss=π_ss,
                C_ss = C_ss,
                q_ss = q_ss,
                k_ss = k_ss,
                ρ_ss = ρ_ss,
                i_ss = i_ss,
                ℓ_ss = k_ss/k_l,
                ι_ss = 0.0,
                v_ss = v_ss,
                s_ss = s_ss)
    end

    function u_0(p)
        @unpack q_ss,   C_ss, k_ss, π_ss, i_ss, ℓ_ss, ι_ss,v_ss = SS(p)
        @unpack init_ρ = p
    return ([q_ss,C_ss,k_ss,π_ss,init_ρ,i_ss,v_ss,v_ss])
    end

    p  =    (σ=params.σ,
            ϵ          = params.ϵ,
            θ          = params.θ,
            ϕ          = params.ϕ,
            ψ          = params.ψ,
            ρ̄         = params.ρ̄,
            θᵨ         = params.θᵨ,
            θᵢ         = params.θᵢ,
            κ          = params.κ,
            δ          = params.δ,
            A          = params.A,
            χₙ         = params.χₙ,
            γ          = params.γ,
            α          = params.α,
            init_ρ     = params.init_ρ,
            S          = params.S,
            i_target   = params.i_target,
            ind_Taylor = params.ind_Taylor,
            ϕ_FTPL     = params.ϕ_FTPL,
            s₀         = params.s₀)



    function bc1!(residual,u,p,t)
            @unpack  q_ss,   C_ss, k_ss, π_ss, ρ_ss, i_ss,v_ss = SS(p)
            @unpack  init_ρ = p
            residual[1] = u[end][1]- q_ss
            residual[2] = u[end][2]- C_ss
            residual[3] = u[1][3]- k_ss
            residual[4] = u[1][7]- v_ss
            residual[5] = u[1][5]- init_ρ
            residual[6] = u[1][6]- i_ss
            residual[7] = u[end][7]- v_ss
            residual[8] = u[end][8]- v_ss
    end


    bvp1 = TwoPointBVProblem(NK!, bc1!, u_0(p), (0.0,params.T),(p))

    u = solve(bvp1, MIRK4(), dt=params.dt)
    
    function result(u,p)
        @unpack δ,s₀, S     = p
        @unpack q_ss, C_ss, k_ss, π_ss, i_ss, ℓ_ss, ι_ss, ρ_ss,v_ss,s_ss = SS(p)
                q  = @view u[1,:][:]
                C  = @view u[2,:][:]
                k  = @view u[3,:][:]
                π  = @view u[4,:][:]
                i  = @view u[6,:][:]
                v  = @view u[7,:][:]
                vˡ = @view u[8,:][:]

                n  = size(u)[1]

            sol1              = similar(zeros(n+3,size(u)[2]))
        sol1[1:n-3,:] = @view u[2:n-2,:]

        @unpack ι, ℓ, w, ν_k, χ, Y = static_funct(p)
        
        sol1[n-2,:] = ι.(q) .+ δ
        sol1[n-1,:] = ℓ.(C,k,q)
        sol1[n,:]   = Y.(C,k,q)
        sol1[n+1,:] = i.-π
        sol1[n+2,:] = v
        sol1[n+3,:] = s₀.+S*vˡ
        
        SS_vec = similar(sol1[:,1])
        SS_vec[1]     = C_ss
        SS_vec[2]     = k_ss
        SS_vec[3]     = π_ss
        SS_vec[4]     = ρ_ss
        SS_vec[5]     = i_ss
        SS_vec[6]     = ι_ss+ δ
        SS_vec[7]     = ℓ_ss
        SS_vec[8]     = Y.(C_ss,k_ss,q_ss)
        SS_vec[9]     = i_ss-π_ss
        SS_vec[10]    = v_ss
        SS_vec[11]    = s_ss


        return(SS_vec=SS_vec,sol1=sol1)
    end


    @unpack sol1, SS_vec = result(u,p)
    return (sol=sol1,SS=SS_vec,t=u.t)
end

@unpack T, dt = pp

function plot_IRF_quad_FTPL(;var =["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s"],
                        solution,T_end=T)
    N_end = T_end/dt+1
    val   = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s"]
    pos   = (zeros(length(var)))
    for k in 1: length(var)
        pos[k] = findfirst(isequal(var[k]),val)
    end
    pos = round.(Int, pos)
    val = val[pos]
    lab = [latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab = reshape(lab,(1,length(val)))

    SS  = solution.SS
    dev = ((solution.sol.-SS)./SS)*100

     
    pp = [dev[k,1:round(Int,N_end)] for k in pos]

    p=plot(solution.t[1:round(Int,N_end)],pp,
        label          = lab,
        xlabel         = L"t",
        legendfontsize = 8,
        ylabel         = L"\%",
        legend         = :outertopright)
    display(p)
    return(p)
end


function compute_dev_quad_FTPL(;solution,n,T)
        SS  = solution.SS[n]
        dev = (((@view solution.sol[n,:]).-SS)./SS)*100
        N   = T/dt+1
        cum = sum(@view dev[1:floor(Int,N)])
        return (cum)
end


function plot_θ_cum_quad_FTPL(;var="Y",θ_range=range(.1,500,length=10),
                T_range=[T],κ_range=[3,30,300],ind_Taylor=pp.ind_Taylor,ϕ_FTPL=pp.ϕ_FTPL)
    val = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s"]
    n   = findfirst(isequal(var), val)
    N   = length(T_range)*length(κ_range)
    lab = [latexstring("\$T={$(T)},\\kappa={$(κ)}\$") for (T,κ) in Iterators.product(T_range, κ_range)][:]
    lab = reshape(lab,1,N)

    y = similar(zeros(length(θ_range),N))
    j = 0
    for θ in θ_range
        j = j+1
        k = 0
        for κ in κ_range
        solution = solve_system_quad_FTPL(;params=define_env(θ=θ,κ=κ,T=T,N_t=T/dt,ϕ_FTPL=ϕ_FTPL,ind_Taylor=ind_Taylor))
            for T in T_range
              k    = k+1
            y[j,k] = compute_dev_quad_FTPL(;solution=solution,n=n,T=T)
            end
        end
    end
    lines = [:dash for k in 1:N]
    for   k in 1: floor(Int,N/2)
    lines[2*k] =:solid
    end
    lines = reshape(lines,1,N)
    p     = plot(θ_range,
            y, 
            label          = lab,
            xlabel         = L"\theta",
            ylabel         = latexstring("\$\\sum_{t=0}{T}\\widehat{{$(val[n])}}_{t}\\left(\\%,\\phi=$(ϕ_FTPL)\\right)\$"),
            legendfontsize = 7,
            palette        = palette([:blue,:red],N),
            linestyle      = lines,
            legend         = :outertopright)
    savefig(p,"theta_cum_$(val[n])_$(T_range[1])_quad.svg")
    display(p)
end

