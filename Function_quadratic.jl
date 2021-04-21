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
                     κ      =   0.1,
                     δ      =   .1)
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
                     dt     =   dt,
                     init_ρ=init_ρ)
    return params
end

function solve_system(;params)
    function w(ℓ,C)
        @unpack σ,ψ =   params
        return (ℓ^(ψ)*C^(1/σ))
    end

    function ι(q)
        @unpack κ=params
        return((q-1)/k)
    end

    function ℓ(C,k,ι)
        @unpack κ,α=params
        dk  =   ι+κ*(ι)^(2)/2
        return(k*(C/k+dk)^(1/(1-α)))
    end

    function ν_k(q,w,ℓ,k)
        @unpack α=params
        return((1/q)*(α/(1-α))*w*(ℓ/k))
    end
        
    function χ(w,q,ν_k)
        @unpack α=params
        return((w/(1-α))^(1-α)*(q*ν_k/α)^(α))
    end
    function NK_natural()
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ =   params
        function obj(u)
            F   =   zeros(6)
            π=u[1]
            C=u[2]
            q=u[3]
            k=u[4]
            ρ=u[5]
            i=u[6]

            F[1]=δ+(i-π)+(ι(q)+κ*(ι(q))/2)/q-ι(q)
            F[1]=F[1]-ν_k(q,w(ℓ,C),ℓ(C,k,ι),k)

            F[2]=i-π-ρ
            
            F[3]=ι(q)-δ-(ι(q)+κ*(ι(q))/2)/q

            F[4]=(1-σ)*(i-π)+σ*ρ

            F[5]=ρ̄-ρ

            F[6]=i-ϕ*π
        end
    end

    function NK!(du,u,p,t)
        @unpack σ,ϵ,θ,σ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ =   params
        x    =   u[1]

        π   =   u[2]

        i   =   u[3]

        ρ   =   u[4]

        du[1]  =   σ*(i-π-ρ)*x
        
        du[2]  =   -(ϵ-1.0)/θ*(x^(1.0\σ+ψ)-1.0)+π*((1.0-σ)*(i-π)+σ*ρ)

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

    SS_vec = [x_ss,π_ss,i_ss,ρ_ss]

    init    =   [x_ss,π_ss,i_ss,init_ρ]
    tspan   =   (0.0,T)


    bvp1 = TwoPointBVProblem(NK_Rote!, bc1!, init, tspan)

    sol1 = solve(bvp1, MIRK4(), dt=dt)
    return (sol=sol1,SS=SS_vec,initial=init,t=sol1.t)
end

function plot_IRF(;pos =[1,2,3,4],solution)
    val =["x","\\pi","i","\\rho"]
    val =val[pos]
    lab=[latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab=reshape(lab,(1,length(val)))

    SS  =   solution.SS
    dev =   ((solution.sol[1:end].-SS)./SS)*100

     
    pp = [dev'[:,k] for k in pos]

    p=plot(solution.t,pp,
        label=lab,
        xlabel=L"t", 
        ylabel=L"\%",
        legend=:topright)
    display(p)
    return(p)
end


function compute_dev(;θ,T)
        solution=solve_system(;params=define_env(θ=θ))
        SS  =   solution.SS[1]
        dev =   ((solution.sol[1,1:end].-SS)./SS)*100
        cum = sum(dev[1:floor(Int,T)])
        return (impact=dev[1],cum=cum)
end

function plot_θ_impact(;theta_range=range(10^(-3),500,length=10))
    p=plot(theta_range,
            [compute_dev(;θ=θ,T=T).impact for θ in theta_range], 
            xlabel=L"\theta", 
            ylabel=L"\%",
            legend=:bottomright,
            palette = palette([:blue,:red], length(theta_range)))
    savefig(p,"theta.svg")
    display(p)
end

function plot_θ_cum(;theta_range=range(1,500,length=10),
                T_range=[1,2,4,10,20,50])
    p=plot(theta_range,
            [[compute_dev(;θ=θ,T=T).cum for θ in theta_range] for T in T_range], 
            label=[latexstring("\$T={$(T)}\$") for T in T_range'],
            xlabel=L"\theta", 
            ylabel=L"\sum_{t=1}^{T}\widehat{{x}}_{t}(\%)",
            legend=:topleft,
            palette = palette([:blue,:red], length(T_range)))
    savefig(p,"theta_cum.svg")
    display(p)
end
