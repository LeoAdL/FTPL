function solve_system_quad(;params)
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
            @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ =   params
            F   =   zeros(6)
            π=u[1]
            C=u[2]
            q=u[3]
            k=u[4]
            ρ=u[5]
            i=u[6]

            F[1]=δ+(i-π)+(ι(q)+κ*(ι(q))/2)/q-ι(q)
            F[1]=F[1]-ν_k(q,w(ℓ(C,k,ι),C),ℓ(C,k,ι),k)

            F[2]=i-π-ρ
            
            F[3]=ι(q)-δ-(ι(q)+κ*(ι(q))/2)/q

            F[4]=(1-σ)*(i-π)+σ*ρ

            F[5]=ρ̄-ρ

            F[6]=i-ϕ*π
            return(F)
        end
        u = nlsolve(obj, zeros(6))

        π=u[1]
        C=u[2]
        q=u[3]
        k=u[4]
        ρ=u[5]
        i=u[6]

        ℓ=ℓ(C,k,ι)
        w=w(ℓ,C)
        ν=ν_k(q,w,ℓ,k)
        χₙ=χ(w,q,ν_k)

        return(π=π,
                C=C,
                q=q,
                k=k,
                ρ=ρ,
                i=i,
                χₙ=χₙ)
    end

    function NK!(du,u,p,t)
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ =   params
        χₙ=NK_natural().χₙ
            q=u[1]
            C=u[2]
            k=u[3]
            π=u[4]
            ρ=u[5]
            i=u[6]

            du[1]=δ+(i-π)+(ι(q)+κ*(ι(q))/2)/q-ι(q)
            du[1]=du[1]-ν_k(q,w(ℓ(C,k,ι),C),ℓ(C,k,ι),k)
            du[1]=du[1]*q

            du[2]=σ*C*(i-π-ρ)
            
            du[3]=(ι(q)-δ-(ι(q)+κ*(ι(q))/2)/q)*k

            du[4]=+((1-σ)*(i-π)+σ*ρ)*π

            du[4]=du[4]-(ϵ-1)/(θ)*(χ(w(ℓ(C,k,ι),C),q,ν_k(q,w(ℓ(C,k,ι),C),ℓ(C,k,ι),k))/χₙ-1)

            du[5]=θᵨ*(ρ̄-ρ)

            du[6]=-θᵢ*(i-ϕ*π)
        
    end

    function SS()
    @unpack T = params
    init    =   zeros(6)
    tspan   =   (0.0,T)
    SS=ODEProblem(NK!, init, tspan)
    u = solve(SS,SSRootfind())
            q_ss=u[1]
            C_ss=u[2]
            k_ss=u[3]
            π_ss=u[4]
            ρ_ss=u[5]
            i_ss=u[6]
       return(π_ss=π_ss,
                C_ss=C_ss,
                q_ss=q_ss,
                k_ss=k_ss,
                ρ_ss=ρ_ss,
                i_ss=i_ss)    
    end

    function bc1!(residual,u,p,t)
        @unpack q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss= SS()
        @unpack init_ρ  =   params
        residual[1] =   u[end][1]- q_ss
        residual[2] =   u[end][2]- C_ss
        residual[3] =   u[end][3]- k_ss
        residual[4] =   u[end][4]- π_ss
        residual[5] =   u[1][5]- init_ρ
        residual[6] =   u[1][6]- i_ss
    end

    @unpack q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss= SS()
    @unpack T,dt,init_ρ = params

    SS_vec = [q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss]

    init    =   [q_ss,C_ss,k_ss,π_ss,init_ρ,i_ss]
    tspan   =   (0.0,T)


    bvp1 = TwoPointBVProblem(NK!, bc1!, init, tspan)

    sol1 = solve(bvp1, MIRK4(), dt=dt)
    return (sol=sol1,SS=SS_vec,initial=init,t=sol1.t)
end

function plot_IRF_quad(;pos =[1,2,3,4,5,6],solution)
    val =["q","Y","k","\\pirho","i",]
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


function compute_dev_quad(;θ,T)
        solution=solve_system(;params=define_env(θ=θ))
        SS  =   solution.SS[2]
        dev =   ((solution.sol[1,1:end].-SS)./SS)*100
        cum = sum(dev[1:floor(Int,T)])
        return (impact=dev[2],cum=cum)
end

function plot_θ_impact_quad(;theta_range=range(10^(-3),500,length=10))
    p=plot(theta_range,
            [compute_dev(;θ=θ,T=T).impact for θ in theta_range], 
            xlabel=L"\theta", 
            ylabel=L"\%",
            legend=:bottomright,
            palette = palette([:blue,:red], length(theta_range)))
    savefig(p,"theta.svg")
    display(p)
end

function plot_θ_cum_quad(;theta_range=range(1,500,length=10),
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
