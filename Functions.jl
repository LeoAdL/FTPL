function solve_system(;params)
    function NK_Rote!(du,u,p,t)
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   p
        x    =   u[1]

        π   =   u[2]

        i   =   u[3]

        ρ   =   u[4]

        du[1]  =   σ*(i-π-ρ)*x
        
        du[2]  =   -(ϵ-1.0)/θ*(x^(1.0\σ+ψ)-1.0)+π*((1.0-σ)*(i-π)+σ*ρ)

        du[3]  =   -θᵢ*(i-ϕ*π)

        du[4]  =   -θᵨ*(ρ-ρ̄)
        
    end

    function SS(p)
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   p

        ρ_ss    =   ρ̄
        π_ss    =   ρ̄/(ϕ-1.0)
        i_ss    =   ρ̄*ϕ/(ϕ-1.0)
        x_ss    =   (1.0+θ/(ϵ-1.0)*π_ss*(σ*ρ̄+(1-σ)*(i_ss-π_ss)))^(1.0/(1.0/σ+ψ))

        return (ρ_ss=ρ_ss,π_ss=π_ss,i_ss=i_ss,x_ss=x_ss)
    end
    @unpack T,dt,init_ρ,σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,γ,init_ρ = params

    p  =(σ=σ,ϵ=ϵ,θ=θ,ϕ=ϕ,ψ=ψ,ρ̄=ρ̄,θᵨ=θᵨ,θᵢ=θᵢ,γ=γ,init_ρ=init_ρ)


    function bc1!(residual,u,p,t)
        @unpack ρ_ss,π_ss,i_ss,x_ss= SS(p)
        @unpack init_ρ  =   p
        residual[1] =   u[end][1]- x_ss
        residual[2] =   u[end][2]- π_ss
        residual[3] =   u[1][3]- i_ss
        residual[4] =   u[1][4]- init_ρ
    end

    @unpack ρ_ss,π_ss,i_ss,x_ss= SS(p)
    @unpack T,dt,init_ρ = params

    SS_vec = [x_ss,π_ss,i_ss,ρ_ss]

    u0    =   [x_ss,π_ss,i_ss,init_ρ]
    tspan   =   (0.0,T)
    
    bvp1 = TwoPointBVProblem(NK_Rote!, bc1!, u0, tspan,(p))

    sol1 = solve(bvp1, MIRK4(), dt=dt)
    return (sol=sol1,SS=SS_vec,t=sol1.t)
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
