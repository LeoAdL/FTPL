function define_env(;ρ̄      =  2.2/100,
                     θᵨ     =   .22,
                     σᵨ     =   0.003,
                     θᵢ      =   1.0,
                     θ      =   100.0,
                     ϵ      =   11,
                     ψ      =   1/2.0,
                     ϕ      =   1.25,
                     γ      =   1.0,
                     T      =   50.0,
                     N_t    =   100.0)
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
                     dt     =   dt,
                     init_ρ=init_ρ)
    return params
end

function solve_system(;params)
    function NK_Rote!(du,u,p,t)
        @unpack σ,ϵ,θ,σ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   params
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
    sol1[1,:] = sol1[1,:].-1.0
    return (sol=sol1*100.0,SS=SS_vec,initial=init,t=sol1.t)
end

function plot_IRF(;pos =[1,2,3,4],solution)
    val =["x","\\pi","i","\\rho"]
    val =val[pos]
    lab=[latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab=reshape(lab,(1,length(val)))

    pp =[(solution.sol)'[:,k] for k in pos]
    p=plot(solution.t,pp,
        label=lab,
        xlabel=L"t", 
        ylabel=L"\%",
        legend=:topright)
    display(p)
    return(p)
end


function compute_dev(;θ,T)
        solution=solve_system(;params=define_env(θ=θ,T=T))
        solll     =   solution.sol[1,:]
        cum = sum(solll)
        return (impact=solll[1],cum=cum/T)
end

function plot_θ(;range=.1:19.9:100.0)
    lab=[latexstring("\$\\theta={$(θ)}\$") for θ in range']
    p=plot(compute_dev(10).t,
            [compute_dev(θ).dev for θ in range], 
            label=lab,
            xlabel=L"t", 
            ylabel=L"\%",
            legend=:topright,
            palette = palette([:blue,:red], length(range)))
    savefig(p,"theta.svg")
    display(p)
end
