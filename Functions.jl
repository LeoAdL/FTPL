function solve_system(;params)
    function NK_Rote!(du,u,p,t)
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ =   p
        x    =   u[1]

        π   =   u[2]

        i   =   u[3]

        ρ   =   u[4]

        du[1]  =   σ*(i-π-ρ)*x
        
        du[2]  =   -(ϵ-1.0)/θ*(x^(1.0/σ+ψ)-1.0)+π*((1.0-σ)*(i-π)+σ*ρ)

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
    
    function u_0(p)
    @unpack π_ss,i_ss,x_ss = SS(p)
    @unpack init_ρ  =   p
    return ([x_ss,π_ss,i_ss,init_ρ])
    end


    p  =(σ=params.σ,ϵ=params.ϵ,θ=params.θ,
        ϕ=params.ϕ,ψ=params.ψ,ρ̄=params.ρ̄,θᵨ=params.θᵨ,θᵢ=params.θᵢ,
        γ=params.γ,init_ρ=params.init_ρ)


    function bc1!(residual,u,p,t)
        @unpack ρ_ss,π_ss,i_ss,x_ss= SS(p)
        @unpack init_ρ  =   p
        residual[1] =   u[end][1]- x_ss
        residual[2] =   u[end][2]- π_ss
        residual[3] =   u[1][3]- i_ss
        residual[4] =   u[1][4]- init_ρ
    end

    
    function SS_vec(p)
    @unpack ρ_ss,π_ss,i_ss,x_ss= SS(p)
    return([x_ss,π_ss,i_ss,ρ_ss])
    end

    @unpack T,dt,init_ρ = params
    bvp1 = TwoPointBVProblem(NK_Rote!, bc1!,u_0(p), (0.0,T),(p))
    sol1 = solve(bvp1, MIRK4(), dt=dt)
    return (sol=sol1,SS=SS_vec(p),t=sol1.t)
end

@unpack T,ϕ,dt   = pp

function plot_IRF(;var =["x","\\pi","i","\\rho"],
                        solution,T_end=T)
    N_end=T_end/dt+1
    val =["x","\\pi","i","\\rho"]
    pos=(zeros(length(var)))
    for k in 1:length(var)
        pos[k]=findfirst(isequal(var[k]),val)
    end
    pos=round.(Int, pos)
    val =val[pos]
    lab=[latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab=reshape(lab,(1,length(val)))

    SS  =   solution.SS
    dev =   ((solution.sol.-SS)./SS)*100

     
    pp = [dev[k,1:round(Int,N_end)] for k in pos]

    p=plot(solution.t[1:round(Int,N_end)],pp,
        label=lab,
        xlabel=L"t",
        legendfontsize=8,
        ylabel=L"\%",
        legend=:outertopright)
    display(p)
    return(p)
end


function compute_dev_quad(;solution,n,T)
        SS  =   solution.SS[n]
        dev =   (((@view solution.sol[n,:]).-SS)./SS)*100
        N=T/dt+1
        cum = sum(@view dev[1:floor(Int,N)])
        return (cum)
end

function plot_θ_cum(;var="x",θ_range=range(.1,500,length=10),ϕ=ϕ,
                T_range=[0,T])
    val =["x","\\pi","i","\\rho"]
    n=findfirst(isequal(var), val)
    N   = length(T_range)  
    lab=[latexstring("\$T={$(T)}\$") for T in T_range]
    lab=reshape(lab,1,N)

    y = similar(zeros(length(θ_range),N))
    j=0
    for θ in θ_range
        j=j+1
        k=0
        solution=solve_system(;params=define_env(θ=θ,T=T,N_t=T/dt,ϕ=ϕ))
            for T in T_range
            k = k+1
            y[j,k] = compute_dev_quad(;solution=solution,n=n,T=T)
        end
    end
    p=plot(θ_range,
            y, 
            label=lab,
            xlabel=L"\theta", 
            ylabel=latexstring("\$\\sum_{t=0}{T}\\widehat{{$(val[n])}}_{t}\\left(\\%,\\phi=$(ϕ)\\right)\$"),
            legendfontsize=7,
            palette = palette([:blue,:red],N),
            legend=:outertopright)
    savefig(p,"theta_cum_$(val[n])_$(T_range[1]).svg")
    display(p)
end

