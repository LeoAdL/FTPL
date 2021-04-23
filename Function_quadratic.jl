function solve_system_quad(;params)
    function ι(q)
        @unpack κ=params
        return((q-1.0)/κ)
    end

    function ℓ(C,k,q)
        @unpack κ,α,A=params
        return(k*((C/k+ι(q)+κ*(ι(q))^(2.0)/2.0)/A)^(1.0/(1.0-α)))
    end

    function w(C,k,q)
        @unpack γ,ψ =   params
        return (ℓ(C,k,q)^(ψ)*C^(γ))
    end


    function ν_k(C,k,q)
        @unpack α=params
        return((1.0/q)*(α/(1.0-α))*w(C,k,q)*(ℓ(C,k,q)/k))
    end
        
    function χ(C,k,q)
        @unpack α=params
        return((w(C,k,q)/(1.0-α))^(1.0-α)*(q*ν_k(C,k,q)/α)^(α))
    end

    function Y(C,k,q)
        @unpack α,A=params
        return(A*(k)^(α)*(ℓ(C,k,q))^(1-α))
    end


    function NK!(du,u,p,t)
        @unpack σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ,A =   params
            χₙ=A*(ϵ-1.0)/ϵ
            q=u[1]
            C=u[2]
            k=u[3]
            π=u[4]
            ρ=u[5]
            i=u[6]

            du[1]=q*(i-π+(ι(q)+κ*(ι(q))^(2.0)/2.0)/q-ι(q)-ν_k(C,k,q))
            
            du[2]=σ*C*(i-π-ρ)
            
            du[3]=ι(q)*k

            du[4]=π*((1.0-σ)*(i-π)+σ*ρ)-((ϵ-1)/θ)*(χ(C,k,q)/χₙ-1.0)

            du[5]=-θᵨ*(ρ-ρ̄)

            du[6]=-θᵢ*(i-ϕ*π)
        
    end

    function SS()
        @unpack α,γ,σ,ϵ,θ,ϕ,ψ,ρ̄,θᵨ,θᵢ,κ,δ,A =   params

                    k_c=(α/(ρ̄))*((ϵ-1.0)/ϵ)*(1.0+θ/(ϵ-1.0)*ρ̄^(2)/(ϕ-1.0))
    
                    k_l=(A*k_c)^(1.0/(1.0-α))

            q_ss=1.0
            k_ss=(ρ̄*((1-α)/α)*(k_l)^(1+ψ)*(k_c)^(γ))^(1/(ψ+γ))
            C_ss=k_ss/k_c
            π_ss=ρ̄/(ϕ-1.0)
            ρ_ss=ρ̄
            i_ss=ϕ*π_ss
       return(π_ss=π_ss,
                C_ss=C_ss,
                q_ss=q_ss,
                k_ss=k_ss,
                ρ_ss=ρ_ss,
                i_ss=i_ss,
                ℓ_ss=k_ss/k_l,
                ι_ss=0.0)    
    end

    function bc1!(residual,u,p,t)
        @unpack q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss= SS()
        @unpack init_ρ  =   params
        residual[1] =   u[end][1]- q_ss
        residual[2] =   u[end][2]- C_ss
        residual[3] =   u[1][3]- k_ss
        residual[4] =   u[end][4]- π_ss
        residual[5] =   u[1][5]- init_ρ
        residual[6] =   u[1][6]- i_ss
    end

    @unpack q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss,ℓ_ss,ι_ss= SS()
    @unpack T,dt,init_ρ,α = params

    SS_vec = [q_ss,C_ss,k_ss,π_ss,ρ_ss,i_ss]

    u0    =   [q_ss,C_ss,k_ss,π_ss,init_ρ,i_ss]
    tspan   =   (0.0,T)

    bvp1 = TwoPointBVProblem(NK!, bc1!, u0, tspan)

    u = solve(bvp1, MIRK4(), dt=dt)

    q=u[1,:][:]
    C=u[2,:][:]
    k=u[3,:][:]

    sol1  =   similar(zeros(size(u)[1]+2,size(u)[2]))
    sol1[1:5,:]  = u[2:end,:]

    @unpack δ,A =   params

    sol1[6,:] =   ι.(q) .+ δ
    sol1[7,:] =   ℓ.(C,k,q)
    sol1[8,:] =   Y.(C,k,q)

    SS_vec = [C_ss,k_ss,π_ss,init_ρ,i_ss,ι_ss+ δ,
                ℓ_ss,
                A*(ℓ_ss)^(1.0-α)*(k_ss)^(α)]
    return (sol=sol1,SS=SS_vec,initial=init,t=u.t)
end

function plot_IRF_quad(;pos =[1,2,3,4,5,6,7,8],solution)
    val =["C","k","\\pi","\\rho","i","\\iota","\\ell","Y"]
    val =val[pos]
    lab=[latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab=reshape(lab,(1,length(val)))

    SS  =   solution.SS
    dev =   ((solution.sol.-SS)./SS)*100

     
    pp = [dev[k,:] for k in pos]

    p=plot(solution.t,pp,
        label=lab,
        xlabel=L"t", 
        ylabel=L"\%",
        legend=:topright)
    display(p)
    return(p)
end


function compute_dev_quad(;θ,T,κ)
        solution=solve_system_quad(;params=define_env(θ=θ,κ=κ,N_t=50))
        SS  =   solution.SS[end]
        dev =   ((solution.sol[end,:].-SS)./SS)*100
        cum = sum(dev[1:floor(Int,T)])
        return (impact=dev[1],cum=cum)
end

@unpack T = define_env()

function plot_θ_cum_quad(;theta_range=range(.1,500,length=10),
                T_range=[T],κ_range=[.5,5,50,500.0,5*10.0^6.0])
    N   = length(T_range)*length(κ_range)    
    lab=[latexstring("\$T={$(T)},\\kappa={$(κ)}\$") for (T,κ) in Iterators.product(T_range, κ_range)][:]
    lab=reshape(lab,1,N)
    y = [[compute_dev_quad(;θ=θ,T=T,κ=κ).cum for θ in theta_range] for (T,κ) in Iterators.product(T_range, κ_range)]
    y = reshape(y,1,N)[:]
    y_pp= y[1]

    for n in 2:N
        y_pp=hcat(y_pp,y[n])
    end
    p=plot(theta_range,
            y_pp, 
            label=lab,
            xlabel=L"\theta", 
            ylabel=L"\sum_{t=1}{T}\widehat{{Y}}_{t}(\%)",
            legendfontsize=7,
            palette = palette([:blue,:red],N),
            legend=:bottomright)
    savefig(p,"theta_cum_$(T_range[1]).svg")
    display(p)
end
