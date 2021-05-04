function solve_system(;params)

    function NK_FTPL!(du,u,p,t)
        @unpack σ, ϵ, θ, ϕ, ψ, ρ̄, θᵨ, θᵢ, S, ϕ_FTPL, ind_Taylor, long_term = p
        
        x = u[1]

        π = u[2]

        i = u[3]

        ρ = u[4]

        v = u[5]

        s = u[6]

        Q = u[7]

        y = 1/Q

        du[1] = σ*(i-π-ρ)*x
        
        du[2] = -(ϵ-1.0)/θ*(x^(1.0/σ+ψ)-1.0)+π*((1.0-σ)*(i-π)+σ*ρ)

        du[3] = -θᵢ*(i-ϕ_FTPL*π-ρ̄)*ind_Taylor

        du[4] = -θᵨ*(ρ-ρ̄)

        du[5] = v*(i-π) -s 

        du[6] = S*du[5]

        du[7] = Q*(i-y)*long_term
    end

    function SS(p)
        @unpack σ, ϵ, θ, ϕ, ψ, ρ̄, θᵨ, θᵢ, A, S, γ, ind_Taylor, i_target, ϕ_FTPL, s₀ = p

        Yₙ = A^(1+ψ/(ψ+γ))*((ϵ-1)/ϵ)^(1/(ψ+γ))

        ρ_ss = ρ̄
        i_ss = ρ̄*ind_Taylor +(1-ind_Taylor)*i_target
        π_ss = i_ss - ρ_ss
        Q_ss = 1/i_ss
        x_ss = (1.0+θ/(ϵ-1.0)*π_ss*(σ*ρ̄+(1-σ)*(i_ss-π_ss)))^(1.0/(1.0/σ+ψ))
        Y_ss = x_ss*Yₙ
        v_ss = s₀*Y_ss/(i_ss-π_ss-S)
        s_ss = s₀*Y_ss + S*v_ss
        return (ρ_ss=ρ_ss,π_ss=π_ss,i_ss=i_ss,
                x_ss=x_ss,v_ss=v_ss,Yₙ=Yₙ,s_ss=s_ss,Q_ss=Q_ss)
    end


    function u_0(p)
    @unpack π_ss,   i_ss, x_ss, v_ss,s_ss,Q_ss = SS(p)
    @unpack init_ρ = p
    return ([x_ss,π_ss,i_ss,init_ρ,v_ss,s_ss,Q_ss])
    end


    p  =(σ=params.σ,
        ϵ          = params.ϵ,
        θ          = params.θ,
        ϕ          = params.ϕ,
        ψ          = params.ψ,
        ρ̄         = params.ρ̄,
        θᵨ         = params.θᵨ,
        θᵢ         = params.θᵢ,
        γ          = params.γ,
        init_ρ     = params.init_ρ,
        S          = params.S,
        A          = params.A,
        i_target   = params.i_target,
        ind_Taylor = params.ind_Taylor,
        ϕ_FTPL     = params.ϕ_FTPL,
        s₀         = params.s₀,
        long_term  = params.long_term)


    function bc1!(residual,u,p,t)
        @unpack  ρ_ss,   π_ss, i_ss, x_ss, v_ss,s_ss,Q_ss = SS(p)
        @unpack  init_ρ = p
        residual[1]     = u[end][1]- x_ss
        residual[2]     = u[end][5]- v_ss
        residual[3]     = u[end][6]- s_ss
        residual[4]     = u[1][3]- i_ss
        residual[5]     = u[1][4]- init_ρ
        residual[6]     = u[1][5]- v_ss*(1+(u[1][7]-Q_ss)/Q_ss)
        residual[7]     = u[end][7]- Q_ss
    end

    function SS_vec(p)
    @unpack ρ_ss, π_ss, i_ss, x_ss, v_ss, s_ss,Q_ss = SS(p)
    return([x_ss,π_ss,i_ss,ρ_ss,v_ss,s_ss,Q_ss])
    end

    @unpack T,    dt, init_ρ = params
            bvp1 = TwoPointBVProblem(NK_FTPL!, bc1!,u_0(p), (0.0,T),(p))
            u    = solve(bvp1, MIRK4(), dt=dt)
    
    sol=similar(u)
    sol[1:size(u)[1],:]=@view u[1:size(u)[1],:]
    sol[end,:]=1.0./u[end,:]

    sol[2,:]=sol[2,:] .+1.0

    Vec_ss=SS_vec(p)
    SS=similar(Vec_ss)
    SS[1:size(u)[1],:]=@view Vec_ss[1:size(u)[1],:]
    SS[end,:]=1.0./Vec_ss[end,:]
    SS[2]   =   SS[2] +1.0
    return (sol=sol,SS=SS,t=u.t)
end

@unpack T, ϕ, dt = pp

function plot_IRF_FTPL(;var =["x","\\pi","i","\\rho","v","s","y"],
                        solution,T_end=T)
    N_end = T_end/dt+1
    val   = ["x","\\pi","i","\\rho","v","s","y"]
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
        legend         = :outertopright,
        palette        = :tab20)
    display(p)
    return(p)
end


function compute_dev_FTPL(;solution,n,T)
        SS  = solution.SS[n]
        dev = (((@view solution.sol[n,:]).-SS)./SS)*100
        N   = T/dt+1
        cum = sum(@view dev[1:floor(Int,N)])
        return (cum)
end

function plot_θ_cum_FTPL(;var="x",θ_range=range(1,500,length=20),ϕ=ϕ,
                T_range=[0,T],ind_Taylor=pp.ind_Taylor,ϕ_FTPL=pp.ϕ_FTPL,
                long_term=pp.long_term,T=pp.T,dt=pp.dt)
    val = ["x","\\pi","i","\\rho","v","s","y"]
    n   = findfirst(isequal(var), val)
    N   = length(T_range)
    lab = [latexstring("\$T={$(T)}\$") for T in T_range]
    lab = reshape(lab,1,N)

    y = similar(zeros(length(θ_range),N))
    j = 0
    for θ in θ_range
        j        = j+1
        k        = 0
        solution = solve_system(;params=define_env(θ=θ,T=T,N_t=T/dt,ϕ_FTPL=ϕ_FTPL,ind_Taylor=ind_Taylor,long_term=long_term))
            for T in T_range
              k    = k+1
            y[j,k] = compute_dev_FTPL(;solution=solution,n=n,T=T)
        end
    end
    p=plot(θ_range,
            y, 
            label          = lab,
            xlabel         = L"\theta",
            ylabel         = latexstring("\$\\sum_{t=0}{T}\\widehat{{$(val[n])}}_{t}\\left(\\%,\\phi=$(ϕ_FTPL*ind_Taylor)\\right)\$"),
            legendfontsize = 7,
            palette        = palette([:blue,:red],N),
            legend         = :outertopright)
    savefig(p,"theta_cum_$(ϕ)_FTPL.svg")
    display(p)
end

