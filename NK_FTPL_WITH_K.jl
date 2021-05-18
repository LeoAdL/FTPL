function solve_system_quad_FTPL(;params)

    function static_funct(p)
        @unpack σ, ψ, A, α, γ,κ = p
        function ι(q)
            return((q-1.0)/κ)
        end

        function ℓ(C,k,q)
            return(k*((C/k+ι(q))/A)^(1.0/(1.0-α)))
        end

        function w(C,k,q)
            return (ℓ(C,k,q)^(ψ)*C^(γ))
        end


        function ν_k(C,k,q)
            return((1.0/q)*(α/(1.0-α))*w(C,k,q)*(ℓ(C,k,q)/k))
        end
            
        function χ(C,k,q)
            return((w(C,k,q)/(1.0-α))^(1.0-α)*(q*ν_k(C,k,q)/α)^(α))
        end

        function Y(C,k,q)
            return(A*(k)^(α)*(ℓ(C,k,q))^(1-α))
        end
        return (ι=ι,ℓ=ℓ,w=w,ν_k=ν_k,χ=χ,Y=Y)
    end





    function NK!(du,u,p,t)
        @unpack σ, ϵ, θ, ψ, ρ̄, θᵨ, θᵢ, κ, δ, A, χₙ, S, ind_Taylor, i_target, ϕ_FTPL,long_term = p
         @unpack ι, ℓ, w, ν_k, χ,  Y  = static_funct(p)
                 q  = u[1]
                 C  = u[2]
                 k  = u[3]
                 π  = u[4]
                 ρ  = u[5]
                 i  = u[6]
                 v  = u[7]
                 s  = u[8]
                 Q  = u[9]

                 y=1.0/Q



            du[1] = q*(i-π+(ι(q)+κ*(ι(q))^(2.0)/2.0)/q-ι(q)-ν_k(C,k,q))
            
            du[2] = σ*C*(i-π-ρ)
            
            du[3] = ι(q)*k

            du[4] = π*((1.0-σ)*(i-π)+σ*ρ)-((ϵ-1.0)/θ)*(χ(C,k,q)/χₙ-1.0)
            
            du[5] = -θᵨ*(ρ-ρ̄)

            du[6] = -θᵢ*(i-ϕ_FTPL*π-ρ̄)*ind_Taylor

            du[7] = v*(i-π) -s 

            du[8] = S*du[7]

            du[9] = Q*(i-y)*long_term
        
    end

    function SS(p)
        @unpack α, γ, σ, ϵ, θ, ψ, ρ̄, θᵨ, θᵢ, κ, δ, A, S, ind_Taylor, i_target, ϕ_FTPL, s₀ = p
        @unpack Y  = static_funct(p)

            q_ss = 1.0
            ρ_ss = ρ̄
            i_ss = ρ̄*ind_Taylor +(1-ind_Taylor)*i_target
            π_ss = i_ss - ρ_ss
            Q_ss = 1/i_ss
            k_c  = (α/(ρ̄))*((ϵ-1.0)/ϵ)*(1.0+θ/(ϵ-1.0)*ρ̄*π_ss)
            k_l  = (A*k_c)^(1.0/(1.0-α))
            k_ss = (ρ̄*((1-α)/α)*(k_l)^(1+ψ)*(k_c)^(γ))^(1/(ψ+γ))
            C_ss = k_ss/k_c

            Y_ss = Y(C_ss,k_ss,q_ss)
            v_ss = s₀*Y_ss/(i_ss-π_ss-S)
            s_ss = s₀*Y_ss + S*v_ss



       return(π_ss=π_ss,
                C_ss = C_ss,
                q_ss = q_ss,
                k_ss = k_ss,
                ρ_ss = ρ_ss,
                i_ss = i_ss,
                ℓ_ss = k_ss/k_l,
                ι_ss = 0.0,
                v_ss = v_ss,
                s_ss = s_ss,
                Q_ss = Q_ss)
    end

    function u_0(p)
        @unpack q_ss,   C_ss, k_ss, π_ss, i_ss,v_ss,s_ss,Q_ss = SS(p)
        @unpack init_ρ = p
    return ([q_ss,C_ss,k_ss,π_ss,init_ρ,i_ss,v_ss,s_ss,Q_ss])
    end

    p  =    (σ=params.σ,
            ϵ          = params.ϵ,
            θ          = params.θ,
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
            s₀         = params.s₀,
            long_term  = params.long_term)



    function bc1!(residual,u,p,t)
            @unpack  q_ss,   C_ss, k_ss, π_ss, ρ_ss, i_ss,v_ss,s_ss,Q_ss = SS(p)
            @unpack  init_ρ,long_term = p
            residual[1]     = u[end][1]- q_ss
            residual[2]     = u[end][2]- C_ss
            residual[3]     = u[end][7]- v_ss
            residual[4]     = u[end][8]- s_ss
            residual[5]     = u[1][3]- k_ss
            residual[6]     = u[1][5]- init_ρ
            residual[7]     = u[1][6]- i_ss
            residual[8]     = u[1][7]- v_ss*(1+(u[1][9]-Q_ss)/Q_ss)
            residual[9]     = u[end][9]- Q_ss
    end

    
    bvp1 = TwoPointBVProblem(NK!, bc1!, u_0(p), (0.0,params.T),(p))

    u = solve(bvp1, MIRK4(), dt=params.dt)
    
    function result(u,p)
        @unpack δ    = p
        @unpack q_ss, C_ss, k_ss, π_ss, i_ss, ℓ_ss, ι_ss, ρ_ss,v_ss,s_ss,Q_ss = SS(p)
                q  = @view u[1,:][:]
                C  = @view u[2,:][:]
                k  = @view u[3,:][:]
                π  = @view u[4,:][:]
                i  = @view u[6,:][:]
                v  = @view u[7,:][:]
                s  = @view u[8,:][:]
                Q  = @view u[9,:][:]

                n  = size(u)[1]

            sol1              = similar(zeros(n+3,size(u)[2]))
        sol1[1:n-3,:] = @view u[2:n-2,:]

        @unpack ι, ℓ, w, ν_k, χ, Y = static_funct(p)
        sol1[3,:]   = sol1[3,:].+1
        sol1[n-3,:] = ι.(q) .+ δ
        sol1[n-2,:] = ℓ.(C,k,q)
        sol1[n-1,:] = Y.(C,k,q)
        sol1[n,:]   = i.-π
        sol1[n+1,:] = v
        sol1[n+2,:] = s
        sol1[n+3,:] = 1.0./Q
        
        SS_vec = similar(sol1[:,1])
        SS_vec[1]     = C_ss
        SS_vec[2]     = k_ss
        SS_vec[3]     = π_ss +1.0
        SS_vec[4]     = ρ_ss
        SS_vec[5]     = i_ss
        SS_vec[6]     = ι_ss+ δ
        SS_vec[7]     = ℓ_ss
        SS_vec[8]     = Y.(C_ss,k_ss,q_ss)
        SS_vec[9]     = i_ss-π_ss
        SS_vec[10]    = v_ss
        SS_vec[11]    = s_ss
        SS_vec[12]    = 1.0/Q_ss


        return(SS_vec=SS_vec,sol1=sol1)
    end


    @unpack sol1, SS_vec = result(u,p)
    return (sol=sol1,SS=SS_vec,t=u.t)
end

@unpack T, dt = pp

function plot_IRF_quad_FTPL(;var =["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s","y"],
                        solution,T_end=T)
    N_end = T_end/dt+1
    val   = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s","y"]
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


function compute_dev_quad_FTPL(;solution,n,T)
        SS  = solution.SS[n]
        dev = (((@view solution.sol[n,:]).-SS)./SS)*100
        N   = T/dt+1
        cum = sum(@view dev[1:floor(Int,N)])
        return (cum)
end


function plot_θ_cum_quad_FTPL(;var="Y",θ_range=range(.1,500,length=5),
                T_range=[0,T],κ_range=[3,30,300],ind_Taylor=pp.ind_Taylor,ϕ_FTPL=pp.ϕ_FTPL,params=pp,long_term=pp.long_term)
    @unpack T, dt = pp
    val = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s","y"]
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
        solution = solve_system_quad_FTPL(;params=define_env(θ=θ,κ=κ,T=T,N_t=T/dt,ϕ_FTPL=ϕ_FTPL,ind_Taylor=ind_Taylor,long_term=long_term))
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
    savefig(p,"theta_cum_$(ϕ_FTPL)_quad_FTPL_K.svg")
    display(p)
end

function plot_all_longterm(;var =["C","k","\\pi","i","\\iota","\\ell","Y","r","v","s","y"],pp=pp)
    

    val   = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s","y"]
    pos   = (zeros(length(var)))
    for k in 1: length(var)
        pos[k] = findfirst(isequal(var[k]),val)
    end
    pos = round.(Int, pos)
    val = val[pos]
    lab = [latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab = reshape(lab,(1,length(val)))



    p=plot(layout = length(var),title= lab,size = (1800,1000),palette= :Dark2_4)

    style=[:dash, :dot, :dash, :dot]
    j=0
    for Taylor in [0.0,1.0]
        for l in [0.0,1.0]
            j=j+1
            param=define_env(T=pp.T,N_t=pp.T/pp.dt,ind_Taylor=1.0,ϕ_FTPL=pp.ϕ_FTPL*Taylor,S=0.0
                ,long_term=l)            
            solution =solve_system_quad_FTPL(params=param)
            SS  = solution.SS
            dev = ((solution.sol.-SS)./SS)*100 
            plot = [dev[k,:] for k in pos]
            label=latexstring("\$LD=$(param.long_term),\\phi=$(param.ϕ_FTPL)\$")
            plot!(solution.t,plot,
            xlabel         = L"t",
            legendfontsize = 7,
            ylabel         = L"\%",
            label          = label,
            legend         = :outertopright,
            linestyle      = style[j])
        end
    end
    display(p)
    return(p)
    savefig(p,"long_term_FTPL_with_K.svg")
end

function plot_all_S(;var =["C","k","\\pi","i","\\iota","\\ell","Y","r","v","s","y"],pp=pp)
    

    val   = ["C","k","\\pi","\\rho","i","\\iota","\\ell","Y","r","v","s","y"]
    pos   = (zeros(length(var)))
    for k in 1: length(var)
        pos[k] = findfirst(isequal(var[k]),val)
    end
    pos = round.(Int, pos)
    val = val[pos]
    lab = [latexstring("\$\\widehat{{$(u)}}_{t}\$") for u in val]
    lab = reshape(lab,(1,length(val)))


    p=plot(layout = length(var),title= lab,size = (1800,1000),palette= :Dark2_4)

    style=[:dash, :dot, :dash, :dot]
    j=0
    for Taylor in [0.0,1.0]
        for S in [0.0,1.0]
            j=j+1
            param=define_env(T=pp.T,N_t=pp.T/pp.dt,ind_Taylor=Taylor,ϕ_FTPL=pp.ϕ_FTPL*Taylor,S=S
                ,long_term=0.0)            
            solution =solve_system_quad_FTPL(params=param)
            SS  = solution.SS
            dev = ((solution.sol.-SS)./SS)*100 
            plot = [dev[k,:] for k in pos]
            label=latexstring("\$S=$(param.S),\\phi=$(param.ϕ_FTPL)\$")
            plot!(solution.t,plot,
            xlabel         = L"t",
            legendfontsize = 7,
            ylabel         = L"\%",
            label          = label,
            legend         = :outertopright,
            linestyle      = style[j])
        end
    end
    display(p)
    return(p)
    savefig(p,"S_shape_with_K.svg")
end