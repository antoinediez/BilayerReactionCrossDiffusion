using StaticArrays

function RH_table(a)
    n = length(a) - 1
    np = n+1
    n2 = ceil(Int,(n+1)/2)
    T = zeros(Num,n+1,n2)
    
    for i in 1:np
        for j in 1:n2
            if i==1
                try
                    T[i,j] = a[np-2*(j-1)]
                catch e
                    println(e)
                end
            elseif i==2
                try
                    T[i,j] = a[np-1-2*(j-1)]
                catch e
                end
            else
                # if !(T[i-1,1]==0.0)
                    try
                        T[i,j] = -(T[i-2,1]*T[i-1,j+1] - T[i-1,1]*T[i-2,j+1])/T[i-1,1]
                    catch e
                    end
                # end
            end
        end
    end
    return T
end

function RH_stable(a)
    T = RH_table(a)
    T1 = T[:,1]
    return sum((T1[1] .* T1) .< 0) == 0
end

function RH5(a)
	a0_eval,a1_eval,a2_eval,a3_eval,a4_eval = a 

    C4 = Float64(a4_eval > 0.0)
	C3 = Float64(a3_eval > 0.0)
	C2 = Float64(a2_eval > 0.0)
	C1 = Float64(a1_eval > 0.0)
	C0 = Float64(a0_eval > 0.0)
	C5 = Float64(a4_eval*a3_eval > a2_eval)
	Δ4 = -a0_eval^2 + a0_eval*a2_eval*a3_eval - a2_eval^2*a1_eval - a4_eval^2*a1_eval^2 + 2.0*a0_eval*a1_eval*a4_eval + a1_eval*a2_eval*a3_eval*a4_eval - a3_eval^2*a0_eval*a4_eval
	C6 = Float64(Δ4 > 0.0)

	return C0*C1*C2*C3*C4*C5*C6
end

function compute_RH_fun(n)
    @variables a[1:n]
    RH_sb = RH_stable(a)
    RH_fun = eval(build_function(RH_sb,a;expression=Val{false}))
    return RH_fun
end

function unstable_modes(X_all,a_fun,var...;scale="log",RH_fun=nothing)
    if isnothing(RH_fun)
        RH_fun = compute_RH_fun(length(a))
    end
    stab = [RH_fun(a_fun(x,var...)) for x in X_all]
    nb = sum(1.0 .- stab)
    if scale=="log"
        return log(1.0 + nb)
    else
        return nb
    end
end



            