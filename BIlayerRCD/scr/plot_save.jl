using CairoMakie
using ProgressMeter

function init_plot_2D(
    c,uB,vB,uS,vS;
    dx,dy,
    save_video=false,
    time=0.0,
    range_uS=nothing,
    range_vS=nothing,
    range_uB=nothing,
    range_vB=nothing,
    range_c=nothing,
    colormap=:greys,
    fontsize_label=16,
    fontsize_title=24
)
    round_time = nice_float2string(time,2)
    Nx = size(c)[1]
    Lx = dx * Nx
    Ny = size(c)[2]
    Ly = dy * Ny    

    xx = dx .* collect(0:(Nx-1)) .+ dx/2
    yy = dy .* collect(0:(Ny-1)) .+ dy/2

    fig = Figure(resolution=(800,1000),fontsize=fontsize_label)

    ax_uS = Axis(fig[1, 1], title="uS" * "\ntime=$round_time",titlesize=fontsize_title)
    ax_vS = Axis(fig[1, 3], title="vS"*"\ntime=$round_time",titlesize=fontsize_title)

    xlims!(ax_uS,0,Lx)
    @views hm_uS = lines!(ax_uS,xx,uS,color=:black)
    if ~isnothing(range_uS)
        ylims!(ax_uS,range_uS[1],range_uS[2])
    end

    xlims!(ax_vS,0,Lx)
    @views hm_vS = lines!(ax_vS,xx,vS,color=:black)
    if ~isnothing(range_vS)
        ylims!(ax_vS,range_vS[1],range_vS[2])
    end

    ax_uB = Axis(fig[2, 1], title="uB"*"\ntime=$round_time",titlesize=fontsize_title)
    ax_uB.aspect = AxisAspect(1)
    xlims!(ax_uB,0,Lx)
    ylims!(ax_uB,0,Ly)

    if isnothing(range_uB)
        @views hm_uB = heatmap!(ax_uB,xx,yy,uB,colormap=colormap)
    else
        @views hm_uB = heatmap!(ax_uB,xx,yy,uB,colorrange=range_uB,colormap=colormap)
    end
    Colorbar(fig[2,2],hm_uB,height=250)

    ax_vB = Axis(fig[2, 3], title="vB"*"\ntime=$round_time",titlesize=fontsize_title)
    ax_vB.aspect = AxisAspect(1)
    xlims!(ax_vB,0,Lx)
    ylims!(ax_vB,0,Ly)
    if isnothing(range_vB)
        @views hm_vB = heatmap!(ax_vB,xx,yy,vB,colormap=colormap)
    else
        @views hm_vB = heatmap!(ax_vB,xx,yy,vB,colorrange=range_vB,colormap=colormap)
    end
    Colorbar(fig[2,4],hm_vB,height=250)

    ax_c = Axis(fig[3, 1], title="cells\ntime=$round_time",titlesize=fontsize_title)
    ax_c.aspect = AxisAspect(1)
    xlims!(ax_c,0,Lx)
    ylims!(ax_c,0,Ly)
    if isnothing(range_c)
        @views hm_c = heatmap!(ax_c,xx,yy,c,colormap=colormap)
    else
        @views hm_c = heatmap!(ax_c,xx,yy,c,colorrange=range_c,colormap=colormap)
    end
    Colorbar(fig[3,2],hm_c,height=250)

    rowsize!(fig.layout,1,Auto(0.5))

    if save_video
        stream = VideoStream(fig,framerate=40)
        recordframe!(stream)
        return fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS, stream
    else
        return fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS, nothing
    end
end


function init_plot_1D(
    c,uB,vB,uS,vS;
    dx,dy=nothing,
    save_video=false,
    time=0.0,
    range_uS=nothing,
    range_vS=nothing,
    range_uB=nothing,
    range_vB=nothing,
    range_c=nothing,
    colormap=nothing,
    fontsize_label=16,
    fontsize_title=24
)
    round_time = nice_float2string(time,2)
    Nx = size(uB)[1]
    Lx = dx * Nx

    xx = dx .* collect(0:(Nx-1)) .+ dx/2

    fig = Figure(resolution=(1100,800),fontsize=fontsize_label)

    ax_uS = Axis(fig[1, 1], title="uS\ntime=$round_time",titlesize=fontsize_title)
    ax_vS = Axis(fig[1, 2], title="vS\ntime=$round_time",titlesize=fontsize_title)
    ax_c = Axis(fig[2, 3], title="c\ntime=$round_time",titlesize=fontsize_title)
    ax_uB = Axis(fig[2, 1], title="uB\ntime=$round_time",titlesize=fontsize_title)
    ax_vB = Axis(fig[2, 2], title="vB\ntime=$round_time",titlesize=fontsize_title)

    xlims!(ax_uS,0,Lx)
    hm_uS = lines!(ax_uS,xx,uS,color=:black)

    xlims!(ax_uB,0,Lx)
    hm_uB = lines!(ax_uB,xx,uB,color=:black)

    xlims!(ax_vS,0,Lx)
    hm_vS = lines!(ax_vS,xx,vS,color=:black)

    xlims!(ax_vB,0,Lx)
    hm_vB = lines!(ax_vB,xx,vB,color=:black)

    xlims!(ax_c,0,Lx)
    hm_c = lines!(ax_c,xx,c,color=:black)

    if isnothing(range_uS)
        range_uS = nice_range(uS)
    end
    ylims!(ax_uS,range_uS[1],range_uS[2])
    
    if isnothing(range_vS)
        range_vS = nice_range(vS)
    end
    ylims!(ax_vS,range_vS[1],range_vS[2])

    if isnothing(range_uB)
        range_uB = nice_range(uB)
    end
    ylims!(ax_uB,range_uB[1],range_uB[2])        

    if isnothing(range_vB)
        range_vB = nice_range(vB)
    end
    ylims!(ax_vB,range_vB[1],range_vB[2])        

    if isnothing(range_c)
        range_c = nice_range(c)
    end
    ylims!(ax_c,range_c[1],range_c[2])


    if save_video
        stream = VideoStream(fig,framerate=40)
        recordframe!(stream)
        return fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS, stream
    else
        return fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS, nothing
    end
end


function update_plot_2D!(
    c,uB,vB,uS,vS;
    time,
    fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS,
    range_uS=nothing,
    range_vS=nothing,
    stream=nothing,
)    

    round_time = nice_float2string(time,2)

    @views hm_uS.input_args[2][] = uS
    @views hm_vS.input_args[2][] = vS

    if isnothing(range_uS)
        range_uS = nice_range(uS)
    end
    ylims!(ax_uS,range_uS[1],range_uS[2])

    
    if isnothing(range_vS)
        range_vS = nice_range(vS)
    end
    ylims!(ax_vS,range_vS[1],range_vS[2])

    ax_uS.title = "uS\ntime=$round_time"
    ax_vS.title = "vS\ntime=$round_time"

    hm_uB[3] = uB
    ax_uB.title = "uB\ntime=$round_time"

    hm_vB[3] = vB
    ax_vB.title = "vB\ntime=$round_time"

    hm_c[3] = c
    ax_c.title = "cells\ntime=$round_time"

    if !isnothing(stream)
        recordframe!(stream)
    end
end

function update_plot_1D!(
    c,uB,vB,uS,vS;
    time,
    fig, ax_uS, ax_vS,ax_uB, ax_vB, ax_c, hm_uB, hm_vB, hm_uS, hm_vS, hm_c,
    stream=nothing,
    range_uS=nothing,
    range_vS=nothing,
    range_uB=nothing,
    range_vB=nothing,
    range_c=nothing
)    
    round_time = nice_float2string(time,2)

    hm_uS.input_args[2][] = uS
    hm_vS.input_args[2][] = vS
    hm_uB.input_args[2][] = uB
    hm_vB.input_args[2][] = vB
    hm_c.input_args[2][] = c
    # ylims!(ax_u,0.9*min(minimum(uS),minimum(uB)),1.1*max(maximum(uS),maximum(uB)))
    # ylims!(ax_v,0.9*min(minimum(vS),minimum(vB)),1.1*max(maximum(vS),maximum(vB)))
    # ylims!(ax_c,0.9*min(minimum(c),minimum(c)),1.1*max(maximum(c),maximum(c)))
    if isnothing(range_uS)
        range_uS = nice_range(uS)
    end
    ylims!(ax_uS,range_uS[1],range_uS[2])

    
    if isnothing(range_vS)
        range_vS = nice_range(vS)
    end
    ylims!(ax_vS,range_vS[1],range_vS[2])

    if isnothing(range_uB)
        range_uB = nice_range(uB)
    end
    ylims!(ax_uB,range_uB[1],range_uB[2])        

    if isnothing(range_vB)
        range_vB = nice_range(vB)
    end
    ylims!(ax_vB,range_vB[1],range_vB[2])        

    if isnothing(range_c)
        range_c = nice_range(c)
    end
    ylims!(ax_c,range_c[1],range_c[2])

    ax_uS.title = "uS\ntime=$round_time"
    ax_vS.title = "vS\ntime=$round_time"
    ax_uB.title = "uB\ntime=$round_time"
    ax_vB.title = "vB\ntime=$round_time"
    ax_c.title = "c\ntime=$round_time"

    if !isnothing(stream)
        recordframe!(stream)
    end
end


function nice_float2string(x,K::Int)
    integer_part = trunc(Int,x)
    y = x - integer_part
    y10K = trunc(Int,y*10^K)
    if y10K >= 10
        decimal_part = rpad(y10K,K,"0")
    else
        decimal_part = lpad(y10K,K,"0")
    end
    return "$(integer_part).$decimal_part"
end


function plot_save_sol(
    sol,init,c_index,uB_index,vB_index,uS_index,vS_index,dx,dy;
    dt=0.1,
    range_uS=nothing,
    range_vS=nothing,
    range_uB=nothing,
    range_vB=nothing,
    range_c=nothing,
    colormap=:greys,
    fontsize_label=16,
    fontsize_title=24,
    dir=pwd(),video_name="funny_video")

    dim = length(size(init)) - 1
    if dim == 2
        init_plot = init_plot_2D
        update_plot! = update_plot_2D!
    else
        init_plot = init_plot_1D
        update_plot! = update_plot_1D!
    end
    
    println("Init plot...")
    
    uS = @view init[uS_index...]
    vS = @view init[vS_index...]
    uB = @view init[uB_index...]
    vB = @view init[vB_index...]
    c = @view init[c_index...]

    fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS, stream = init_plot(
        c,uB,vB,uS,vS;
        dx,dy,
        save_video=true,
        time=0.0,
        range_uB=range_uB,
        range_vB=range_vB,
        range_c=range_c,
        range_uS=range_uS,
        range_vS=range_vS,
        colormap=colormap,
        fontsize_label=fontsize_label,
        fontsize_title=fontsize_title,
    )

    T = sol.t[end]
    K = floor(Int,T/dt)

    println("Plotting...")

    @showprogress for k in 1:K
        U = sol(k*dt)
        time = k*dt
        uS = @view U[uS_index...]
        vS = @view U[vS_index...]
        uB = @view U[uB_index...]
        vB = @view U[vB_index...]
        c = @view U[c_index...]
        update_plot!(
        c,uB,vB,uS,vS;
        time,
        fig, ax_c, ax_uB, ax_vB, ax_uS, ax_vS, hm_c, hm_uB, hm_vB, hm_uS, hm_vS,
        range_uS=range_uS,
        range_vS=range_vS,
        stream=stream)
    end

    println("Save plot...")
    
    save(dir*"/"*video_name*".mp4", stream)
end


function init_directory(;simu_name="simu")
    if ispath(simu_name)
        k=1
        while ispath(simu_name*"_$k")
            k+=1
        end
        mkdir(simu_name*"_$k")
        return simu_name*"_$k"
    else
        mkdir(simu_name)
        return simu_name
    end
end

function nice_range(u)
    if abs(minimum(u) - maximum(u))<1e-6
        return (minimum(u)-1e-3,maximum(u)+1e-3)
    else
        return (minimum(u) - 0.05*(maximum(u)-minimum(u)), maximum(u) + 0.05*(maximum(u)-minimum(u)))
    end
end