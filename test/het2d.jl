using WaveFront
using Test
using CairoMakie

function tt_ana_layered(y_coords, velocity_layers, source_y)

    ny = length(y_coords)
    tt = zeros(ny)
    
    iy_source = findmin(abs.(y_coords .- source_y))[2]
    
    for iy in 1:ny
        iy_min = min(iy, iy_source)
        iy_max = max(iy, iy_source)
        
        tsum = 0.0
        for k in iy_min:iy_max-1
            thickness = abs(y_coords[k+1] - y_coords[k])
            v_layer = velocity_layers[k]
            tsum += thickness / v_layer
        end
        last_thickness = abs(y_coords[iy] - y_coords[iy_max])
        tsum += last_thickness / velocity_layers[iy]
        tt[iy] = tsum
    end
    
    return tt
end

function test_het2d(plotit=false)

    MAX_ERROR = 0.1

    h = 10
    x_coords = 0:h:1000
    y_coords = 0:h:1000

    velocity_layers = zeros(length(y_coords))
    for i in eachindex(y_coords)
        velocity_layers[i] = rand(1000:4000)
    end

    velocity = repeat(velocity_layers', length(x_coords), 1)

    source_x = x_coords[Int(floor(length(x_coords)/2))]
    source_y = 0
    source = (source_x, source_y)

    grid = Grid2D(x_coords, y_coords, velocity)

    tt_num = fast_sweeping(grid, [source], verbose=false)
    tt_ana = tt_ana_layered(y_coords, velocity_layers, source_y)
    max_error = maximum(abs.(tt_num[Int(floor(length(x_coords)/2)), :] .- tt_ana))

    @test max_error < MAX_ERROR
    
    if plotit 

        fig = Figure(size=(500,500))

        ax1 = Axis(fig[1,1], xlabel="x", ylabel="y", title="Velmod", yreversed=true)
        heatmap!(ax1, x_coords, y_coords, velocity)

        ax2 = Axis(fig[2,1], xlabel="x", ylabel="y", title="Travel Time", yreversed=true)
        contourf!(ax2, x_coords, y_coords, tt_num)

        ax3 = Axis(fig[3,1], xlabel="y", ylabel="Travel Time", title="Profile")
        num = lines!(ax3, y_coords, tt_num[Int(floor(length(x_coords)/2)), :], label="Numerical")
        ana = lines!(ax3, y_coords, tt_ana, label="Analytical", linestyle=:dash)

        Legend(fig[3, 2],[ana, num],["ana", "num"])
        display(fig)
    end

end

test_het2d()