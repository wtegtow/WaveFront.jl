using WaveFront
using Test
using GLMakie

function tt_ana_layered_3d(z_coords, velocity_layers, source_z)
    nz = length(z_coords)
    tt = zeros(nz)

    iz_source = findmin(abs.(z_coords .- source_z))[2]

    for iz in 1:nz
        iz_min = min(iz, iz_source)
        iz_max = max(iz, iz_source)

        tsum = 0.0
        for k in iz_min:iz_max-1
            thickness = abs(z_coords[k+1] - z_coords[k])
            v_layer = velocity_layers[k]
            tsum += thickness / v_layer
        end
        last_thickness = abs(z_coords[iz] - z_coords[iz_max])
        tsum += last_thickness / velocity_layers[iz]
        tt[iz] = tsum
    end

    return tt
end

function test_het3d(;plotit=false, verbose=false)

    MAX_ERROR_PERCENT = 5

    h = 10
    x_coords = 0:h:1000
    y_coords = 0:h:1000
    z_coords = 0:h:1000

    velocity_layers = rand(1000:4000, length(z_coords))
    velocity = zeros(length(x_coords), length(y_coords), length(z_coords))
    for iz in 1:length(z_coords)
        velocity[:, :, iz] .= velocity_layers[iz]
    end

    source_x = x_coords[Int(floor(length(x_coords)/2))]
    source_y = y_coords[Int(floor(length(y_coords)/2))]
    source_z = 0
    source = (source_x, source_y, source_z)

    grid = Grid3D(x_coords, y_coords, z_coords, velocity)

    tt_num = fast_sweeping(grid, [source], ϵ=0.1, verbose=false)
    tt_ana = tt_ana_layered_3d(z_coords, velocity_layers, source_z)

    ix = Int(floor(length(x_coords)/2))
    iy = Int(floor(length(y_coords)/2))

    tt_profile_num = tt_num[ix, iy, :]

    max_error_abs = maximum(abs.(tt_profile_num .- tt_ana))
    max_tt_ana = maximum(tt_ana)

    max_error_percent = 100 * max_error_abs / max_tt_ana

    if verbose println("Max Error %: ", max_error_percent) end

    @test max_error_percent < MAX_ERROR_PERCENT

    if plotit
        Makie.inline!(true)

        fig = Figure(size=(600,600))

        ax1 = Axis(fig[1,1], xlabel="x", ylabel="z", title="Velocity Slice", yreversed=true)
        heatmap!(ax1, x_coords, z_coords, velocity[:, iy, :])

        ax2 = Axis(fig[2,1], xlabel="x", ylabel="z", title="Travel Time Slice", yreversed=true)
        heatmap!(ax2, x_coords, z_coords, tt_num[:, iy, :])

        ax3 = Axis(fig[3,1], xlabel="z", ylabel="Travel Time", title="1D Profile")
        num = lines!(ax3, z_coords, tt_profile_num)
        ana = lines!(ax3, z_coords, tt_ana, linestyle=:dash)

        Legend(fig[3,2], [ana, num], ["Ana", "Num"])
        display(fig)
    end

end

test_het3d(;plotit=false, verbose=false)