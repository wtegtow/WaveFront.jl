using WaveFront
using Test
#using GLMakie 

function tt_ana_hom(x_coords, y_coords, velocity, source)
    nx, ny = length(x_coords), length(y_coords)
    x0, y0 = source
    tt = zeros(nx, ny)
    for ix in 1:nx, iy in 1:ny
        x, y = x_coords[ix], y_coords[iy]
        dist = sqrt((x - x0)^2 + (y - y0)^2)
        tt[ix, iy] = dist / velocity[ix, iy]
    end
    return tt
end

function test_hom2d(;plotit=false)

    MAX_ERROR = 0.1

    x_coords = 0:5:1000 
    y_coords = 0:5:1000
    source_coord = (x_coords[Int(floor(length(x_coords)/2))], 
                    y_coords[Int(floor(length(y_coords)/2))])
    velocity = ones(length(x_coords), length(y_coords)) .* 500;

    grid = Grid2D(x_coords, y_coords, velocity)

    tt_num = fast_sweeping(grid, [source_coord], verbose=false);
    tt_ana = tt_ana_hom(x_coords, y_coords, velocity, source_coord)
    misfit = tt_num .- tt_ana
    max_error = maximum(abs.(tt_num .- tt_ana))

    @test max_error < MAX_ERROR

    if plotit 

        #Makie.inline!(true)

        fig = Figure()
        ax1 = Axis(fig[1,1], yreversed=true, title="Eikonal", xlabel="x", ylabel="y")
        ax2 = Axis(fig[1,2], yreversed=true, title="Analytical", xlabel="x", ylabel="y")
        ax3 = Axis(fig[1,3], yreversed=true, title="Missfit", xlabel="x", ylabel="y")
        
        contourf!(ax1, grid.x_coords, grid.y_coords, tt_num, levels=10)
        contourf!(ax2, grid.x_coords, grid.y_coords, tt_ana, levels=10)
        contourf!(ax3, grid.x_coords, grid.y_coords, misfit, levels=10)

        #Colorbar(fig[2,1], limits=(minimum(tt_num), maximum(tt_num)), vertical=false, label="travel time")
        #Colorbar(fig[2,2], limits=(minimum(tt_ana), maximum(tt_ana)), vertical=false, label="travel time")
        #Colorbar(fig[2,3], limits=(minimum(misfit), maximum(misfit)), vertical=false, label="travel time")

        display(fig)

    end

end 

test_hom2d(plotit=false)