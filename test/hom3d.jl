using WaveFront
#using GLMakie


function tt_ana_hom_3d(x_coords, y_coords, z_coords, velocity, source)
    nx, ny, nz = length(x_coords), length(y_coords), length(z_coords)
    x0, y0, z0 = source
    tt = zeros(nx, ny, nz)
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        x, y, z = x_coords[ix], y_coords[iy], z_coords[iz]
        dist = sqrt((x - x0)^2 + (y - y0)^2 + (z - z0)^2)
        tt[ix, iy, iz] = dist / velocity[ix, iy, iz]
    end
    return tt
end

function test_hom2d(;plotit=false, verbose=false)

    MAX_ERROR_PERCENT = 5

    h = 20
    x_coords = 0:h:1000 
    y_coords = 0:h:1000
    z_coords = 0:h:1000
    source_coords = [(500, 500, 500)]
    vel0 = 1
    velocity = fill(vel0, length(x_coords), length(y_coords), length(z_coords)) 

    grid = Grid3D(x_coords, y_coords, z_coords, velocity)

    tt_num = fast_sweeping(grid, source_coords; ϵ=1e-4, verbose=false)
    tt_ana = tt_ana_hom_3d(x_coords, y_coords, z_coords, velocity, source_coords[1])

    max_error_abs = maximum(abs.(tt_num .- tt_ana))
    max_tt_ana = maximum(tt_ana)
    max_error_percent = 100 * max_error_abs / max_tt_ana
    if verbose println("Max Error % ", max_error_percent) end

    @test max_error_percent < MAX_ERROR_PERCENT

    if plotit  
        #Makie.inline!(true)
        fig = Figure(size=(800,800)) 

        midz = Int(floor(length(z_coords)/2))
        midy = Int(floor(length(y_coords)/2))
        midx = Int(floor(length(x_coords)/2))

        xy = [tt_num[:,:,midz], tt_ana[:,:,midz], abs.(tt_num[:,:,midz] .- tt_ana[:,:,midz])]
        xz = [tt_num[:,midy,:], tt_ana[:,midy,:], abs.(tt_num[:,midy,:] .- tt_ana[:,midy,:])]
        yz = [tt_num[midx,:,:], tt_ana[midx,:,:], abs.(tt_num[midx,:,:] .- tt_ana[midx,:,:])]

        titles = ["Num", "Ana", "Missfit"]
        cmap = :viridis
        for i in 1:3
            # XY
            ax_xy = Axis(fig[i,1], title="XY "*titles[i], yreversed=true)
            hm_xy = contourf!(ax_xy, x_coords, y_coords, xy[i], colormap=cmap)
            Colorbar(fig[i,2], hm_xy, width=15, label="s")

            # XZ
            ax_xz = Axis(fig[i,3], title="XZ "*titles[i], yreversed=true)
            hm_xz = contourf!(ax_xz, x_coords, z_coords, xz[i], colormap=cmap)
            Colorbar(fig[i,4], hm_xz, width=15, label="s")

            # YZ
            ax_yz = Axis(fig[i,5], title="YZ "*titles[i], yreversed=true)
            hm_yz = contourf!(ax_yz, y_coords, z_coords, yz[i], colormap=cmap)
            Colorbar(fig[i,6], hm_yz, width=15, label="s")
        end

        display(fig)

    end 
end 

test_hom2d(plotit=false, verbose=false)

