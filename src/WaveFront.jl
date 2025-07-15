module WaveFront
export fast_sweeping, Grid2D, Grid3D

# ====================================================
# 2D
# ====================================================

struct Grid2D{V<:AbstractVector, M<:AbstractMatrix}
    x_coords::V
    y_coords::V
    velocity::M
end

function fast_sweeping(grid::Grid2D, sources_phys; ϵ = 1e-6, max_iter=50, verbose=false)

    x_coords = grid.x_coords
    Δx = x_coords[2] - x_coords[1]
    nx = length(x_coords)

    y_coords = grid.y_coords
    Δy = y_coords[2] - y_coords[1]
    ny = length(y_coords)

    C = grid.velocity
    BIG_VAL = 1e10
    T = fill(BIG_VAL, nx, ny)
    @assert size(T) == size(C) 

    sources = [(argmin(abs.(s[1] .- x_coords)),
                argmin(abs.(s[2] .- y_coords))) for s in sources_phys]

    for s in sources
        T[s...] = 0.0
    end

    sweeping_orders = (
        (1:1:nx , 1:1:ny), 
        (nx:-1:1, 1:1:ny), 
        (1:1:nx, ny:-1:1),
        (nx:-1:1, ny:-1:1)
    )

    converged = false
    for iter in 1:max_iter

        T_old = copy(T)

        for (x_order, y_order) in sweeping_orders

            for i in x_order, j in y_order

                if (i,j) in sources
                    continue
                end

                T_xm = i == 1  ? BIG_VAL : T[i-1, j]
                T_xp = i == nx ? BIG_VAL : T[i+1, j]
                T_ym = j == 1  ? BIG_VAL : T[i, j-1]
                T_yp = j == ny ? BIG_VAL : T[i, j+1]

                T_xmin = min(T_xm, T_xp)
                T_ymin = min(T_ym, T_yp)
                
                # assign v <= u 
                if T_xmin <= T_ymin 
                    v = T_xmin 
                    Δv = Δx
                    u = T_ymin 
                    Δu = Δy
                else 
                    v = T_ymin 
                    Δv = Δy
                    u = T_xmin 
                    Δu = Δx
                end 

                c_loc = C[i,j]

                T_star = v + (Δv/c_loc)
                
                if T_star <= u 
                    T_new = T_star 
                else 
                    # solve quadratic equation
                    a = Δv^2 + Δu^2
                    b = -2*(Δv^2*u + Δu^2*v)
                    c_quad = Δv^2 * u^2 + Δu^2 * v^2 - (Δu^2 * Δv^2 / c_loc^2)
                    T_p = (-b + sqrt(b^2 - 4*a*c_quad)) / (2a)
                    T_m = (-b - sqrt(b^2 - 4*a*c_quad)) / (2a)
                    T_new = max(T_p, T_m)
                end

                T[i,j] = T_new
                
            end # i,j end 

        end # sweep order end 

        if sum(abs.(T - T_old))  < ϵ
            if verbose println("Solution converged after $(iter) iterations.") end 
            converged = true 
            break
        end
    end
    if !converged println("Solution not converged.") end 
    return T 
end

# ====================================================
# 3D
# ====================================================

struct Grid3D{V<:AbstractVector, A<:AbstractArray{<:Real,3}}
    x_coords::V
    y_coords::V
    z_coords::V
    velocity::A
end

function fast_sweeping(grid::Grid3D, sources_phys; ϵ = 0.1, max_iter=50, verbose=false)

    x_coords = grid.x_coords
    Δx = x_coords[2] - x_coords[1]
    nx = length(x_coords)

    y_coords = grid.y_coords
    Δy = y_coords[2] - y_coords[1]
    ny = length(y_coords)

    z_coords = grid.z_coords
    Δz = z_coords[2] - z_coords[1]
    nz = length(z_coords)

    C = grid.velocity
    BIG_VAL = 1e10
    T = fill(BIG_VAL, nx, ny, nz)
    @assert size(T) == size(C)

    sources = [(argmin(abs.(s[1] .- x_coords)),
                argmin(abs.(s[2] .- y_coords)),
                argmin(abs.(s[3] .- z_coords))) for s in sources_phys]

    for s in sources
        T[s...] = 0.0
    end

    sweeping_orders = (
        (1:1:nx, 1:1:ny, 1:1:nz),    
        (1:1:nx, 1:1:ny, nz:-1:1),   
        (1:1:nx, ny:-1:1, 1:1:nz),   
        (1:1:nx, ny:-1:1, nz:-1:1),  
        (nx:-1:1, 1:1:ny, 1:1:nz),   
        (nx:-1:1, 1:1:ny, nz:-1:1),  
        (nx:-1:1, ny:-1:1, 1:1:nz),  
        (nx:-1:1, ny:-1:1, nz:-1:1) 
    )

    converged = false
    for iter in 1:max_iter

        T_old = copy(T)

        for (x_order, y_order, z_order) in sweeping_orders

            for i in x_order, j in y_order, k in z_order 

                if (i,j,k) in sources 
                    continue
                end

                T_xm = i == 1   ? BIG_VAL : T[i-1, j, k]
                T_xp = i == nx  ? BIG_VAL : T[i+1, j, k]
                T_ym = j == 1   ? BIG_VAL : T[i, j-1, k]
                T_yp = j == ny  ? BIG_VAL : T[i, j+1, k]
                T_zm = k == 1   ? BIG_VAL : T[i, j, k-1]
                T_zp = k == nz  ? BIG_VAL : T[i, j, k+1]
                
                Tx_min = min(T_xm, T_xp)
                Ty_min = min(T_ym, T_yp)
                Tz_min = min(T_zm, T_zp)

                # assign w <= v <= u
                vals = [Tx_min, Ty_min, Tz_min]
                increments = [Δx, Δy, Δz]   
                perm = sortperm(vals)

                w, v, u = vals[perm]
                Δw, Δv, Δu = increments[perm]
                c_loc = C[i, j, k]
                
                # find solution
                sol_found = false 
                T_new = nothing

                # solve equation with 1 unknown 
                if !sol_found
                    T_star = w + Δw / c_loc
                    if T_star <= v
                        T_new = T_star
                        sol_found = true 
                    end 
                end 

                # solve equation with 2 unknowns
                if !sol_found
                    a = Δv^2 + Δw^2 
                    b = -2*(Δv^2*w + Δw^2*v)
                    c_quad = Δv^2 * w^2 + Δw^2 * v^2 - (Δv^2 * Δw^2 / c_loc)

                    T_p = (-b + sqrt(b^2 - 4*a*c_quad)) / (2a)
                    T_m = (-b - sqrt(b^2 - 4*a*c_quad)) / (2a)
                    T_star = max(T_p, T_m)

                    if T_star <= u
                        T_new = T_star
                        sol_found = true 
                    end
                end

                # solve equation with 3 unknowns
                if !sol_found

                    a = 1.0 / Δw^2 + 1.0 / Δv^2 + 1.0 / Δu^2
                    b = -2.0 * (w / Δw^2 + v / Δv^2 + u / Δu^2)
                    c_quad = (w^2 / Δw^2 + v^2 / Δv^2 + u^2 / Δu^2) - 1.0 / c_loc^2

                    D = b^2 - 4*a*c_quad

                    if D < 0
                        T_new = max(w, v, u) 
                    else
                        T_new = (-b + sqrt(D)) / (2a)
                    end
                end

                T[i,j,k] = T_new 

            end # i,j,k end

        end # sweep end 

        if maximum(abs.(T - T_old)) < ϵ
            if verbose println("Solution converged after $(iter) iterations.") end 
            converged = true 
            break
        end
    end # iter end 
    if !converged println("Solution not converged.") end 
    return T 
end # fun end 

end # module end 