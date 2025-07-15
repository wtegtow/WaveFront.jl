module WaveFront

export fast_sweeping, Grid2D

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
    T = fill(Inf, nx, ny)
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

    for iter in 1:max_iter

        T_old = copy(T)

        Threads.@threads for (x_order, y_order) in sweeping_orders

            for i in x_order, j in y_order

                if (i,j) in sources
                    continue
                end

                T_xm = i == 1  ? Inf : T[i-1, j]
                T_xp = i == nx ? Inf : T[i+1, j]
                T_ym = j == 1  ? Inf : T[i, j-1]
                T_yp = j == ny ? Inf : T[i, j+1]

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

            end 
        end 

        if sum(abs.(T - T_old))  < ϵ
            if verbose println("Solution converged after $(iter) iterations.") end 
            break
        end
    end
    return T 
end

# ====================================================
# 3D
# ====================================================

struct Grid3D{V<:AbstractVector, M<:AbstractMatrix}
    x_coords::V
    y_coords::V
    z_coords::V
    velocity::M
end



end