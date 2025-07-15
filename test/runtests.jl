using WaveFront
using Test

@testset "WaveFront.jl" begin
    include(joinpath(@__DIR__, "hom2d.jl"))
    include(joinpath(@__DIR__, "het2d.jl"))

    include(joinpath(@__DIR__, "hom3d.jl"))
    include(joinpath(@__DIR__, "het3d.jl"))
end
