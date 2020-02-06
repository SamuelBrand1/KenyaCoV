using Test
using KenyaCoV
include("test_functions.jl");
@testset "Simulation tests" begin
    @test test_no_infecteds()
end
