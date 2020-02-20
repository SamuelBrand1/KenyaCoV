using Test
import KenyaCoV
include("test_functions.jl");
@testset "Simulation tests" begin
    @test test_no_infecteds()
    @test stays_nonnegative()
end
