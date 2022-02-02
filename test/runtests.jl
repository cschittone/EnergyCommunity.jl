using EnergyCommunity
using Test, FileIO, GLPK, MathOptInterface

const MOI = MathOptInterface

# EC groups
const EC_GROUPS = [EnergyCommunity.GroupCO(), EnergyCommunity.GroupNC()]


include("tests.jl")

@testset "EnergyCommunity tests" begin

    # Loop over group types
    for group in EC_GROUPS

        _base_test(GPLK.Optimizer, group)

    end

end