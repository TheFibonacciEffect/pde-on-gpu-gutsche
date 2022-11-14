using Test

push!(LOAD_PATH, "../src")
using PorousConvection2D

@test PorousConvection.TEST_VAR == true