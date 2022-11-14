using Test

push!(LOAD_PATH, "../src")
using PorousConvection

@test PorousConvection.TEST_VAR == true