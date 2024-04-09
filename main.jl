using Distributions
using ForwardDiff
using JLD2
using LinearAlgebra
using NLsolve
using ProgressMeter
using QuadGK
using Random
using SpecialFunctions

## PREPARE DIRECTORIES
dirs = ["Data", "Data/Frozen", "Data/Deflected", "Data/Thermal", "Data/Joint", "Data/Phase"]
[isdir(d) ? nothing : mkdir(d) for d in dirs]

@inline function derivatives(Φ, s, μ, ρ, dρ, σ, dσ; frozen = false)
    dist = sqrt.((ρ .- s) .^ 2 + σ .^ 2)
    dΦ = [ForwardDiff.derivative(r -> Φ(r), d) for d in dist]

    ddσ = -(4 * π^2 / μ) .* (dΦ .* σ ./ dist)
    ddρ = -(4 * π^2) .* (ρ + dΦ .* (ρ .- s) ./ dist)

    if frozen
        return (0 .* dρ, 0 .* ddρ, dσ, ddσ)
    else
        (dρ, ddρ, dσ, ddσ)
    end
end

function RKstep(Φ, s, μ, current_state, δτ; frozen = false)
    k1 = derivatives(Φ, s, μ, current_state...; frozen)
    k2 = derivatives(Φ, s, μ, (current_state .+ k1 .* (δτ / 4))...; frozen)
    k3 = derivatives(Φ, s, μ, (current_state .+ (k1 .+ k2) .* (δτ / 8))...; frozen)
    k4 = derivatives(Φ, s, μ, (current_state .+ k3 .* δτ .- k2 .* (δτ / 2))...; frozen)
    k5 = derivatives(
        Φ,
        s,
        μ,
        (current_state .+ k1 .* (δτ * 3 / 16) .+ k4 .* (δτ * 9 / 16))...;
        frozen,
    )
    k6 = derivatives(
        Φ,
        s,
        μ,
        (
            current_state .- k1 .* (3 / 7 * δτ) .+ k2 .* (2 / 7 * δτ) .+
            k3 .* (12 / 7 * δτ) .- k4 .* (12 / 7 * δτ) .+ k5 .* (8 / 7 * δτ)
        )...;
        frozen,
    )
    res =
        current_state .+
        (δτ / 90) .* (7 .* k1 .+ 32 .* k3 .+ 12 .* k4 .+ 32 .* k5 .+ 7 .* k6)
    return res
end

function solverRK(Φ, s, μ, τmax, δτ, curr_state; show_prog = true, frozen = false)
    nτs = (floor(τmax / δτ) |> Int) + 1   # Number of time steps
    res = zeros(4, nτs)
    res[:, 1] = [x for x in curr_state]

    prog = Progress(nτs)

    for ii = 2:nτs
        curr_state = RKstep(Φ, s, μ, curr_state, δτ; frozen = frozen)
        res[:, ii] = [x for x in curr_state]
        if show_prog
            next!(prog)
        end
    end
    return res
end

# function passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
#     σ = curr_state[3]
#     dσ = curr_state[4]

#     while (dσ > 0) & (σ < 0)
#         curr_state = RKstep(Φ, s, μ, curr_state, δτ; frozen = frozen)
#         σ = curr_state[3]
#         dσ = curr_state[4]
#     end
#     return (σ > 0)
# end


function passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
    σ = curr_state[3]
    dσ = curr_state[4]

    passed = 0

    while length(σ) > 0

        curr_state = RKstep(Φ, s, μ, curr_state, δτ; frozen = frozen)
        ρ = curr_state[1]
        dρ = curr_state[2]
        σ = curr_state[3]
        dσ = curr_state[4]

        # Find how many particles have passed the oscillator
        passed = passed + sum(σ .> 0)

        # Keep only the particles that have not bounced back or passed
        idx = findall(x -> x < 0, σ .* dσ)

        σ = σ[idx]
        dσ = dσ[idx]
        ρ = ρ[idx]
        dρ = dρ[idx]
        curr_state = (ρ, dρ, σ, dσ)
    end
    return (passed)
end
