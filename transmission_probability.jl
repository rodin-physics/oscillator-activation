include("main.jl")

# Simulation parameters
seed_number = 150

# High temperature
ωT_min = 4
ωT_max = 150
nPts = 200

ωTs_high = range(ωT_min, ωT_max, length = nPts)

# Low temperature
ωT_min = 1
ωT_max = 4
nPts = 50

ωTs_low = range(ωT_min, ωT_max, length = nPts)

ωTs = vcat(ωTs_low, ωTs_high)
nPts = length(ωTs)

nTrials = 10000  # Number of runs at each pamaters combination
δτ = 1e-3 / 2

σ = -20
impact_parameters = [0, 2]
μs = [1 / 5, 1, 5]
params = [(s, μ) for s in impact_parameters, μ in μs] |> vec

# Interaction function
Φ0 = 10
λ = 2
@inline function Φ(r)
    res = Φ0 * exp(-r^2 / 2 / λ^2)
    return res
end

## FROZEN MASS
println("Frozen mass")

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Frozen/Frozen_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")
        transmission = zeros(nPts)
        prog = Progress(nPts)
        Threads.@threads for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            rand_energy = -ωT * log.(1 .- rand(Float64, nTrials))
            init_speed = sqrt.(rand_energy * 8 * π^2 / μ)
            curr_state = (zeros(nTrials), zeros(nTrials), σ .* ones(nTrials), init_speed)
            transmission[ii] = passCheck(Φ, s, μ, δτ, curr_state; frozen = true) / nTrials

            next!(prog)
            GC.safepoint()
        end
        save_object(
            "Data/Frozen/Frozen_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ωTs, transmission, s, μ),
        )
    end
end

## DEFLECTED MASS
println("Deflected mass")

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Deflected/Deflected_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")
        transmission = zeros(nPts)
        prog = Progress(nPts)
        Threads.@threads for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            rand_energy = -ωT * log.(1 .- rand(Float64, nTrials))
            init_speed = sqrt.(rand_energy * 8 * π^2 / μ)
            curr_state = (zeros(nTrials), zeros(nTrials), σ .* ones(nTrials), init_speed)
            transmission[ii] = passCheck(Φ, s, μ, δτ, curr_state; frozen = false) / nTrials

            next!(prog)
            GC.safepoint()
        end
        save_object(
            "Data/Deflected/Deflected_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ωTs, transmission, s, μ),
        )
    end
end

# THERMAL MASS
println("Thermal mass")

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Thermal/Thermal_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")
        transmission = zeros(nPts)
        prog = Progress(nPts)
        Threads.@threads for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            rand_energy = -ωT * log.(1 .- rand(Float64, nTrials))
            init_speed = sqrt.(rand_energy * 8 * π^2 / μ)

            η = 1e-12
            n = rand(Geometric(1 - exp(-1 / ωT) - η), nTrials)
            ζ = sqrt.(n .+ 1 / 2) .* √(2)
            ϕ = 2 * π * rand(Float64, nTrials)

            ρ0 = ζ .* sin.(ϕ)
            dρ0 = 2 * π * ζ .* cos.(ϕ)

            curr_state = (ρ0, dρ0, σ .* ones(nTrials), init_speed)

            transmission[ii] = passCheck(Φ, s, μ, δτ, curr_state; frozen = false) / nTrials

            next!(prog)
            GC.safepoint()
        end
        save_object(
            "Data/Thermal/Thermal_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ωTs, transmission, s, μ),
        )
    end
end

# THERMAL MASS JOINT
println("Thermal mass joint")

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Joint/Joint_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")
        transmission = zeros(nPts)
        prog = Progress(nPts)
        Threads.@threads for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            rand_energy = -ωT * log.(1 .- rand(Float64, nTrials))
            nMax = [floor(Int, r) for r in rand_energy]
            n = [rand(0:n) for n in nMax]

            particle_energy = rand_energy .- n
            init_speed = sqrt.(particle_energy * 8 * π^2 / μ)

            ζ = sqrt.(n .+ 1 / 2) .* √(2)
            ϕ = 2 * π * rand(Float64, nTrials)

            ρ0 = ζ .* sin.(ϕ)
            dρ0 = 2 * π * ζ .* cos.(ϕ)

            curr_state = (ρ0, dρ0, σ .* ones(nTrials), init_speed)

            transmission[ii] = passCheck(Φ, s, μ, δτ, curr_state; frozen = false) / nTrials

            next!(prog)
            GC.safepoint()
        end
        save_object(
            "Data/Joint/Joint_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ωTs, transmission, s, μ),
        )
    end
end
