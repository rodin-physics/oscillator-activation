include("main.jl")

# Simulation parameters
seed_number = 150
ωT_min = 1
ωT_max = 150
nPts = 100
ωTs = range(ωT_min, ωT_max, length = nPts)

nTrials = 1000  # Number of runs at each pamaters combination
δτ = 1e-3 / 2

σ = -20
impact_parameters = [0, 1, 2, 4]
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
        for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            passes = zeros(Int, nTrials)
            prog = Progress(nTrials)
            Threads.@threads for run = 1:nTrials

                rand_energy = -ωT * log(1 - rand())
                init_speed = √(rand_energy * 8 * π^2 / μ)
                curr_state = (0, 0, σ, init_speed)
                res = passCheck(Φ, s, μ, δτ, curr_state; frozen = true)
                passes[run] = res
                next!(prog)
                GC.safepoint()
            end
            transmission[ii] = sum(passes) ./ nTrials
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
        for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            passes = zeros(Int, nTrials)
            prog = Progress(nTrials)
            Threads.@threads for run = 1:nTrials

                rand_energy = -ωT * log(1 - rand())
                init_speed = √(rand_energy * 8 * π^2 / μ)
                curr_state = (0, 0, σ, init_speed)
                res = passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
                passes[run] = res
                next!(prog)
                GC.safepoint()
            end
            transmission[ii] = sum(passes) ./ nTrials
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
        for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            passes = zeros(Int, nTrials)
            prog = Progress(nTrials)
            Threads.@threads for run = 1:nTrials
                rand_energy = -ωT * log(1 - rand())
                init_speed = √(rand_energy * 8 * π^2 / μ)

                η = 1e-12
                n = rand(Geometric(1 - exp(-1 / ωT) - η))
                ζ = √(n + 1 / 2) * √(2)
                ϕ = 2 * π * rand()

                ρ0 = ζ * sin(ϕ)
                dρ0 = 2 * π * cos(ϕ)

                curr_state = (ρ0, dρ0, σ, init_speed)
                res = passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
                passes[run] = res
                next!(prog)
                GC.safepoint()
            end
            transmission[ii] = sum(passes) ./ nTrials
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
        for ii in eachindex(ωTs)
            ωT = ωTs[ii]
            passes = zeros(Int, nTrials)
            prog = Progress(nTrials)
            Threads.@threads for run = 1:nTrials
                rand_energy = -ωT * log(1 - rand())
                nMax = floor(Int, rand_energy)
                n = rand(0:nMax)

                particle_energy = rand_energy - n
                init_speed = √(particle_energy * 8 * π^2 / μ)

                ζ = √(n + 1 / 2) * √(2)
                ϕ = 2 * π * rand()

                ρ0 = ζ * sin(ϕ)
                dρ0 = 2 * π * cos(ϕ)
                curr_state = (ρ0, dρ0, σ, init_speed)
                res = passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
                passes[run] = res
                next!(prog)
                GC.safepoint()
            end
            transmission[ii] = sum(passes) ./ nTrials
        end
        save_object(
            "Data/Joint/Joint_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ωTs, transmission, s, μ),
        )
    end
end
