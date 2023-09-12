include("main.jl")

# Simulation parameters
seed_number = 150

δτ = 1e-3 / 2

σ = -20
impact_parameters = [0, 2]
# impact_parameters = [0, 1, 2, 4]
μs = [1]
# impact_parameters = [0]
# μs = [1]
params = [(s, μ) for s in impact_parameters, μ in μs] |> vec

nPts = 200
en_min_particle = 1 / 10
en_max_particle = 12

en_min_osc = 0 / 10
en_max_osc = 20

ens_particle = range(en_min_particle, en_max_particle, length = nPts)
ens_osc = range(en_min_osc, en_max_osc, length = nPts)
# ϕs = range(0, 2 * π, length = nPts)
ϕs = (0:8) .* π / 4

# Interaction function
Φ0 = 10
λ = 2
@inline function Φ(r)
    res = Φ0 * exp(-r^2 / 2 / λ^2)
    return res
end

idx =
    [
        (ii, jj, kk) for ii in eachindex(ens_particle), jj in eachindex(ens_osc),
        kk in eachindex(ϕs)
    ] |> vec

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Phase/Phase_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")

        transmission = zeros(Int, (length(ens_particle), length(ens_osc), length(ϕs)))

        prog = Progress(length(idx))

        Threads.@threads for run in idx
            ii = run[1]
            jj = run[2]
            kk = run[3]
            particle_energy = ens_particle[ii]
            osc_energy = ens_osc[jj]
            ϕ = ϕs[kk]

            init_speed = √(particle_energy * 8 * π^2 / μ)

            ζ = √(osc_energy) * √(2)

            ρ0 = ζ * sin(ϕ)
            dρ0 = 2 * π * cos(ϕ)
            curr_state = (ρ0, dρ0, σ, init_speed)
            τmax = 2 * σ / init_speed |> abs
            res = passCheck(Φ, s, μ, δτ, curr_state; frozen = false)
            transmission[ii, jj, kk] = res
            next!(prog)
            GC.safepoint()
        end
        save_object(
            "Data/Phase/Phase_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2",
            (ens_particle, ens_osc, ϕs, transmission, s, μ),
        )
    end
end
