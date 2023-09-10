include("main.jl")

# Simulation parameters
seed_number = 150

δτ = 1e-3 / 2

σ = -20
impact_parameters = [0, 1, 2, 4]
μs = [1 / 5, 1, 5]
# impact_parameters = [0]
# μs = [1]
params = [(s, μ) for s in impact_parameters, μ in μs] |> vec

nPts = 100
en_min = 1
en_max = 12
ens_particle = range(en_min, en_max, length = nPts)
ens_osc = range(en_min, en_max, length = nPts)
ϕs = range(0, 2 * π, length = nPts)


# Interaction function
Φ0 = 10
λ = 2
@inline function Φ(r)
    res = Φ0 * exp(-r^2 / 2 / λ^2)
    return res
end

idx = [(ii, jj, kk) for ii = 1:nPts, jj = 1:nPts, kk = 1:nPts] |> vec

for p in params
    s = p[1]
    μ = p[2]
    Random.seed!(seed_number)
    if !isfile("Data/Phase/Phase_s$(s)_Φ$(Φ0)_λ$(λ)_μ$(μ).jld2")

        transmission = zeros(Int, (nPts, nPts, nPts))

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


# tt = load_object("Data/Phase/Phase_s0_Φ10_λ2_μ1.jld2")
# r = tt[end-2]
# heatmap(r[:,:,90])
# [:,:,50]
