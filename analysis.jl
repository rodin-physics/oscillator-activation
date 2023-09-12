using CurveFit
using JLD2
## LOAD THE DATA
deflected = readdir("Data/Deflected/"; join = true)
frozen = readdir("Data/Frozen/"; join = true)
thermal = readdir("Data/Thermal/"; join = true)
joint = readdir("Data/Joint/"; join = true)

deflected = load_object.(deflected)
frozen = load_object.(frozen)
thermal = load_object.(thermal)
joint = load_object.(joint)

impact_parameters = unique([x[3] for x in deflected])
masses = unique([x[4] for x in deflected])

s0Frozen = filter(x -> x[3] == 0, frozen)
s2Frozen = filter(x -> x[3] == 2, frozen)

s0Deflected = filter(x -> x[3] == 0, deflected)
s2Deflected = filter(x -> x[3] == 2, deflected)

s0Thermal = filter(x -> x[3] == 0, thermal)
s2Thermal = filter(x -> x[3] == 2, thermal)

s0Joint = filter(x -> x[3] == 0, joint)
s2Joint = filter(x -> x[3] == 2, joint)

fits0Frozen = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s0Frozen]
fits2Frozen = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s2Frozen]

fits0Deflected = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s0Deflected]
fits2Deflected = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s2Deflected]

fits0Thermal = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s0Thermal]
fits2Thermal = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s2Thermal]

fits0Joint = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s0Joint]
fits2Joint = [power_fit(1 ./ x[1][2:end], -log10.(x[2][2:end])) for x in s2Joint]
