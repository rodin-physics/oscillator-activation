using ForwardDiff
using NLsolve
include("plotting_settings.jl")

set_theme!(publication_theme)

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

## PLOT PARAMETERS
colors = [my_vermillion, my_sky]
impact_labs = ["Head-on", "Offset"]
symbols = ['◯', '×', '+']
labs_mass = ["Light", "Medium", "Heavy"]

sizes = [18, 36, 36]

nPts = 500
ωTs = range(s0Frozen[1][1][1], s0Frozen[1][1][end], length = nPts)
# Interaction function
Φ0 = 10
λ = 2
@inline function Φ(r)
    res = Φ0 * exp(-r^2 / 2 / λ^2)
    return res
end

group_marker = [
    MarkerElement(
        marker = symbols[ii],
        color = :black,
        strokecolor = :transparent,
        markersize = sizes[ii],
    ) for ii = 1:3
]
group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors]
group_activation =
    [LineElement(color = :black, linewidth = 4, linestyle = s) for s in [nothing, :dash]]

# STATIONARY OSCILLATOR
fig = Figure(resolution = (1200, 1800))
gdata = fig[1, 1] = GridLayout()
glegend = fig[2, 1] = GridLayout()

axtop = Axis(
    gdata[1, 1],
    ylabel = "Transmission probability",
    xlabel = L"\omega_T",
    title = "Frozen oscillator",
    titlefont = :boldlatex,
)

axbottom = Axis(
    gdata[2, 1],
    ylabel = "Transmission probability",
    xlabel = L"\omega_T",
    title = "Deflected oscillator",
    titlefont = :boldlatex,
)

text!(
    axtop,
    0.05,
    0.95,
    text = L"(\mathrm{a})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)
text!(
    axbottom,
    0.05,
    0.95,
    text = L"(\mathrm{b})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)

rowsize!(fig.layout, 2, Auto(1 / 8))

insettop =
    Axis(fig, bbox = BBox(650, 1130, 1210, 1570), xlabel = L"1/\omega_T", yscale = log10)
insetbottom =
    Axis(fig, bbox = BBox(650, 1130, 410, 770), xlabel = L"1/\omega_T", yscale = log10)

ylims!(axtop, (-0.05, 1.05))
ylims!(axbottom, (-0.05, 1.05))

ylims!(insettop, (1e-3, 2))
ylims!(insetbottom, (1e-3, 2))

Legend(
    glegend[1, 1],
    [group_marker, group_color, group_activation],
    [labs_mass, impact_labs, ["Frozen", "Relaxed"]],
    ["Mass", "Impact", "Activation"],
    halign = :center,
    valign = :center,
    tellheight = false,
    tellwidth = false,
    # margin = (10, 10, 10, 10),
    framevisible = false,
    titlefont = :boldlatex,
    orientation = :horizontal,
)

for ii in eachindex(s0Frozen)
    d = s0Frozen[ii]
    scatter!(
        axtop,
        d[1],
        d[2],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insettop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Frozen)
    d = s2Frozen[ii]
    scatter!(
        axtop,
        d[1],
        d[2],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insettop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s0Deflected)
    d = s0Deflected[ii]
    scatter!(
        axbottom,
        d[1],
        d[2],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insetbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Deflected)
    d = s2Deflected[ii]
    scatter!(
        axbottom,
        d[1],
        d[2],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insetbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

lines!(axtop, ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axtop, ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)

lines!(insettop, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(insettop, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)


function fn(ρ)
    der = ForwardDiff.derivative(x -> Φ(x), ρ[1] - 2)
    return [-ρ[1] - der]
end

sol = (nlsolve(fn, [0.0]).zero)[1]
activation = Φ(sol - 2) + 1 * sol^2 / 2

lines!(axbottom, ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axbottom, ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)
lines!(
    axbottom,
    ωTs,
    exp.(-activation ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(insetbottom, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(insetbottom, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)
lines!(
    insetbottom,
    1 ./ ωTs,
    exp.(-activation ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

fig

save("Stationary_Oscillator.pdf", fig)

# MOVING OSCILLATOR
fig = Figure(resolution = (1200, 1800))
gdata = fig[1, 1] = GridLayout()
glegend = fig[2, 1] = GridLayout()

axtop = Axis(
    gdata[1, 1],
    ylabel = "Transmission probability",
    xlabel = L"\omega_T",
    title = "Separate distributions",
    titlefont = :boldlatex,
)

axbottom = Axis(
    gdata[2, 1],
    ylabel = "Transmission probability",
    xlabel = L"\omega_T",
    title = "Joint distributions",
    titlefont = :boldlatex,
)

text!(
    axtop,
    0.05,
    0.95,
    text = L"(\mathrm{a})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)
text!(
    axbottom,
    0.05,
    0.95,
    text = L"(\mathrm{b})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)

rowsize!(fig.layout, 2, Auto(1 / 8))

insettop =
    Axis(fig, bbox = BBox(650, 1130, 1210, 1570), xlabel = L"1/\omega_T", yscale = log10)
insetbottom =
    Axis(fig, bbox = BBox(650, 1130, 410, 770), xlabel = L"1/\omega_T", yscale = log10)

ylims!(axtop, (-0.05, 1.05))
ylims!(axbottom, (-0.05, 1.05))

ylims!(insettop, (1e-3, 2))
ylims!(insetbottom, (1e-3, 2))

Legend(
    glegend[1, 1],
    [group_marker, group_color, group_activation],
    [labs_mass, impact_labs, ["Frozen", "Relaxed"]],
    ["Mass", "Impact", "Activation"],
    halign = :center,
    valign = :center,
    tellheight = false,
    tellwidth = false,
    # margin = (10, 10, 10, 10),
    framevisible = false,
    titlefont = :boldlatex,
    orientation = :horizontal,
)

for ii in eachindex(s0Thermal)
    d = s0Thermal[ii]
    scatter!(
        axtop,
        d[1],
        d[2],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insettop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Thermal)
    d = s2Thermal[ii]
    scatter!(
        axtop,
        d[1],
        d[2],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insettop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s0Joint)
    d = s0Joint[ii]
    scatter!(
        axbottom,
        d[1],
        d[2],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insetbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Joint)
    d = s2Joint[ii]
    scatter!(
        axbottom,
        d[1],
        d[2],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        insetbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

function fn(ρ)
    der = ForwardDiff.derivative(x -> Φ(x), ρ[1] - 2)
    return [-ρ[1] - der]
end

sol = (nlsolve(fn, [0.0]).zero)[1]
activationOffset = Φ(sol - 2) + 1 * sol^2 / 2

function fn(ρ)
    der = ForwardDiff.derivative(x -> Φ(x), ρ[1])
    return [-ρ[1] - der]
end

sol = (nlsolve(fn, [-4.0]).zero)[1]
activationHeadOn = Φ(sol) + 1 * sol^2 / 2

lines!(
    axtop,
    ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    axtop,
    ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)
lines!(axtop, ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axtop, ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)

lines!(
    insettop,
    1 ./ ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    insettop,
    1 ./ ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(insettop, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(insettop, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)

lines!(
    axbottom,
    ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    axbottom,
    ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(axbottom, ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axbottom, ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)


lines!(
    insetbottom,
    1 ./ ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    insetbottom,
    1 ./ ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(insetbottom, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(insetbottom, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)
fig

save("Moving_Oscillator.pdf", fig)

## PHASE SPACE

## LOAD THE DATA
phase = readdir("Data/Phase/"; join = true)
phase = load_object.(phase)

impact_parameters = unique([x[end-1] for x in phase])
masses = unique([x[end] for x in phase])

s0Phase = filter(x -> x[end-1] == 0, phase)
s2Phase = filter(x -> x[end-1] == 2, phase)
# using GLMakie
# GLMakie.activate!()

# d = s0Phase[1]
# cmap = :Hiroshige
# volume(
#     d[1],
#     d[2],
#     d[3],
#     d[4],
#     algorithm = :iso,
#     isorange = 0.05,
#     isovalue = 0,
#     colormap = cmap,
#     alpha = 0.8
# )
# # volume()

# r = LinRange(-1, 1, 100)
# cube = [(x .^ 2 + y .^ 2 + z .^ 2) for x in r, y in r, z in r]
# contour(cube, alpha = 0.5)
