using ColorSchemes
using ForwardDiff
using JLD2
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
fig = Figure(size = (1200, 1800))
gdata = fig[1, 1] = GridLayout()
glegend = fig[2, 1] = GridLayout()

axtop = Axis(
    gdata[1, 1],
    ylabel = "Transmission probability",
    xlabel = L"1/\omega_T",
    title = "Frozen oscillator",
    titlefont = :boldlatex,
    yscale = log10,
)

axbottom = Axis(
    gdata[2, 1],
    ylabel = "Transmission probability",
    xlabel = L"1/\omega_T",
    title = "Deflected oscillator",
    titlefont = :boldlatex,
    yscale = log10,
)

text!(
    axtop,
    0.85,
    0.95,
    text = L"(\mathrm{a})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)
text!(
    axbottom,
    0.85,
    0.95,
    text = L"(\mathrm{b})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)

rowsize!(fig.layout, 2, Auto(1 / 8))

insettop = Axis(
    fig,
    bbox = BBox(150, 510, 1110, 1380),
    xlabel = L"\omega_T",
    xaxisposition = :top,
    yaxisposition = :right,
)

ylims!(insettop, (-0.05, 1.05))

ylims!(axtop, (1e-4, 2))
ylims!(axbottom, (1e-4, 2))

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
        insettop,
        d[1],
        d[2],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii] / 2,
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axtop,
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
        insettop,
        d[1],
        d[2],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii] / 2,
    )
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axtop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s0Deflected)
    d = s0Deflected[ii]
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Deflected)
    d = s2Deflected[ii]
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

lines!(insettop, ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(insettop, ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)

lines!(axtop, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axtop, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)


function fn(ρ)
    der = ForwardDiff.derivative(x -> Φ(x), ρ[1] - 2)
    return [-ρ[1] - der]
end

sol = (nlsolve(fn, [0.0]).zero)[1]
activation = Φ(sol - 2) + 1 * sol^2 / 2

lines!(axbottom, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axbottom, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)
lines!(
    axbottom,
    1 ./ ωTs,
    exp.(-activation ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

fig

save("Stationary_Oscillator.pdf", fig)

# MOVING OSCILLATOR
fig = Figure(size = (1200, 1800))
gdata = fig[1, 1] = GridLayout()
glegend = fig[2, 1] = GridLayout()

axtop = Axis(
    gdata[1, 1],
    ylabel = "Transmission probability",
    xlabel = L"1/\omega_T",
    title = "Separate distributions",
    titlefont = :boldlatex,
    yscale = log10,
)

axbottom = Axis(
    gdata[2, 1],
    ylabel = "Transmission probability",
    xlabel = L"1/\omega_T",
    title = "Joint distributions",
    titlefont = :boldlatex,
    yscale = log10,
)

text!(
    axtop,
    0.85,
    0.95,
    text = L"(\mathrm{a})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)
text!(
    axbottom,
    0.85,
    0.95,
    text = L"(\mathrm{b})",
    align = (:left, :top),
    space = :relative,
    fontsize = 36,
)

rowsize!(fig.layout, 2, Auto(1 / 8))

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
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axtop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Thermal)
    d = s2Thermal[ii]
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axtop,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[2],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s0Joint)
    d = s0Joint[ii]
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axbottom,
        1 ./ d[1][idx],
        d[2][idx],
        color = colors[1],
        marker = symbols[ii],
        markersize = sizes[ii],
    )

end

for ii in eachindex(s2Joint)
    d = s2Joint[ii]
    idx = findall(x -> x > 0, d[2])
    scatter!(
        axbottom,
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
    1 ./ ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    axtop,
    1 ./ ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(axtop, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axtop, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)

lines!(
    axbottom,
    1 ./ ωTs,
    exp.(-activationHeadOn ./ ωTs),
    color = colors[1],
    linewidth = 4,
    linestyle = :dash,
)
lines!(
    axbottom,
    1 ./ ωTs,
    exp.(-activationOffset ./ ωTs),
    color = colors[2],
    linewidth = 4,
    linestyle = :dash,
)

lines!(axbottom, 1 ./ ωTs, exp.(-Φ(0) ./ ωTs), color = colors[1], linewidth = 4)
lines!(axbottom, 1 ./ ωTs, exp.(-Φ(2) ./ ωTs), color = colors[2], linewidth = 4)
fig

save("Moving_Oscillator.pdf", fig)

## PHASE SPACE

## LOAD THE DATA

s0Phase = load_object("Data/Phase/Phase_s0_Φ10_λ2_μ1.jld2")
s2Phase = load_object("Data/Phase/Phase_s2_Φ10_λ2_μ1.jld2")
my_cmap = ColorScheme(range(my_vermillion, my_blue, length = 2))

fig = Figure(size = (2400, 600))

ϕs_lab = [
    L"\phi = 0",
    L"\pi / 4",
    L" \pi / 2",
    L"3\pi / 4",
    L"\pi",
    L"5 \pi / 4",
    L"3\pi / 2",
    L"7 \pi / 4",
]

main_grid = fig[1, 1] = GridLayout()

gs0 = main_grid[1, 1] = GridLayout()
gs2 = main_grid[2, 1] = GridLayout()

axtop = [Axis(gs0[1, ii], xaxisposition = :top) for ii in eachindex(ϕs_lab)]
axbottom = [Axis(gs2[1, ii]) for ii in eachindex(ϕs_lab)]

lab_top = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
lab_bottom = ["(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)"]

for ii in eachindex(ϕs_lab)
    contourf!(axtop[ii], s0Phase[1], s0Phase[2], s0Phase[4][:, :, ii], colormap = my_cmap)
    contourf!(
        axbottom[ii],
        s2Phase[1],
        s2Phase[2],
        s2Phase[4][:, :, ii],
        colormap = my_cmap,
    )

    axtop[ii].xticks = [0, 5, 10]
    axbottom[ii].xlabel = ϕs_lab[ii]

    hidexdecorations!(axbottom[ii], label = false)

    text!(
        axtop[ii],
        0.75,
        0.95,
        text = lab_top[ii],
        align = (:left, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )

    text!(
        axbottom[ii],
        0.75,
        0.95,
        text = lab_bottom[ii],
        align = (:left, :top),
        space = :relative,
        fontsize = 36,
        font = :latex,
        color = :black,
    )
end

for ii in eachindex(ϕs_lab)[1:end-1]
    hideydecorations!(axtop[ii], label = false)
    hideydecorations!(axbottom[ii], label = false)
end
axtop[1].ylabel = "Head-on"
axbottom[1].ylabel = "Offset"

axtop[1].xlabel = L"\mathcal{E}_\mathrm{particle}"
axtop[end].ylabel = L"\mathcal{E}\mathrm{osc}"
axbottom[end].ylabel = ""

axtop[end].yaxisposition = :right
axbottom[end].yaxisposition = :right

axtop[end].yticks = [0, 10, 20]
axbottom[end].yticks = [0, 10, 20]

fig
save("Phase_Space.pdf", fig)
