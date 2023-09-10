using CairoMakie
using Colors
using JLD2

## PLOTTING 
## Friendly colors
my_red = colorant"rgba(204, 121, 167, 1.0)"
my_vermillion = colorant"rgba(213, 94, 0, 1.0)"
my_orange = colorant"rgba(230, 159, 0, 1.0)"
my_yellow = colorant"rgba(240, 228, 66, 1.0)"
my_green = colorant"rgba(0, 158, 115, 1.0)"
my_sky = colorant"rgba(86, 180, 233, 1.0)"
my_blue = colorant"rgba(0, 114, 178, 1.0)"
my_black = colorant"rgba(0, 0, 0, 1.0)"
set_theme!()
Makie.to_font("CMU Serif Bold")

publication_theme = Theme(
    fonts = (; latex = "CMU Serif", boldlatex = "CMU Serif Bold"),
    fontsize = 36,
    figure_padding = 10,
    # palette = (
    #     color = [
    #         my_red,
    #         my_vermillion,
    #         my_orange,
    #         my_yellow,
    #         my_green,
    #         my_sky,
    #         my_blue,
    #         my_black,
    #     ]
    # ),
    Axis = (
        titlefont = :latex,
        xgridvisible = false,
        ygridvisible = false,
        xlabelfont = :latex,
        ylabelfont = :latex,
        xticklabelfont = :latex,
        yticklabelfont = :latex,
        xticklabelsize = 36,
        yticklabelsize = 36,
        titlesize = 36,
        xlabelsize = 36,
        ylabelsize = 36,
    ),
    Legend = (labelsize = 36, labelfont = :latex),
)
