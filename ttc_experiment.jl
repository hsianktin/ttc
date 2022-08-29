# a schematic model for the LacZ complemention system

T = 10 # delay
f(x) = exp(-(x - T)^2/2)/√(2π)
using NumericalIntegration
function ∫fdx(f,x,dx=1e-1)
    lin_space = collect(0:dx:x)
    if length(lin_space) > 1
    return integrate(lin_space,f.(lin_space))
    else
        return 0
    end
end
function ∫gdx(f,x,dx=1e-1)
    lin_space = collect(0:dx:x)
    if length(lin_space) > 1
    return sum(f.(lin_space))*dx
    else
        return 0
    end
end
F(x) = ∫fdx(f,x)
G(x) = ∫fdx(F,x,1e-1)

X = collect(0:0.1:20)
fX = f.(X)
FX = F.(X)
@time GX = G.(X)

using DataFrames
plot_df = DataFrame(
    x = X,
    f = fX,
    F = FX,
    G = GX
)
CSV.write("fig/sample_df.csv", plot_df)
h(x) = exp(-(x - T)^2)/√(π)
g(x) = (h(x) + h(x-5))/2
Fg(x) = ∫fdx(g, x)
FG(x) = ∫fdx(Fg, x)
plot_df = DataFrame(
    x = X,
    f = g.(X),
    F = Fg.(X),
    G = FG.(X)
)
CSV.write("fig/sample_distinct_bimodal.csv", plot_df)

g(x) = (f(x) + f(x+3))/2
Fg(x) = ∫fdx(g, x) 
FG(x) = ∫fdx(Fg, x,0.1)

plot_df = DataFrame(
    x = X,
    f = g.(X)./Fg(20),
    F = Fg.(X)./ ∫fdx(g,20),
    G = FG.(X)./ ∫fdx(g,20),
)
CSV.write("fig/sample_adjacent_bimodal.csv", plot_df)

using PGFPlotsX
axis1 = @pgf Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "time \$t\$",
        ylabel = "linearized signal \$S(t)\$",
        # grid = "major",
        legend_pos  = "south east"
    },
    Plot(
        Table({x = "x", y = "G", "col sep"="comma"}, "sample_df.csv")
        ),
    Plot(
        Table({x = "x", y = "G", "col sep"="comma"}, "sample_distinct_bimodal.csv")
        ),
    Plot(
        Table({x = "x", y = "G", "col sep"="comma"}, "sample_adjacent_bimodal.csv")
        ),
    )
    
axis2 = @pgf Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "time \$t\$",
        ylabel = " \$\\mathrm{Pr}(T \\leq t)\$",
        # grid = "major",
        legend_pos  = "south east"
    },
    Plot(
        Table({x = "x", y = "F", "col sep"="comma"}, "sample_df.csv")
        ),
    Plot(
        Table({x = "x", y = "F", "col sep"="comma"}, "sample_distinct_bimodal.csv")
        ),
    Plot(
        Table({x = "x", y = "F", "col sep"="comma"}, "sample_adjacent_bimodal.csv")
        ),
)
axis3 = @pgf Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "time \$t\$",
        ylabel = " \$\\rho_T(t)\$",
        # grid = "major",
        legend_pos  = "south east"
    },
    Plot(
        Table({x = "x", y = "f", "col sep"="comma"}, "sample_df.csv")
        ),
    Plot(
        Table({x = "x", y = "f", "col sep"="comma"}, "sample_distinct_bimodal.csv")
        ),

    Plot(
        Table({x = "x", y = "f", "col sep"="comma"}, "sample_adjacent_bimodal.csv")
        ),
)

pgfsave("fig/sample_plot_1.tex", axis1)
pgfsave("fig/sample_plot_2.tex", axis2)
pgfsave("fig/sample_plot_3.tex", axis3)
