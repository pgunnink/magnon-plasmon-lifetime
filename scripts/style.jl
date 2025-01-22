using PyPlot
using DrWatson
plt.style.use(joinpath(@__DIR__, "prl.mplstyle"))

save_dir = joinpath(papersdir(), "figures/")
column_width = 3 + 3 / 8
fig_width = column_width
fig_height = fig_width / 1.618
