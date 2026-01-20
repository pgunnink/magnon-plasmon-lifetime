using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
includet("MnF2_3D_def.jl")

##
p = ParamsMnF()
pp = ParamsPlasmons()


kunit_display = u"μm^-1"
Γ_display = u"meV"

lim_min = 1e-3u"meV"
lim_max = 1.5Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[1] + 1.5Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[2]


##

NE = 200
Nk_dos = 200
end_point = pi / p.a

krange_dos_right = range(0 / p.a, stop=end_point, length=Nk_dos / 2 |> Int)
krange_dos_left = range(-end_point, stop=0 / p.a, length=Nk_dos / 2 |> Int)

Erange = range(lim_min, stop=lim_max, length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos_right)
        res_loc = dos2_adaptive([k, 0k, 0k], E, x -> Ek_MnF(x, p), p.a; rtol=1e-3, maxevals=Int(1e5), N=2, d=3, η=1e-1u"meV")
        res[i, j+Nk_dos/2|>Int] = res_loc[1]
    end

    for (j, k) in enumerate(krange_dos_left)
        res_loc = dos2_adaptive([0k, k, 0k], E, x -> Ek_MnF(x, p), p.a; rtol=1e-3, maxevals=Int(1e5), N=2, d=3, η=1e-1u"meV")
        res[i, j] = res_loc[1]
    end
end
NE_cut = 100
Erange_cut = range(lim_min, stop=lim_max, length=NE_cut)
res_cut = zeros(typeof(1.0dosunit), NE_cut)

@showprogress for (i, E) in enumerate(Erange_cut)
    res_cut[i] = dos2_adaptive([0 / p.a, 0 / p.a, 0 / p.a], E, x -> Ek_MnF(x, p), p.a; rtol=1e-5, maxevals=Int(1e5), N=2, d=3, η=1e-1u"meV")[1]
end

##
fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), width_ratios=[3, 1], height_ratios=[1, 20])
ax = axs[2, 1]
Eunit = u"meV"
Nk = 100
krange = range(0 / p.a, stop=end_point, length=Nk)

splitting_colors = ["lightcoral", "deepskyblue"]
splitting_colors = ["white", "white"]

for i in 1:1
    to_plot = [Ek_MnF([k, 0k, 0k], p)[i] .|> Eunit .|> ustrip for k in krange]
    ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color=splitting_colors[i])
end

krange = range(0 / p.a, stop=-end_point, length=Nk)
for i in 1:1
    to_plot = [Ek_MnF([0k, k, 0k], p)[i] .|> Eunit .|> ustrip for k in krange]
    ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color=splitting_colors[i])
end
toplot = res .|> dosunit |> ustrip
@info maximum(toplot)
fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(1e-2, maximum(toplot))

cb = ax.pcolormesh([krange_dos_left; krange_dos_right] ./ (end_point) .|> upreferred .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="nearest", linewidth=0, rasterized=true)
cb.set_edgecolor("face")
clb = fig.colorbar(cb, cax=axs[1, 1], orientation="horizontal", location="top")
clb.ax.set_title(L"$D_{\mathbf{q}}(\epsilon)$ $(\mathrm{meV}^{-1}$)")
ax.set_xlabel(L"k")
ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
ax.set_ylabel(L"$\hbar\omega$ (meV)")
ax.set_xlim(upreferred.(extrema(krange ./ (pi / p.a)))...)
ax.set_xticks([-1, 0, 1], [L"(0,\pi/a,0)", L"\Gamma", L"(\pi/a,0,0)"])
krange_plasmons = range(-0.1u"nm^-1", 0.1u"nm^-1", length=5000)
ax.plot(krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, 0k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--", lw=1)
# ax.axhline((Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[1] + Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[2]) |> u"meV" |> ustrip, color="black", linestyle="dashed")
ax = axs[2, 2]
ax.plot(res_cut .|> dosunit |> ustrip, Erange_cut .|> u"meV" .|> ustrip)
# ax.set_xscale("log")
ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xlabel(L"$D_{\mathbf{q}}(\epsilon)$ $(\mathrm{meV}^{-1}$)")
ax.margins(y=0)
ax.xaxis.set_label_position("top")
ax.xaxis.tick_top()

axs[1, 2].axis("off")
# ax.tick_params(which="both", axis="y", direction="in", pad=-22)
fig.subplots_adjust(wspace=0.05, hspace=0.05)
wsave(joinpath(save_dir, "dispersion_MnF_3D.pdf"), fig)
fig |> display

