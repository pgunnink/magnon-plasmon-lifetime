using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
includet("cuprate_def.jl")

##
p = ParamsCuprate()
pp = ParamsPlasmons()
@info "kF" pp.Ef / (ħ * 1e6u"m/s") |> u"nm^-1"
@info "edge of particle-hole " ħ * ωq_plasmons(pp.Ef / (ħ * 1e6u"m/s"), pp) |> u"meV"
kunit_display = u"μm^-1"
Γ_display = u"meV"
dosunit = u"meV^-1"
Eunit = u"meV"

@info "bottom of band" Ek_cupr([0pi / p.a, 0pi / p.a], p)
lim_min = 1e-3u"meV"
lim_max = 2.5Ek_cupr([pi / p.a, 0pi / p.a], p)

##
NE = 200
Nk_dos = 200
end_point = pi / p.a

krange_dos_right = range(0 / p.a, stop=end_point, length=Nk_dos / 2 |> Int)
krange_dos_left = range(-end_point, stop=0 / p.a, length=Nk_dos / 2 |> Int)

Erange = range(lim_min, stop=lim_max, length=NE)
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos_right)
        res_loc = dos2_adaptive([k, 0k], E, x -> Ek_cupr(x, p), p.a; rtol=1e-3, maxevals=Int(1e5), N=1, d=2, η=1e-2u"meV")
        res[i, j+Nk_dos/2|>Int] = res_loc[1]
    end

    for (j, k) in enumerate(krange_dos_left)
        res_loc = dos2_adaptive([0k, k], E, x -> Ek_cupr(x, p), p.a; rtol=1e-3, maxevals=Int(1e5), N=1, d=2, η=1e-2u"meV")
        res[i, j] = res_loc[1]
    end
end

##
NE_cut = 100
Erange_cut = range(lim_min, stop=lim_max, length=NE_cut)
res_cut = zeros(typeof(1.0dosunit), NE_cut)

@showprogress for (i, E) in enumerate(Erange_cut)
    res_cut[i] = dos2_adaptive([0 / p.a, 0 / p.a, 0 / p.a], E, x -> Ek_cupr(x, p), p.a; rtol=1e-5, maxevals=Int(1e5), N=1, d=3, η=1e-1u"meV")[1]
end

# ##
# NE_cut = 100
# Erange_cut = range(lim_min, stop=lim_max, length=NE_cut)
# res_cut_right = zeros(typeof(1.0dosunit), NE_cut)
# res_cut_left = zeros(typeof(1.0dosunit), NE_cut)
# E_cut_right = zeros(typeof(1.0u"meV"), NE_cut)
# E_cut_left = zeros(typeof(1.0u"meV"), NE_cut)


# # @showprogress for (i, E) in enumerate(Erange_cut)
# @showprogress for (i, k) in enumerate(krange_plasmons)
#     res_cut_left[i] = dos2_adaptive([-0k, -k], ħ * ωq_plasmons([-0k, -k], pp), x -> Ek_cupr(x, p), p.a; rtol=1e-5, maxevals=Int(1e5), N=1, d=2, η=1e-2u"meV")[1]
#     E_cut_left[i] = ħ * ωq_plasmons([-0k, -k], pp)
#     res_cut_right[i] = dos2_adaptive([k, 0k], ħ * ωq_plasmons([k, 0k], pp), x -> Ek_cupr(x, p), p.a; rtol=1e-5, maxevals=Int(1e5), N=1, d=2, η=1e-2u"meV")[1]
#     E_cut_right[i] = ħ * ωq_plasmons([k, 0k], pp)

# end

##
krange_plasmons = range(1e-5u"μm^-1", 5u"μm^-1", length=NE_cut)

fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), width_ratios=[3, 1], height_ratios=[1, 20])
ax = axs[2, 1]
Nk = 100
krange = range(0 / p.a, stop=end_point, length=Nk)
to_plot = [Ek_cupr([k, 0k], p) .|> Eunit .|> ustrip for k in krange]
ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color="white")

krange = range(0 / p.a, stop=-end_point, length=Nk)
to_plot = [Ek_cupr([0k, -k], p) .|> Eunit .|> ustrip for k in krange]
ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color="white")

toplot = res .|> dosunit |> ustrip
@info maximum(toplot)
fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(1, maximum(toplot))

cb = ax.pcolormesh([krange_dos_left; krange_dos_right] ./ (end_point) .|> upreferred .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="nearest", linewidth=0, rasterized=true)
cb.set_edgecolor("face")
clb = fig.colorbar(cb, cax=axs[1, 1], orientation="horizontal", location="top")
clb.ax.set_title(L"$\mathrm{D}_{\mathbf{q}}$ ($ \mathrm{meV}^{-1}$)")
# ax.set_xlabel(L"k")
ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
ax.set_ylabel(L"$\hbar\omega$ (meV)")
ax.set_xlim(upreferred.(extrema(krange ./ (pi / p.a)))...)
ax.set_xticks([-1, 0, 1], [L"(0,\pi/a)", L"\Gamma", L"(\pi/a,0)"])

ax.plot(krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, 0k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--", lw=1)

ax.plot(-krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([-k, -k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--", lw=1)
# ax.axhline(2Ek_cupr([0pi / p.a, pi / p.a], p) |> u"meV" |> ustrip, color="black", linestyle="dashed")

ax = axs[2, 2]
ax.plot(res_cut .|> dosunit |> ustrip, Erange_cut .|> u"meV" .|> ustrip,)
# ax.plot(res_cut_right .|> dosunit |> ustrip, E_cut_right .|> u"meV" .|> ustrip, label=L"\mathbf{q}\parallel(\hat{x})")
ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
# ax.axhline(2Ek_cupr([0pi / p.a, pi / p.a], p) |> u"meV" |> ustrip, color="black", linestyle="dashed")
# ax.set_xscale("log")
ax.set_yticklabels([])
# ax.set_yticks([])
ax.set_xlabel(L"$\mathrm{D}_{\mathbf{q}}$ $(\mathrm{meV}^{-1}$)")
ax.margins(y=0)
ax.xaxis.set_label_position("top")
ax.xaxis.tick_top()

# ax.legend(loc="upper center", fontsize=5)
axs[1, 2].axis("off")
# ax.tick_params(which="both", axis="y", direction="in", pad=-22)
fig.subplots_adjust(wspace=0.05, hspace=0.05)
wsave(joinpath(save_dir, "dispersion_cuprate.pdf"), fig)
fig |> display
##

