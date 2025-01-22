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

Nq_plasmons = 200
qrange_plasmons = range(1e-5u"nm^-1", stop=5u"μm^-1", length=Nq_plasmons)

# qrange_plasmons = range(1e-5u"nm^-1", stop=0.05u"nm^-1", length=Nq_plasmons)

@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
# @info» 1 / κ(qrange_plasmons[end], pp) |> u"Å"
##
NE = 200
Erange = range(0.1u"meV", stop=15u"meV", length=NE)
# Erange = range(1u"meV", stop=20u"meV", length=NE)

spectralunit = u"meV^-1"
res_Γ_right = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
res_Γ_left = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)


res_Γ_right_backwards = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
res_Γ_left_backwards = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)


@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_cupr([-0q, -q], E, p, pp)[1] |> u"meV"
        res_Γ_left[i, j] = Γloc

        Γloc = Γ_cupr_backwards([-0q, -q], E, p, pp)[1] |> u"meV"
        res_Γ_left_backwards[i, j] = Γloc

    end
end
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_cupr([q, 0q], E, p, pp)[1] |> u"meV"
        res_Γ_right[i, j] = Γloc

        Γloc = Γ_cupr_backwards([q, 0q], E, p, pp)[1] |> u"meV"
        res_Γ_right_backwards[i, j] = Γloc

    end
end

##
res_right = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_left = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)

Nz = 1e3
@info Nz * p.c |> u"cm"
# correction = q -> sum(exp(-2κ(q, pp) * n * lattice_c) for n in 1:Nz)
# correction = q -> (1 - exp(-2Nz * κ(q, pp) * p.c)) / (-1 + exp(2κ(q, pp) * p.c))
correction = q -> 1 / (-1 + exp(2κ(q, pp) * p.c))
ratio_left = zeros(Float64, NE, Nq_plasmons)
ratio_right = zeros(Float64, NE, Nq_plasmons)

for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = correction(q) * (res_Γ_right[i, j])# + res_Γ_right_backwards[i, j])
        res_right[i, j] = spectral_plasmons(E, q, pp, Γloc)
        ratio_right[i, j] = abs(Γloc) / E
    end
end

for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = correction(q) * res_Γ_left[i, j]
        res_left[i, j] = spectral_plasmons(E, q, pp, Γloc)
        ratio_left[i, j] = abs(Γloc) / E
    end
end
@info maximum(ratio_left)
@info maximum(ratio_right)
fig, ax = plt.subplots()

fnorm = plt.matplotlib.colors.LogNorm(5e-2, maximum([abs.(res_left) |> ustrip; abs.(res_right) |> ustrip]))
# fnorm = plt.matplotlib.colors.Normalize(1e-2, maximum([abs.(res_left) |> ustrip; abs.(res_right) |> ustrip]))

# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)
cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, abs.(res_right) |> ustrip, norm=fnorm, shading="nearest", linewidth=0, rasterized=true, cmap="Blues_r")
cb.set_edgecolor("face")
cb = ax.pcolormesh(-qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, abs.(res_left) |> ustrip, norm=fnorm, shading="nearest", linewidth=0, rasterized=true, cmap="Blues_r")
cb.set_edgecolor("face")

ax.text(-0.1, -0.2, L"\mathbf q\parallel \hat{y}", transform=ax.transAxes, size=8)
ax.text(0.97, -0.2, L"\mathbf q\parallel \hat{x}", transform=ax.transAxes, size=8)


clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")

top_E = 2Ek_cupr([pi / p.a, 0pi / p.a], p)
# top_q = qrange_plasmons[argmin(abs.(top_E .- Eplasmons_res))]
min_E = 2Ek_cupr([0pi / p.a, 0pi / p.a], p)
# min_q = qrange_plasmons[argmin(abs.(min_E .- Eplasmons_res))]
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)
ax.axhline(min_E |> u"meV" |> ustrip, color="white", lw=0.5)

Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, color="white", linestyle="dotted", lw=1, alpha=1)
Eplasmons_res = [Eq_plasmons([-0q, -q], pp) .|> u"meV" for q in qrange_plasmons]
ax.plot(-qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, linestyle="dotted", color="white", lw=1, alpha=1)
ax.xaxis.set_major_formatter(plt.matplotlib.ticker.FuncFormatter((x, pos) -> @sprintf "%.0f" abs(x)))

ax.set_ylim(extrema(Erange .|> u"meV" .|> ustrip)...)
ax.set_xlabel(L"$q$ ($\mathrm{ \mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "cuprate_spectral_inf.pdf"), fig)
fig
