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

@info 2Ek_MnF([0pi / p.a, 0pi / p.a, 0pi / p.a], p)[1]
@info 2sqrt(p.K * (p.K + 2 * 8 * p.J))

@info 2Ek_MnF([pi / p.a, pi / p.a, pi / p.a], p)[1]
@info 2 * (p.K + 8p.J)
##


Nq_plasmons = 200
qrange_plasmons = range(1e-5u"nm^-1", stop=0.01u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
NE = 200
Erange = range(0.1u"meV", stop=20u"meV", length=NE)
spectralunit = u"meV^-1"
res_Γ_y = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
res_Γ_x = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)


@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_MnF([q, 0.0u"nm^-1", 0.0u"nm^-1"], E, p, pp)[1] |> u"meV"
        res_Γ_x[i, j] = Γloc

    end
end
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_MnF([0.0u"nm^-1", q, 0.0u"nm^-1"], E, p, pp)[1] |> u"meV"
        res_Γ_y[i, j] = Γloc
    end
end

##
res_y = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_x = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)

Nz = 1e3
@info Nz * p.c |> u"cm"
correction = q -> 1 / (-1 + exp(2κ(q, pp) * p.c))
# correction = q -> (1 - exp(-2Nz * κ(q, pp) * p.c)) / (-1 + exp(2κ(q, pp) * p.c))
ratio_x = zeros(Float64, NE, Nq_plasmons)
ratio_y = zeros(Float64, NE, Nq_plasmons)


for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        # Γloc = Nz * res_Γ_y[i, j]
        Γloc = correction(q) * res_Γ_y[i, j]
        res_y[i, j] = spectral_plasmons(E, q, pp, Γloc)
        ratio_y[i, j] = abs(Γloc) / E

    end
end

for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        # Γloc = Nz * res_Γ_x[i, j]
        Γloc = correction(q) * res_Γ_x[i, j]
        res_x[i, j] = spectral_plasmons(E, q, pp, Γloc)
        ratio_x[i, j] = abs(Γloc) / E

    end
end
@info maximum(ratio_x)
@info qrange_plasmons[argmax(ratio_x)[2]] |> u"μm^-1"
@info Erange[argmax(ratio_x)[1]] |> u"meV"

@info maximum(ratio_y)

fig, ax = plt.subplots()

fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum([abs.(res_x) |> ustrip; abs.(res_y) |> ustrip]))

# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)
cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, abs.(res_x) |> ustrip, norm=fnorm, shading="nearest", linewidth=0, rasterized=true, cmap="Blues_r")
cb.set_edgecolor("face")
cb = ax.pcolormesh(-qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, abs.(res_y) |> ustrip, norm=fnorm, shading="nearest", linewidth=0, rasterized=true, cmap="Blues_r")
cb.set_edgecolor("face")

ax.text(-0.1, -0.2, L"\mathbf q\parallel [001]", transform=ax.transAxes, size=8)
ax.text(0.97, -0.2, L"\mathbf q\parallel [0\bar10]", transform=ax.transAxes, size=8)


clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

top_E = Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[1] + Ek_MnF([pi / p.a, 0pi / p.a, pi / p.a], p)[2]
# top_q = qrange_plasmons[argmin(abs.(top_E .- Eplasmons_res))]
min_E = 2Ek_MnF([0pi / p.a, 0pi / p.a, 0pi / p.a], p)[1]
# min_q = qrange_plasmons[argmin(abs.(min_E .- Eplasmons_res))]
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)
ax.axhline(min_E |> u"meV" |> ustrip, color="white", lw=0.5)

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, color="white", linestyle="dotted", lw=1, alpha=1)
Eplasmons_res = [Eq_plasmons([-0q, -q], pp) .|> u"meV" for q in qrange_plasmons]
ax.plot(-qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, linestyle="dotted", color="white", lw=1, alpha=1)

ax.xaxis.set_major_formatter(plt.matplotlib.ticker.FuncFormatter((x, pos) -> @sprintf "%.0f" abs(x)))
ax.set_ylim(extrema(Erange .|> u"meV" .|> ustrip)...)
ax.set_xlabel(L"$q$ ($\mathrm{\mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "MnF_3D_spectral_inf.pdf"), fig)
fig

##