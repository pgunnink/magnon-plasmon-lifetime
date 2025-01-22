using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
include("style.jl")
includet("EP_functions.jl")
Base.@kwdef mutable struct ParamsAFM
    J = 10u"meV"
    K = 1u"meV"
    a = 1u"nm"
    g = 0.1^3 * 8a * Unitful.q
end
A_afm(k, p) = p.J + p.K
B_afm(k, p) = p.J * (cos(p.a * k[1]) + cos(p.a * k[2])) / 2
ϕ_afm(k, p) = atanh(-B_afm(k, p) / A_afm(k, p)) / 2
u_afm(k, p) = cosh(ϕ_afm(k, p))
v_afm(k, p) = sinh(ϕ_afm(k, p))

Ek_afm(k, p) = sqrt(A_afm(k, p)^2 - B_afm(k, p)^2)
ωk_afm(k, p) = Ek_afm(k, p) / ħ
##



p = ParamsAFM()
pp = ParamsPlasmons(Ef=50u"meV", ϵ=10)
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
Eunit = u"meV"
Nk = 100
krange = range(-pi / p.a, stop=pi / p.a, length=Nk)
to_plot = [Ek_afm([k, 0k], p) |> Eunit |> ustrip for k in krange]
ax.plot(krange ./ (pi / p.a) .|> upreferred .|> ustrip, to_plot, color="white")
krange_plasmons = range(-1u"nm^-1", 1u"nm^-1", length=1000)
ax.plot(krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, 0k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--")

lim_min = 1e-3u"meV"
lim_max = 2.5Ek_afm([pi / p.a, 0pi / p.a], p) |> Eunit

NE = 30
Nk_dos = 30
krange_dos = range(-pi / p.a, stop=pi / p.a, length=Nk_dos)
Erange = range(lim_min, stop=lim_max, length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos)
        res_loc = dos2_adaptive([k, 0k], E, x -> Ek_afm(x, p), p.a)
        res[i, j] = res_loc[1]
    end
end

cb = ax.pcolormesh(krange_dos ./ (pi / p.a) .|> upreferred .|> ustrip, Erange .|> u"meV" .|> ustrip, res .|> dosunit |> ustrip)
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"$\mathrm{DOS}_{2}$ ($ \mathrm{meV}^{-1}$)")
ax.set_xlabel(L"k")
ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
ax.set_ylabel(L"$E$ (meV)")
ax.set_xlim(upreferred.(extrema(krange ./ (pi / p.a)))...)
ax.set_xticks([-1, 0, 1], [L"-\pi/a", L"0", L"\pi/a"])

wsave(joinpath(save_dir, "dispersion_afm.pdf"), fig)
fig |> display


##

πq(k, q, p) = [0p.g, 0p.g, 2p.g]

vertex(k, q, p) = u_afm(k, p) * u_afm(q .- k, p) + v_afm(-k, p) * v_afm(q .+ k, p)



Γ_afm(q, E, p, pp; Nk=100) = Γ_square(q, E, k -> Ek_afm(k, p), (k, q) -> vertex(k, q, p), pp, (k, q) -> πq(k, q, p), p.a; Nk=Nk)


Γ_afm(q, E, p, pp) = Γ_adaptive(q, E, k -> Ek_afm(k, p), (k, q) -> vertex(k, q, p), pp, (k, q) -> πq(k, q, p), p.a)



Γ_afm_no_enhancement(q, E, p, pp) = Γ_adaptive(q, E, k -> Ek_afm(k, p), (k, q) -> 1, pp, (k, q) -> πq(k, q, p), p.a)





# function Γ(q::AbstractVector, p, pp; Nk=100)
#     Γ(q, Eq_plasmons(q, pp), p, pp; Nk=Nk)
# end

# function Γ(q::AbstractVector, E, p, pp; Nk=100)
#     krange = range(-pi / p.a, stop=pi / p.a, length=Nk)
#     sample_grid = Iterators.product(krange, krange)
#     res = 0.0im * u"meV^-1"
#     for k in sample_grid
#         res += (u_afm(k, p) * u_afm(q .- k, p) + v_afm(k, p) * v_afm(q .- k, p))^2 / (E - Ek_afm(k, p) - Ek_afm(q .- k, p) + 1e-2u"meV" * 1im)
#     end
#     return (Efield(q, pp)[3] * πq(q, p))^2 * res ./ Nk^2
# end



Nq_plasmons = 50
qrange_plasmons = range(1e-5u"nm^-1", stop=0.1u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
res = @showprogress [Γ_afm([q, 0.0u"nm^-1"], Eq_plasmons([q, 0.0u"nm^-1"], pp), p, pp)[1] |> u"meV" for q in qrange_plasmons]

res_no_enhancement = @showprogress [Γ_afm_no_enhancement([q, 0.0u"nm^-1"], Eq_plasmons([q, 0.0u"nm^-1"], pp), p, pp)[1] |> u"meV" for q in qrange_plasmons]

Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]
##
kunit_display = u"μm^-1"
Γ_display = u"meV"
fig, axs = plt.subplots(2, 1, height_ratios=[1, 2], sharex=true)
ax = axs[2]
# ax.axvline(qrange_plasmons[argmin(abs.(2Ek_afm([0kunit_display, 0kunit_display], p) .- Eplasmons_res))] |> kunit_display .|> ustrip, color="black", linestyle="--")
# ax.axvline(qrange_plasmons[argmin(abs.(2Ek_afm([pi / p.a, 0kunit_display], p) .- Eplasmons_res))] |> kunit_display .|> ustrip, color="black", linestyle="--")

color = "C0"

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, -res .|> imag .|> ustrip, color=color)

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, res .|> real .|> ustrip, color=color, linestyle="--")

ax.set_ylabel(L"$\Gamma$ (meV)", color=color)
ax.tick_params(axis="y", labelcolor=color)

ax2 = ax.twinx()
color = "C2"
ax2.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, "-", color=color)
ax2.set_ylabel(L"$\hbar\omega_q$ (meV)", color=color)

ax2.tick_params(axis="y", labelcolor=color)





# ax.axvline(2p.K |> Γ_display |> ustrip, color="black", linestyle="--")

ax.set_xlabel(L"$k$ ($\mathrm{\mu m}^{-1}$)")

ax = axs[1]
ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, ((-res .|> imag) ./ Eplasmons_res .|> ustrip), color="black")
# ax.set_yscale("log")
# ax.set_ylim(1e-5, res ./ Eplasmons_res .|> ustrip |> maximum)
# ax.axhline(γdefault)
ax.set_ylabel(L"\gamma")
fig.subplots_adjust(hspace=0)
ax.margins(x=0)
# wsave(joinpath(save_dir, "dispersioncut.pdf"), fig)
fig

## now do a spectral plot


Nq_plasmons = 50
qrange_plasmons = range(1e-5u"nm^-1", stop=0.1u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
NE = 50
Erange = range(1u"meV", stop=40u"meV", length=NE)
spectralunit = u"meV^-1"
res = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_Γ = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_afm([q, 0.0u"nm^-1"], E, p, pp)[1] |> u"meV"
        res_Γ[i, j] = Γloc
    end
end
##
Nz = 1
for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Nz * res_Γ[i, j]
        res[i, j] = spectral_plasmons(E, q, pp, Γloc)
    end
end

fig, ax = plt.subplots()
toplot = abs.(res) |> ustrip
fnorm = plt.matplotlib.colors.LogNorm(1e-4, maximum(toplot))
# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)

cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="nearest")
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

top_E = 2Ek_afm([pi / p.a, 0kunit_display], p)
top_q = qrange_plasmons[argmin(abs.(top_E .- Eplasmons_res))]
min_E = 2Ek_afm([0pi / p.a, 0kunit_display], p)
min_q = qrange_plasmons[argmin(abs.(min_E .- Eplasmons_res))]
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)
ax.axhline(min_E |> u"meV" |> ustrip, color="white", lw=0.5)

# ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, color="white")
ax.set_xlabel(L"$q$ ($\mathrm{\mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "plasmon_spectral_Nz_$(Nz).pdf"), fig)
fig


## zoom into top of the band
Nq_zoom = 50
qrange_zoom = range(top_q * 0.8, stop=top_q * 1.2, length=Nq_zoom)
@info Eq_plasmons(qrange_zoom[1], pp) .|> u"meV", Eq_plasmons(qrange_zoom[end], pp) .|> u"meV"
NE = 50
Erange_zoom = range(Eq_plasmons(qrange_zoom[1], pp), stop=Eq_plasmons(qrange_zoom[end], pp), length=NE)

res_zoom = zeros(0.0spectralunit |> typeof, NE, Nq_zoom)
res_zoom_undamped = zeros(0.0spectralunit |> typeof, NE, Nq_zoom)
res_zoom_Γ = zeros(0.0u"meV" |> typeof, NE, Nq_zoom)
@showprogress for (i, E) in enumerate(Erange_zoom)
    for (j, q) in enumerate(qrange_zoom)
        Γloc = Nz * Γ([q, 0.0u"nm^-1"], E, p, pp) |> u"meV"
        res_zoom_Γ[i, j] = Γloc
    end
end
##

for (i, E) in enumerate(Erange_zoom)
    for (j, q) in enumerate(qrange_zoom)
        Γloc = res_zoom_Γ[i, j]
        res_zoom[i, j] = L(E, Eq_plasmons(q, pp), Γloc + γdefault * Eq_plasmons(q, pp))
        res_zoom_undamped[i, j] = L(E, Eq_plasmons(q, pp), γdefault * Eq_plasmons(q, pp))

    end
end
fig, ax = plt.subplots()
toplot = abs.(res_zoom) |> ustrip
fnorm = plt.matplotlib.colors.LogNorm(1e-4, maximum(toplot))
# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)

cb = ax.pcolormesh(qrange_zoom .|> kunit_display .|> ustrip, Erange_zoom .|> u"meV" .|> ustrip, toplot, norm=fnorm)
fig.colorbar(cb, ax=ax)
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_zoom]

ax.plot(qrange_zoom .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, color="white")
fig
