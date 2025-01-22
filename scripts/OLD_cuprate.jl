using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
Base.@kwdef mutable struct ParamsCuprate
    J = 105u"meV"
    K = 3u"meV"
    a = 2.5u"Å"
    g = 0.1^3 * a * Unitful.q
    # g = 4.4e-4u"meV / (kV / cm)"
    Jp = -0.5 * 26.25u"meV"
end
A_cupr(k, p) = 4p.J + p.K + p.Jp * (4 - 2cos(2k[1] * p.a) - 2cos(2k[2] * p.a))
B_cupr(k, p) = 2p.J * (cos(p.a * k[1]) + cos(p.a * k[2]))
ϕ_cupr(k, p) = atanh(-B_cupr(k, p) / A_cupr(k, p)) / 2
u_cupr(k, p) = cosh(ϕ_cupr(k, p))
v_cupr(k, p) = sinh(ϕ_cupr(k, p))

Ek_cupr(k, p) = sqrt(A_cupr(k, p)^2 - B_cupr(k, p)^2)
ωk_cupr(k, p) = Ek_cupr(k, p) / ħ




p = ParamsCuprate()
pp = ParamsPlasmons(Ef=147u"meV", ϵ=10)
##
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
Eunit = u"meV"
Nk = 100
krange = range(-pi / p.a, stop=0pi / p.a, length=Nk)
to_plot = [Ek_cupr([0k, k], p) |> Eunit |> ustrip for k in krange]
ax.plot(krange ./ (pi / p.a) .|> upreferred .|> ustrip, to_plot, color="white")
krange = range(-0pi / p.a, stop=pi / p.a, length=Nk)

to_plot = [Ek_cupr([k, k], p) |> Eunit |> ustrip for k in krange]
ax.plot(krange ./ (pi / p.a) .|> upreferred .|> ustrip, to_plot, color="white")



krange_plasmons = range(-pi / p.a, pi / p.a, length=1000)
ax.plot(krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, 0k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--")



lim_min = 1e-3u"meV"
lim_max = 2.5Ek_cupr([pi / 2p.a, pi / 2p.a], p) |> Eunit

NE = 60
Nk_dos = 60
krange_dos = range(-pi / p.a, stop=pi / p.a, length=Nk_dos)
Erange = range(lim_min, stop=lim_max, length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos[1:Int(Nk_dos // 2)])
        res_loc = dos2_adaptive([k, 0k], E, x -> Ek_cupr(x, p), p.a; rtol=0.01, η=3u"meV", maxevals=Int(1e4))[1]
        # res_loc = dos2_square([k, 0k], E, x -> Ek_cupr(x, p), p.a; Nk=300, η=1u"meV")[1]
        res[i, j] = res_loc[1]
    end
end

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos[Int(Nk_dos // 2)+1:end])
        # res_loc = dos2_adaptive([k, 0k], E, x -> Ek_cupr(x, p), p.a; rtol=0.01, η=3u"meV")[1]
        # res_loc = dos2_square([k, k], E, x -> Ek_cupr(x, p), p.a; Nk=300, η=1u"meV")[1]
        res_loc = dos2_adaptive([k, k], E, x -> Ek_cupr(x, p), p.a; rtol=0.01, η=3u"meV", maxevals=Int(1e4))[1]

        res[i, j+Int(Nk_dos // 2)] = res_loc[1]
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

wsave(joinpath(save_dir, "dispersion_cuprate.pdf"), fig)
fig |> display

##
NE = 20
Nk_dos = 20
krange_dos = range(pi / p.a * 0.45, stop=pi / p.a * 0.55, length=Nk_dos)
Erange = range(450u"meV", stop=550u"meV", length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)
@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos)
        res_loc = dos2_adaptive([k, k], E, x -> Ek_cupr(x, p), p.a; η=1u"meV")[1]
        res[i, j] = res_loc[1]
    end
end

pp = ParamsPlasmons(Ef=105u"meV", ϵ=10)

fig, ax = plt.subplots()
cb = ax.pcolormesh(krange_dos ./ (pi / p.a) .|> upreferred .|> ustrip, Erange .|> u"meV" .|> ustrip, res .|> dosunit |> ustrip)
ax.plot(krange_dos ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, k], pp) |> Eunit |> ustrip for k in krange_dos], color="white", linestyle="--")

clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"$\mathrm{DOS}_{2}$ ($ \mathrm{meV}^{-1}$)")
ax.set_xlabel(L"k")
# ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
ax.set_ylabel(L"$E$ (meV)")
# ax.set_xlim(upreferred.(extrema(krange ./ (pi / p.a)))...)
# ax.set_xticks([-1, 0, 1], [L"-\pi/a", L"0", L"\pi/a"])

fig


##
πq(k, q, p) = 2p.g .* [sin(k[2] * p.a), sin(k[1] * p.a), 0]
vertex(k, q, p) = 1

Γ_cupr(q, E, p, pp; rtol=1e-3, Γη=1e-2u"meV", maxevals=Int(1e3)) = Γ_adaptive(q, E, k -> Ek_cupr(k, p), (k, q) -> vertex(k, q, p), pp, (k, q) -> πq(k, q, p), p.a; η=Γη, rtol=rtol, maxevals=maxevals)


Γ_cupr_square(q, E, p, pp; Nk=100, η=1e-3u"meV") = Γ_square(q, E, k -> Ek_cupr(k, p), (k, q) -> vertex(k, q, p), pp, (k, q) -> πq(k, q, p), p.a; Nk=Nk, η=η)


#Γ_square(q, E, k -> Ek_cupr(k, p), vertex, pp, q -> πq(q, p), p.a; Nk=Nk)
# Γ_afm(q, E, p, pp) =


##
Nq_plasmons = 100
qrange_plasmons = range(0.01u"nm^-1", stop=12u"nm^-1", length=Nq_plasmons)
qrange_plasmons = range(3u"nm^-1", stop=6u"nm^-1", length=Nq_plasmons)
qrange_plasmons = range(pi / p.a * 0.466, stop=pi / p.a * 0.48, length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
# res = @showprogress [Γ_cupr_([q, q], Eq_plasmons([q, q], pp), p, pp; rtol=0.01, Γη=1e-2u"meV")[1] |> u"meV" for q in qrange_plasmons]
res = @showprogress [Γ_cupr_square([q, q], Eq_plasmons([q, q], pp), p, pp; Nk=1000, η=1e-1u"meV")[1] |> u"meV" for q in qrange_plasmons]
Eplasmons_res = [Eq_plasmons([q, q], pp) .|> u"meV" for q in qrange_plasmons]
##

kunit_display = u"μm^-1"
Γ_display = u"meV"
fig, axs = plt.subplots(2, 1, height_ratios=[1, 2], sharex=true)
ax = axs[2]
Nz = 1
color = "C0"

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, -Nz * res .|> imag .|> ustrip, color=color)
ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Nz * res .|> real .|> ustrip, color=color, linestyle="--")

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
ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, ((-Nz * res .|> imag) ./ Eplasmons_res .|> ustrip), color="black")
# ax.set_yscale("log")
# ax.set_ylim(1e-5, res ./ Eplasmons_res .|> ustrip |> maximum)
# ax.axhline(γdefault)
ax.set_ylabel(L"\gamma")
fig.subplots_adjust(hspace=0)
ax.margins(x=0)
wsave(joinpath(save_dir, "dispersioncut_cuprate.pdf"), fig)
fig

##


Nq_plasmons = 100
qrange_plasmons = range(1e-5u"nm^-1", stop=0.3u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
NE = 100
Erange = range(1u"meV", stop=100u"meV", length=NE)
spectralunit = u"meV^-1"
res = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_Γ = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_cupr([q, 0.0u"nm^-1"], E, p, pp)[1] |> u"meV"
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

cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm)
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

top_E = 2Ek_cupr([pi / 2p.a, pi / 2p.a], p)
top_q = qrange_plasmons[argmin(abs.(top_E .- Eplasmons_res))]
min_E = 2Ek_cupr([0pi / p.a, 0kunit_display], p)
min_q = qrange_plasmons[argmin(abs.(min_E .- Eplasmons_res))]
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)
ax.axhline(min_E |> u"meV" |> ustrip, color="white", lw=0.5)

# ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, color="white")
ax.set_xlabel(L"$q$ ($\mathrm{\mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "plasmon_spectral_cuprate.pdf"), fig)
fig
