using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
Base.@kwdef mutable struct ParamsMnF2D
    J = 1u"meV"
    δJ = 0.2u"meV"
    K = 1u"meV"
    a = 5u"Å"
    c = 3u"Å"
    g = 8 * 0.1^3 * a * Unitful.q * cos(deg2rad(128 / 2))
    α = 1e-3
end

kunit_display = u"μm^-1"
Γ_display = u"meV"

A_MnF2D(k, p) = 4p.J + p.K
B_MnF2D(k, p) = 2p.J * (cos(0.5p.a * k[1] - 0.5p.a * k[2]) + cos(0.5p.a * k[1] + 0.5p.a * k[2]))
δA_MnF2D(k, p) = 2p.δJ * (cos(k[1] * p.a) - cos(k[2] * p.a))

ϕ_MnF2D(k, p) = atanh(-B_MnF2D(k, p) / A_MnF2D(k, p)) / 2
u_MnF2D(k, p) = cosh(ϕ_MnF2D(k, p))
v_MnF2D(k, p) = sinh(ϕ_MnF2D(k, p))

Ek_MnF2D(k, p) = sqrt(A_MnF2D(k, p)^2 - B_MnF2D(k, p)^2) * [1, 1] .+ δA_MnF2D(k, p) .* [1, -1]
ωk_MnF2D(k, p) = Ek_MnF2D(k, p) / ħ




p = ParamsMnF2D()
pp = ParamsPlasmons(Ef=10u"meV", ϵ=30)


##
lim_min = 1e-3u"meV"
lim_max = 1.5Ek_MnF2D([pi / p.a, 0pi / p.a], p)[1] + 1.5Ek_MnF2D([pi / p.a, 0pi / p.a], p)[2]


##

NE = 20
Nk_dos = 10
end_point = pi / p.a

krange_dos_right = range(0 / p.a, stop=end_point, length=Nk_dos / 2 |> Int)
krange_dos_left = range(-end_point, stop=0 / p.a, length=Nk_dos / 2 |> Int)

Erange = range(lim_min, stop=lim_max, length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos_right)
        res_loc = dos2_adaptive([k, 0k], E, x -> Ek_MnF2D(x, p), p.a; rtol=1e-3, maxevals=Int(1e4), N=2)
        res[i, j+Nk_dos/2|>Int] = res_loc[1]
    end

    for (j, k) in enumerate(krange_dos_left)
        res_loc = dos2_adaptive([0k, k], E, x -> Ek_MnF2D(x, p), p.a; rtol=1e-3, maxevals=Int(1e4), N=2)
        res[i, j] = res_loc[1]
    end
end

##

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
Eunit = u"meV"
Nk = 100
krange = range(0 / p.a, stop=end_point, length=Nk)
for i in 1:2
    to_plot = [Ek_MnF2D([k, 0k], p)[i] .|> Eunit .|> ustrip for k in krange]
    ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color="C$(i-1)")
end

krange = range(0 / p.a, stop=-end_point, length=Nk)
for i in 1:2
    to_plot = [Ek_MnF2D([0k, k], p)[i] .|> Eunit .|> ustrip for k in krange]
    ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color="C$(i-1)")
end
toplot = res .|> dosunit |> ustrip
@info maximum(toplot)
fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
fnorm = plt.matplotlib.colors.Normalize(1e-2, maximum(toplot))

cb = ax.pcolormesh([krange_dos_left; krange_dos_right] ./ (end_point) .|> upreferred .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm)
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"$\mathrm{DOS}_{2}$ ($ \mathrm{meV}^{-1}$)")
ax.set_xlabel(L"k")
ax.set_ylim(0.0, lim_max |> Eunit |> ustrip)
ax.set_ylabel(L"$E$ (meV)")
ax.set_xlim(upreferred.(extrema(krange ./ (pi / p.a)))...)
ax.set_xticks([-1, 0, 1], [L"-\pi/a", L"0", L"\pi/a"])
krange_plasmons = range(-1u"nm^-1", 1u"nm^-1", length=1000)
ax.plot(krange_plasmons ./ (pi / p.a) .|> upreferred .|> ustrip, [Eq_plasmons([k, 0k], pp) |> Eunit |> ustrip for k in krange_plasmons], color="white", linestyle="--")

wsave(joinpath(save_dir, "dispersion_MnF2D.pdf"), fig)
fig |> display


##

nn = [[0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [-0.5, -0.5, 0.5], [-0.5, 0.5,
    0.5]]
P = [[-0.212466, -0.212466,
        0.95379], [0.212466, -0.212466, -0.95379], [0.212466, 0.212466,
        0.95379], [-0.212466, 0.212466, -0.95379]]
PAB(k, p) = sum([x .* exp(1im * (y .* [p.a, p.a, p.c]) ⋅ [k; 0u"m^-1"]) for (x, y) in zip(P, nn)])

@info PAB([0.0 / p.a, 0.01 / p.a], p)
@info PAB([0.01 / p.a, 0.0 / p.a], p)


vertex(k, q, p) = 1
πq(k, q, p) = [p.g .* (PAB(k, p) .* u_MnF2D(k, p)^2 .+ conj.(PAB(k, p)) .* v_MnF2D(k, p)^2), p.g .* (conj.(PAB(k, p)) .* u_MnF2D(k, p)^2 .+ PAB(k, p) .* v_MnF2D(k, p)^2)]

πq(k, q, p) = [p.g .* [0, -0.849865, 0.0] .* (u_MnF2D(k, p)^2 + v_MnF2D(k, p)^2), p.g .* [0, -0.849865, 0.0] .* (u_MnF2D(k, p)^2 + v_MnF2D(k, p)^2)]


vertexbare(k, p) = p.g .* [
    -0.84 * cos(0.5k[1] * p.a) * cos(0.5k[2] * p.a),
    -3.81im * sin(0.5k[1] * p.a) * sin(0.5k[2] * p.a),
    -0.84im * sin(0.5k[1] * p.a) * cos(0.5k[2] * p.a)
]


πq(k, q, p) = [vertexbare(k .- q, p) .* v_MnF2D(k, p) * v_MnF2D(k .+ q, p), vertexbare(q .- k, p) .* u_MnF2D(k, p) * u_MnF2D(k .+ q, p)]


Γ_MnF2D(q, E, p, pp) = Γ_adaptive(q, E, k -> Ek_MnF2D(k, p), (k, q) -> 1, pp, (k, q) -> πq(k, q, p), p.a; maxevals=Int(1e4), η=0.0001u"meV", N=2)


##
Nq_plasmons = 50
qrange_plasmons = range(1e-5u"nm^-1", stop=0.3u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
NE = 50
Erange = range(2u"meV", stop=15u"meV", length=NE)
spectralunit = u"meV^-1"
res = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_Γ = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_MnF2D([q, 0.0u"nm^-1"], E, p, pp)[1] |> u"meV"
        res_Γ[i, j] = Γloc
    end
end
##
Nz = 100
for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Nz * res_Γ[i, j]
        res[i, j] = spectral_plasmons(E, q, pp, Γloc)
    end
end

fig, ax = plt.subplots()
toplot = abs.(res) |> ustrip
fnorm = plt.matplotlib.colors.LogNorm(1e-2, maximum(toplot))
# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)

cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="gouraud")
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

top_E = Ek_MnF2D([pi / p.a, 0pi / p.a, pi / p.a], p)[1] + Ek_MnF2D([pi / p.a, 0pi / p.a, pi / p.a], p)[2]
# top_q = qrange_plasmons[argmin(abs.(top_E .- Eplasmons_res))]
min_E = 2Ek_MnF2D([0pi / p.a, 0pi / p.a, 0pi / p.a], p)[1]
# min_q = qrange_plasmons[argmin(abs.(min_E .- Eplasmons_res))]
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)
ax.axhline(min_E |> u"meV" |> ustrip, color="white", lw=0.5)

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, "--", color="white")
ax.set_ylim(extrema(Erange .|> u"meV" .|> ustrip)...)
ax.set_xlabel(L"$q$ ($\mathrm{\mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "MnF2D_2D_spectral_Nz_$(Nz).pdf"), fig)
fig

