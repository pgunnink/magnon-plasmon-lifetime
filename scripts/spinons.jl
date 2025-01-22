using Meshes
using ProgressMeter
using PyPlot
using Unitful

includet("plasmons.jl")
include("style.jl")
includet("EP_functions.jl")
includet("spinon_definitions.jl")

kunit_display = u"μm^-1"
Γ_display = u"meV"

## form the BZ, defining a as the nn distance
p = ParamsSpinons()
b1 = [pi, -pi / sqrt(3)] ./ 2p.a
b2 = [0, 4pi / sqrt(3)] ./ 2p.a

BZ_corners = [
    b2 ./ 2, b1 .+ b2 ./ 2, -b2 ./ 2, -b1 .- b2 ./ 2
]
kunit = u"nm^-1"

BZ_square = Ngon([(x[1], x[2]) .|> kunit .|> ustrip for x in BZ_corners]...)
sample_grid = []
y1 = b1 ./ 2
y2 = b2 ./ 2 .+ b1 ./ 2

sample_grid = []
fig, ax = plt.subplots()
for i in -1:0.01:1
    for j in -1:0.01:1
        point = i .* y1 + j .* y2
        push!(sample_grid, point)
    end
end
ax.scatter([x[1] for x in sample_grid] .|> kunit .|> ustrip, [x[2] for x in sample_grid] .|> kunit .|> ustrip)

ax.scatter([x[1] for x in BZ_corners] .|> kunit .|> ustrip, [x[2] for x in BZ_corners] .|> kunit .|> ustrip, color="red")
plt.gca().set_aspect("equal")

Gamma = [0 / p.a, 0 / p.a]
M1 = BZ_corners[1]
K1 = BZ_corners[2]

ax.scatter((Gamma .|> kunit .|> ustrip)..., color="C5")
ax.scatter((M1 .|> kunit .|> ustrip)..., color="C6")
ax.scatter((K1 .|> kunit .|> ustrip)..., color="C7")


Δk = 0.1 / p.a
path = [Gamma, M1, BZ_corners[4], Gamma]
kpoints = []
pos_points = []

for n in 1:length(path)-1
    push!(kpoints, path[n])
    push!(pos_points, length(kpoints))
    p2 = path[n+1]
    p1 = path[n]
    δk = p2 .- p1
    Nk = norm(δk) / Δk
    for i in 1:Nk
        push!(kpoints, p1 + δk / Nk * (i))
    end
end
push!(pos_points, length(kpoints))
@info kpoints[pos_points]
ax.scatter([x[1] for x in kpoints] .|> kunit .|> ustrip, [x[2] for x in kpoints] .|> kunit .|> ustrip, color="black", s=1)


fig
## CALCULATE EIGENFUNCTIONS AT THE K POINTS 

dirac_point = b2 / 4
@info round.(Hi(dirac_point, p) .|> u"meV" .|> ustrip; digits=2)
F = eigen(Hi(dirac_point, p) .|> u"meV" .|> ustrip)
@info F.values
@info L(0, F.values[4] - F.values[3], 1e-2)
@info (F.vectors[:, 3] .|> abs) .^ 2
@info (F.vectors[:, 4] .|> abs) .^ 2
EP_diag = ((F.vectors |> adjoint) * Pij(dirac_point, [1u"V/m", 1u"V/m", 1u"V/m"], p) * F.vectors) .|> u"meV"
@info EP_diag[3, 4] .|> ustrip
##
Ni = 1000
irange = range(-1, stop=1, length=Ni)
res = zeros(Float64, Ni, Ni)
for (i, xi) in enumerate(irange)
    for (j, yi) in enumerate(irange)
        k = xi .* y1 + yi .* y2
        Ek = eigvals(Hi(k, p) .|> u"meV" .|> ustrip)
        res[i, j] = L(0, Ek[4] - Ek[3], 1e-2)
    end
end
@info sum(res) / Ni^2
fig, ax = plt.subplots()
ax.pcolormesh(irange, irange, res)
fig

##
function DOS2_spinons(q, E, p; rtol=1e-3, η=1e-2u"meV", maxevals=Int(1e4), atol=0)
    Estripped = E |> u"meV" |> ustrip
    ηstripped = η |> u"meV" |> ustrip
    g = x -> begin
        k = x[1] .* y1 + x[2] .* y2
        Ek = eigvals(Hi(k, p) .|> u"meV" .|> ustrip)
        Ekq = eigvals(Hi(k .+ q, p) .|> u"meV" .|> ustrip)
        res = 0
        for i in 1:3
            for j in 4:6
                res += L(Estripped, Ek[j] - Ekq[i], ηstripped)
            end
        end
        res
    end
    (1 / 2^2) .* hcubature(g, .-ones(2), ones(2), initdiv=10, rtol=rtol, maxevals=maxevals, atol=atol) .* u"meV^-1"
end


##
res_kpath = zeros(typeof(1.0u"meV"), length(kpoints), 6)
NE = 50
Erange = range(0u"meV", stop=6u"meV", length=NE)
res_dos = zeros(typeof(1.0u"meV^-1"), NE, length(kpoints))
@showprogress for (i, k) in enumerate(kpoints)
    Hloc = Hi(k, p) .|> u"meV" .|> ustrip
    res_kpath[i, :] = eigvals(Hloc) .* u"meV"
    for (j, E) in Erange |> enumerate
        res_loc = DOS2_spinons(k, E, p)
        res_dos[j, i] = res_loc[1]
        # @info res_loc[2] / res_loc[1] |> ustrip
    end
end
##
NE_cut = 100
Erange_cut = range(0u"meV", stop=6u"meV", length=NE_cut)
res_cut = zeros(typeof(1.0u"meV^-1"), NE_cut)

@showprogress for (i, E) in enumerate(Erange_cut)
    res_cut[i] = DOS2_spinons([0 / p.a, 0 / p.a], E, p;)[1]
end

##
fig, axs = plt.subplots(2, 2, figsize=(fig_width, fig_height), width_ratios=[3, 1], height_ratios=[1, 20])
ax = axs[2, 1]
bottom_band = eigvals(Hi([0 / p.a, 0 / p.a], p) .|> u"meV" .|> ustrip)[1]
ax.plot((res_kpath .|> ustrip) .- bottom_band, color="white", alpha=1)
ax.margins(x=0)
ax.axhline(eigvals(Hi(dirac_point, p) .|> u"meV" .|> ustrip)[3] - bottom_band, linestyle="dashed", color="white")

toplot = res_dos .|> u"meV^-1" |> ustrip
@info extrema(res_dos)
fnorm = plt.matplotlib.colors.Normalize(minimum(toplot), maximum(toplot))
# fnorm = plt.matplotlib.colors.LogNorm(minimum(toplot), maximum(toplot))


cb = ax.pcolormesh(1:length(kpoints), Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="gouraud")
ax.set_xticks(pos_points; labels=[L"(0,0)", L"(0, \frac{\pi}{ 2\sqrt{3}a})", L"(-\frac{\pi}{2a},\frac{\sqrt{3}\pi}{2a})", L"(0,0)"])
ax.set_ylabel(L"$\epsilon$ (meV)")

clb = fig.colorbar(cb, cax=axs[1, 1], orientation="horizontal", location="top")
clb.ax.set_title(L"$\mathrm{DOS}_{2}$ ($ \mathrm{meV}^{-1}$)")


ax = axs[2, 2]
ax.plot(res_cut .|> ustrip, Erange_cut .|> u"meV" .|> ustrip)
# ax.set_xscale("log")
ax.set_yticklabels([])
ax.set_yticks([])
ax.set_xlabel(L"$\mathrm{DOS}_2$ $(\mathrm{meV}^{-1}$)")
ax.margins(y=0)
axs[1, 2].axis("off")
fig.subplots_adjust(wspace=0.05, hspace=0.05)


wsave(joinpath(save_dir, "dispersion_spinons.pdf"), fig)


fig |> display

##


##
NE = 10
Erange = range(0.001u"meV", stop=8u"meV", length=NE)
res = @showprogress [DOS2_spinons([0.00 / p.a, 0.0 / p.a], E, p; rtol=1e-4, η=1e-1u"meV") for E in Erange]
fig, ax = plt.subplots()
ax.plot(Erange .|> ustrip, [x[1] for x in res] .|> ustrip)
ax.set_xlabel(L"\hbar\omega")
ax.set_ylabel(L"\mathrm{DOS}_2")
fig

##

function Γ_spinons(q, E, p, pp; η=1e-1u"meV", rtol=1e-3, maxevals=Int(1e4))
    Estripped = E |> u"meV" |> ustrip
    ηstripped = η |> u"meV" |> ustrip
    g = x -> begin
        k = x[1] .* y1 + x[2] .* y2
        Ek, Uk = eigen(Hi(k, p) .|> u"meV" .|> ustrip)
        Ekplusq, Ukplusq = eigen(Hi(k .+ q, p) .|> u"meV" .|> ustrip)
        EP = Pij(k .+ q, Efield(q, E / ħ, pp), p) .|> u"meV" .|> ustrip
        diagEP2 = abs.(adjoint(Ukplusq) * EP * Uk) .^ 2
        res = 0im
        for i in 1:3
            for j in 4:6
                res += 2diagEP2[i, j] / (Estripped + Ek[i] - Ekplusq[j] + 1im * ηstripped)
            end
        end
        res
    end
    (1 / 2^2) .* hcubature(g, .-ones(2), ones(2), rtol=rtol, maxevals=maxevals, norm=imag) .* u"meV"
end


pp = ParamsPlasmons(Ef=50u"meV", ϵ=30)
qtest = [0.001u"nm^-1", 0.001u"nm^-1"]
@time Γ_spinons(qtest, Eq_plasmons(qtest, pp), p, pp; η=1e-1u"meV", maxevals=Int(1e4), rtol=1e-6)

##
@info Eq_plasmons(qtest, pp) |> u"meV"


@info Γ_spinons(qtest, Eq_plasmons(qtest, pp), p, pp; η=1e-1u"meV", maxevals=Int(1e4), rtol=1e-6)
@info DOS2_spinons(qtest, Eq_plasmons(qtest, pp), p; η=1e-1u"meV", maxevals=Int(1e4), atol=0, rtol=1e-6)
##
Nq_plasmons = 50
qrange_plasmons = range(1e-5u"nm^-1", stop=0.01u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
res_full = @showprogress [Γ_spinons([0q, q], Eq_plasmons([q, 0.0u"nm^-1"], pp), p, pp, η=0.1u"meV") .|> u"meV" for q in qrange_plasmons]
res = [x[1] for x in res_full]
@info [x[2] / x[1] |> ustrip for x in res_full]
@info sum([x[2] / x[1] |> ustrip for x in res_full]) / length(res_full)
Eplasmons_res = [Eq_plasmons([0.0u"nm^-1", q], pp) .|> u"meV" for q in qrange_plasmons]

##

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


##

kunit_display = u"μm^-1"
Γ_display = u"meV"

Nq_plasmons = 100
qrange_plasmons = range(1e-5u"nm^-1", stop=0.012u"nm^-1", length=Nq_plasmons)
@info Eq_plasmons(qrange_plasmons[1], pp) .|> u"meV", Eq_plasmons(qrange_plasmons[end], pp) .|> u"meV"
NE = 100
Erange = range(0.001u"meV", stop=9u"meV", length=NE)
spectralunit = u"meV^-1"
res = zeros(0.0spectralunit |> typeof, NE, Nq_plasmons)
res_Γ = zeros(0.0im * u"meV" |> typeof, NE, Nq_plasmons)
@showprogress for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = Γ_spinons([0q, q], E, p, pp; η=0.01u"meV", maxevals=Int(1e4))[1] |> u"meV"
        res_Γ[i, j] = Γloc
    end
end
##
Nz = 1e7
@info Nz * p.a |> u"cm"
# correction = q -> sum(exp(-2κ(q, pp) * n * lattice_c) for n in 1:Nz)
correction = q -> (1 - exp(-2Nz * κ(q, pp) * p.a)) / (-1 + exp(2κ(q, pp) * p.a))
for (i, E) in enumerate(Erange)
    for (j, q) in enumerate(qrange_plasmons)
        Γloc = correction(q) * res_Γ[i, j]
        res[i, j] = spectral_plasmons(E, q, pp, Γloc)
    end
end

fig, ax = plt.subplots()
toplot = abs.(res) |> ustrip
fnorm = plt.matplotlib.colors.LogNorm(1e-3, maximum(toplot))
# fnorm = plt.matplotlib.colors.Normalize(extrema(toplot)...)

cb = ax.pcolormesh(qrange_plasmons .|> kunit_display .|> ustrip, Erange .|> u"meV" .|> ustrip, toplot, norm=fnorm, shading="nearest")
clb = fig.colorbar(cb, ax=ax)
clb.ax.set_title(L"A_p")
Eplasmons_res = [Eq_plasmons([q, 0.0u"nm^-1"], pp) .|> u"meV" for q in qrange_plasmons]

top_E = 6u"meV"
ax.axhline(top_E |> u"meV" |> ustrip, color="white", lw=0.5)

ax.plot(qrange_plasmons .|> kunit_display .|> ustrip, Eplasmons_res .|> Γ_display .|> ustrip, linestyle="dotted", color="white")
ax.set_ylim(extrema(Erange .|> u"meV" .|> ustrip)...)
ax.set_xlabel(L"$q$ ($\mathrm{\mu m}^{-1}$)")
ax.set_ylabel(L"$\hbar\omega$ (meV)")
wsave(joinpath(save_dir, "spinon_spectral_Nz_$(Nz).pdf"), fig)
fig

