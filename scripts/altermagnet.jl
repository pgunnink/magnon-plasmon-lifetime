using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
include("style.jl")
includet("EP_functions.jl")
Base.@kwdef mutable struct ParamsAM
    J = 10u"meV"
    δJ = 4u"meV"
    K = 1u"meV"
    a = 1u"nm"
    g = 0.1^3 * a * Unitful.q
end
A_am(k, p) = 4p.J + p.K
δA_am(k, p) = 2p.δJ * (cos(k[1] * p.a) - cos(k[2] * p.a))
B_am(k, p) = 2p.J * (cos(0.5p.a * k[1] - 0.5p.a * k[2]) + cos(0.5p.a * k[1] + 0.5p.a * k[2]))
ϕ_am(k, p) = atanh(-B_am(k, p) / A_am(k, p)) / 2
u_am(k, p) = cosh(ϕ_am(k, p))
v_am(k, p) = sinh(ϕ_am(k, p))

Ek_am(k, p) = sqrt(A_am(k, p)^2 - B_am(k, p)^2) * [1, 1] .+ [+δA_am(k, p), -δA_am(k, p)]
ωk_am(k, p) = Ek_am(k, p) / ħ

##

p = ParamsAM()
pp = ParamsPlasmons(Ef=10u"meV")



##
lim_min = 1e-3u"meV"
lim_max = 1.5Ek_am([pi / p.a, 0pi / p.a], p)[1] + 1.5Ek_am([pi / p.a, 0pi / p.a], p)[2]

NE = 100
Nk_dos = 10
end_point = 0.1pi / p.a

krange_dos_right = range(0 / p.a, stop=end_point, length=Nk_dos / 2 |> Int)
krange_dos_left = range(-end_point, stop=0 / p.a, length=Nk_dos / 2 |> Int)

Erange = range(lim_min, stop=lim_max, length=NE)
dosunit = u"meV^-1"
res = zeros(typeof(1.0dosunit), NE, Nk_dos)

@showprogress for (i, E) in enumerate(Erange)
    for (j, k) in enumerate(krange_dos_right)
        res_loc = dos2_adaptive([k, 0k], E, x -> Ek_am(x, p), p.a; rtol=1e-2)
        res[i, j+Nk_dos/2|>Int] = res_loc[1]
    end

    for (j, k) in enumerate(krange_dos_left)
        res_loc = dos2_adaptive([k, k], E, x -> Ek_am(x, p), p.a; rtol=1e-2)
        res[i, j] = res_loc[1]
    end
end


##
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
Eunit = u"meV"
Nk = 100
krange = range(0 / p.a, stop=end_point, length=Nk)
for i in 1:2
    to_plot = [Ek_am([k, 0k], p)[i] .|> Eunit .|> ustrip for k in krange]
    ax.plot(krange ./ (end_point) .|> upreferred .|> ustrip, to_plot, color="C$(i-1)")
end

krange = range(0 / p.a, stop=-end_point, length=Nk)
for i in 1:2
    to_plot = [Ek_am([k, k], p)[i] .|> Eunit .|> ustrip for k in krange]
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

wsave(joinpath(save_dir, "dispersion_am.pdf"), fig)
fig |> display
