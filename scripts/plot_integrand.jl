using PyPlot

includet("functions.jl")


## plot integrand
q = 0.00002k_unit
ω = 0.05u"GHz"
Λ = 1 / sqrt(p.Js / (k * p.T)) |> k_unit
Nk = 500
k_range = range(-0.1pi / p.a |> k_unit |> ustrip, stop=0.1pi / p.a |> k_unit |> ustrip, length=Nk)
dk = k_range[2] - k_range[1]
res = zeros(ComplexF64, Nk, Nk)
res_isotropic = zeros(ComplexF64, Nk, Nk)

for (i, kx) in k_range |> enumerate
    for (j, ky) in k_range |> enumerate
        res[i, j] = integrand_full(ω, [0q, q], p)([kx, ky])
        res_isotropic[i, j] = integrand_isotropic_full(ω, [0q, q], p)([kx, ky])
    end
end
@info Nk
@info sum(res_isotropic) * dk^2
@info sum(res) * dk^2
##
fnorm = plt.Normalize(extrema(res_isotropic .|> real)...)
fig, axs = plt.subplots(1, 2)
axs[1].pcolormesh(k_range, k_range, real.(res_isotropic), norm=fnorm)
axs[2].pcolormesh(k_range, k_range, real.(res), norm=fnorm)
for ax in axs
    ax.set_aspect("equal")
end
fig
