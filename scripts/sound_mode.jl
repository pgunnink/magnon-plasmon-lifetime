includet("functions.jl")


##
q = 0.0002k_unit
Nω = 10
ωrange = range(0.03u"GHz", stop=30u"GHz", length=Nω)

RPA = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nω)
χ0_res = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nω)
RPA_iso = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nω)
χ0_res_iso = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nω)



@showprogress for j in 1:Nω
    ωj = ωrange[j]
    y = χ0_manual_integration([q, 0q], ωj, ParamsYIG()) .* χ_unit
    χ0_res[j] = y
    RPA[j] = y ./ (1 - y * Jz(q, ParamsYIG()))

    y = χ0_iso_manual_integration([q, 0q], ωj, ParamsYIG()) .* χ_unit
    χ0_res_iso[j] = y
    RPA_iso[j] = y ./ (1 - y * Jz(q, ParamsYIG()))
end

##
fig, ax = plt.subplots()
# ax.plot(ωrange .|> ustrip, real.(RPA .|> ustrip), label="full", color="C0")
# ax.plot(ωrange .|> ustrip, imag.(RPA .|> ustrip), "--", color="C0")
# ax.plot(ωrange .|> ustrip, imag.(RPA_iso .|> ustrip), "--", color="C1")
ax.plot(ωrange .|> ustrip, real.(RPA_iso .|> ustrip), label="full", color="C1")

# ax.axhline(0)
plt.legend()
fig
##

fig, axs = plt.subplots(1, 2)
axs[1].plot(ωrange .|> ustrip, real.(χ0_res .|> ustrip), label="full", color="C0")
axs[1].plot(ωrange .|> ustrip, real.(χ0_res_iso .|> ustrip), label="iso", color="C1")
axs[1].set_yscale("log")
axs[1].axhline(Jz(q, p)^-1 |> χ_unit |> ustrip)
plt.legend()

axs[2].plot(ωrange .|> ustrip, -imag.(χ0_res .|> ustrip), label="full", color="C0")
axs[2].plot(ωrange .|> ustrip, -imag.(χ0_res_iso .|> ustrip), label="iso", color="C1")
axs[2].set_yscale("log")
axs[1].set_title("real")
axs[2].set_title("imag")
# ax.plot(ωrange .|> ustrip, imag.(RPA_iso .|> ustrip), "--", color="C1")
# ax.plot(ωrange .|> ustrip, real.(RPA_iso .|> ustrip), label="full", color="C1")

# ax.axhline(0)
plt.legend()
fig


##
Nω = 10
ωrange = range(0.001u"GHz", stop=0.1u"GHz", length=Nω)
Nq = 10
qrange = range(1e-3u"μm^-1", stop=1e-2u"μm^-1 ", length=Nq)


res = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nq, Nω)


RPA = zeros(typeof((0.0 + 0.0im) * u"meV^-1"), Nq, Nω)

@showprogress for i in 1:Nq
    qi = qrange[i]
    for j in 1:Nω
        ωj = ωrange[j]
        res[i, j] = χ0([qi, 0qi], ωj, ParamsYIG()) * χ_unit
        RPA[i, j] = res[i, j] ./ (1 - res[i, j] * Jz(qi, ParamsYIG()))
    end
end

