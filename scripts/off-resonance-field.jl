includet("plasmons.jl")

Nq_plasmons = 200
qrange_plasmons = range(1e-5u"nm^-1", stop=0.02u"nm^-1", length=Nq_plasmons)
@info ωq_plasmons(qrange_plasmons[1], pp) .|> u"THz", ωq_plasmons(qrange_plasmons[end], pp) .|> u"THz"

Nω = 100
ωrange = range(0.1u"THz", stop=40u"THz", length=Nω)

pp = ParamsPlasmons()
res = zeros(0.0u"V/m" |> typeof, Nω, Nq_plasmons)
for (i, ω) in enumerate(ωrange)
    for (j, q) in enumerate(qrange_plasmons)
        res[i, j] = E(q, ω * (1 + 1e-4 * 1im), pp) |> real
    end
end

fig, ax = plt.subplots()
cb = ax.pcolormesh(qrange_plasmons .|> ustrip, ωrange .|> u"THz" .|> ustrip, res |> ustrip .|> real)
fig