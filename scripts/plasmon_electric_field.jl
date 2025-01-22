includet("plasmons.jl")
using PyPlot

pp = ParamsPlasmons()

fig, ax = plt.subplots()

N = 100
qrange = range(0.01u"nm^-1", stop=100u"nm^-1", length=N)

ax.plot(qrange .|> ustrip, [Efield(q, pp) ./ ωq_plasmons(q, pp) |> ustrip for q in qrange])

fig