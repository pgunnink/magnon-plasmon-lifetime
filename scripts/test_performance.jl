using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
Base.@kwdef mutable struct ParamsCuprate
    J = 130u"meV"
    K = 0.001u"meV"
    a = 2.5u"Å"
    g = 0.1^3 * a * Unitful.q
    # g = 4.4e-4u"meV / (kV / cm)"
    Jp = 32.5u"meV"
    J0 = J + K + 4Jp
end



A_cupr(k, p) = p.J0 - p.Jp * (2cos(2k[1] * p.a) + 2cos(2k[2] * p.a))
B_cupr(k, p) = p.J * (cos(p.a * k[1]) + cos(p.a * k[2])) / 2

ϕ_cupr(k, p) = atanh(-B_cupr(k, p) / A_cupr(k, p)) / 2
u_cupr(k, p) = cosh(ϕ_cupr(k, p))
v_cupr(k, p) = sinh(ϕ_cupr(k, p))

Ek_cupr(k, p) = sqrt(A_cupr(k, p)^2 - B_cupr(k, p)^2)
ωk_cupr(k, p) = Ek_cupr(k, p) / ħ


##
p = ParamsCuprate()
qtest = 1u"nm^-1"
@time Ek_cupr([qtest, 0qtest], p)

##
@time L(Ek_cupr([qtest, 2qtest], p), Ek_cupr([qtest, 0qtest], p), 1e-4u"meV")
kunit = u"nm^-1"
resunit = u"meV^-1"
E = 1u"meV"
g = k -> ustrip(resunit,
    L(E, Ek_cupr(k .* kunit, p) + Ek_cupr(qtest .- k .* kunit, p)) / (2pi)^2
)


##

function f_perform(q, E, f; rtol=1e-3)
    kunit = u"nm^-1"
    resunit = u"meV^-1"
    k -> ustrip(resunit,
        L(E, f(k .* kunit) + f(q .- k .* kunit)) / (2pi)^2
    )
end
ftest = f_perform([qtest, qtest], 1u"meV", x -> Ek_cupr(x, p))
@time ftest([2.0, 3.0])