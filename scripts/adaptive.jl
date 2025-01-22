using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
using HCubature

includet("plasmons.jl")
include("style.jl")
includet("EP_functions.jl")
Base.@kwdef mutable struct ParamsAFM
    J = 10u"meV"
    K = 1u"meV"
    a = 1u"nm"
    g = 0.1^3 * a * Unitful.q
end
A_afm(k, p) = p.J + p.K
B_afm(k, p) = p.J * (cos(p.a * k[1]) + cos(p.a * k[2])) / 2
ϕ_afm(k, p) = atanh(-B_afm(k, p) / A_afm(k, p)) / 2
u_afm(k, p) = cosh(ϕ_afm(k, p))
v_afm(k, p) = sinh(ϕ_afm(k, p))

Ek_afm(k, p) = sqrt(A_afm(k, p)^2 - B_afm(k, p)^2)
ωk_afm(k, p) = Ek_afm(k, p) / ħ


q = [0.01 / p.a, 0.04 / p.a]
E = 30u"meV"
ηsmear = 1e-1u"meV"

L(x, y, η=ηsmear) = 1 / π * η / ((x - y)^2 + (η)^2)

function dos2_adaptive(q, E, f, a)
    kunit = u"nm^-1"
    resunit = u"meV^-1"
    g = k -> ustrip(resunit,
        L(E, f(k .* kunit) + f(q .- k .* kunit)) / (2pi)^2
    )
    lim = ustrip(kunit, pi / a)
    hcubature(g, [-lim, -lim], [lim, lim])
end


##
@time dos2_square(q, E, x -> Ek_afm(x, p), p.a; Nk=1000)
@time dos2_adaptive(q, E, x -> Ek_afm(x, p), p.a;)


##

πq(q, p) = [0p.g, 0p.g, 2p.g * (cos(q[1] * p.a) + cos(q[2] * p.a))]

vertex(k, q, p) = u_afm(k, p) * u_afm(q .- k, p) + v_afm(.-k, p) * v_afm(q .+ k, p)
##

@info Γ_square(q, E, k -> Ek_afm(k, p), (k, q) -> vertex(k, q, p), pp, q -> πq(q, p), p.a; Nk=500) |> u"meV"
@info Γ_adaptive(q, E, k -> Ek_afm(k, p), (k, q) -> vertex(k, q, p), pp, q -> πq(q, p), p.a)