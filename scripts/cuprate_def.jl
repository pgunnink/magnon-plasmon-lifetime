using Unitful
using Printf

Base.@kwdef mutable struct ParamsCuprate
    # J = 0.1 * 1 / 2 * 105u"meV"
    J = 1u"meV"
    K = 0.01u"meV"
    a = 2.5u"Å"
    bareg = 1e-2Unitful.q
    g = bareg * a * cos(deg2rad(90 / 2))
    # g = 8 * 0.1^3 * Unitful.q * a
    # g = 4.4e-4u"meV / (kV / cm)"
    # Jp = 0.1 * 1 / 2 * 26.25u"meV"
    Jp = 0u"meV"
    α = 5e-2
    c = a  #6u"Å"
    η = (2K + 2 * 4J) * 0.05 # 0.3u"meV"
end
A_cupr(k, p) = 4p.J + p.K + p.Jp * (4 - 2cos(k[1] * p.a - k[2] * p.a) - 2cos(k[1] * p.a + k[2] * p.a))
B_cupr(k, p) = 2p.J * (cos(p.a * k[1]) + cos(p.a * k[2]))
ϕ_cupr(k, p) = atanh(-B_cupr(k, p) / A_cupr(k, p)) / 2
u_cupr(k, p) = cosh(ϕ_cupr(k, p))
v_cupr(k, p) = sinh(ϕ_cupr(k, p))

Ek_cupr(k, p) = sqrt(A_cupr(k, p)^2 - B_cupr(k, p)^2)
ωk_cupr(k, p) = Ek_cupr(k, p) / ħ


vertexbare(k, p) = p.g .* [
    sin(k[2] * p.a), sin(k[1] * p.a), 0
]


πq(k, q, p) = [vertexbare(k .- q, p) .* v_cupr(k, p) * v_cupr(k .+ q, p), vertexbare(q .- k, p) .* u_cupr(k, p) * u_cupr(k .+ q, p)]

Γ_cupr(q, E, p, pp) = Γ_adaptive(q, E, k -> ones(2) .* Ek_cupr(k, p), (k, q) -> 1, pp, (k, q) -> πq(k, q, p), p.a, N=2, maxevals=Int(1e4), η=1e-2u"meV")


Γ_cupr_backwards(q, E, p, pp) = Γ_adaptive_backwards(q, E, k -> ones(2) .* Ek_cupr(k, p), (k, q) -> 1, pp, (k, q) -> πq(k, q, p), p.a, N=2, maxevals=Int(1e4), η=1e-2u"meV")


##
