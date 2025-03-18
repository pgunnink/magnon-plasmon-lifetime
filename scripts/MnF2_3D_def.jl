using Unitful
using Printf
# He = 520kOe\
# 1 Oe = 0.0001T
HE = 520e3 * 0.0001u"T"
HA = 8.2e3 * 0.0001u"T"
Base.@kwdef mutable struct ParamsMnF
    # J = 1u"meV"
    J = HE / (8 / (2Unitful.μB)) |> u"meV"
    δJ = 0 * 0.05u"meV"
    # K = 1u"meV"
    K = HA / (2 / (2Unitful.μB)) |> u"meV"
    a = 5u"Å"
    c = 3u"Å"
    # g = 8 * 0.1^3 * Unitful.q * cos(deg2rad(128 / 2))
    bareg = 1e-2Unitful.q
    g = bareg * cos(deg2rad(119 / 2))
    Nz = 1
    α = 5e-2
    η = (2K + 2 * 8J) * 0.05 #0.5u"meV"
end
A_MnF(k, p) = 8p.J + p.K
B_MnF(k, p) = 8p.J * cos(k[1] * p.a / 2) * cos(k[2] * p.a / 2) * cos(k[3] * p.a / 2)
δA_MnF(k, p) = 2p.δJ * (cos(k[1] * p.a) - cos(k[3] * p.a))

ϕ_MnF(k, p) = atanh(-B_MnF(k, p) / A_MnF(k, p)) / 2
u_MnF(k, p) = cosh(ϕ_MnF(k, p))
v_MnF(k, p) = sinh(ϕ_MnF(k, p))

Ek_MnF(k, p) = sqrt(A_MnF(k, p)^2 - B_MnF(k, p)^2) * [1, 1] .+ δA_MnF(k, p) .* [1, -1]
ωk_MnF(k, p) = Ek_MnF(k, p) / ħ


##


# 100 orientation
vertexbare(k, p) = p.g .* [
    -2im * p.a * sin(0.5k[1] * p.a) * sin(0.5k[2] * p.a) * sin(0.5k[3] * p.a),
    -0.45im * p.a * cos(0.5k[1] * p.a) * cos(0.5k[2] * p.a) * sin(0.5k[3] * p.a),
    -0.45im * p.a * cos(0.5k[1] * p.a) * sin(0.5k[2] * p.a) * cos(0.5k[3] * p.a)
]



πq(k, q, p) = [vertexbare(k .- q, p) .* v_MnF(k, p) * v_MnF(k .+ q, p), vertexbare(q .- k, p) .* u_MnF(k, p) * u_MnF(k .+ q, p)]




πk_independent(k, p) = [vertexbare(k, p) .* v_MnF(k, p) * v_MnF(k, p), vertexbare(.-k, p) .* u_MnF(k, p) * u_MnF(k, p)]

Γ_MnF(q, E, p, pp) = Γ_adaptive(q, E, k -> Ek_MnF(k, p), (k, q) -> 1, pp, (k, q) -> πq(k, q, p), p.a, N=2, maxevals=Int(1e4), η=0.0001u"meV")
