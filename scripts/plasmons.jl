using UnPack
using Unitful
import Unitful.me, Unitful.ħ, Unitful.ϵ0, Unitful.c




Base.@kwdef mutable struct ParamsPlasmons
    ϵ = 10
    Ef = 100u"meV"
    D = Unitful.q^2 * Ef / (π * ħ^2)
    a = 3u"Å"
    γ = 1e-2
    vF = 1e6u"m/s"
    β = vF / sqrt(2)
end


ωq_plasmons(q::Unitful.AbstractQuantity, p) = sqrt(
    1 / 2 * p.D / (2p.ϵ * ϵ0) * sqrt(4q^2 + (p.D / (2ϵ0 * c^2))^2) - 1 / 2 * (p.D / (2p.ϵ * ϵ0 * c))^2
)

# ωq_plasmons(q::Unitful.AbstractQuantity, p) = sqrt(
#     1 / 2 * p.D / (2p.ϵ * ϵ0) * sqrt(
#         4q^2 + (p.D / (2ϵ0 * c^2))^2 - 4p.β^2 * q^2 * p.ϵ / c^2
#     ) -
#     1 / 2 * (p.D / (2p.ϵ * ϵ0 * c))^2 + p.β^2 * q^2
# )
ωq_plasmons(q::AbstractVector, p) = ωq_plasmons(sqrt(q[1]^2 + q[2]^2), p)
Eq_plasmons(q, p) = ħ * ωq_plasmons(q, p)
κ(q::Unitful.AbstractQuantity, p) = sqrt(q^2 - ωq_plasmons(q, p)^2 * p.ϵ / c^2)
κ(q::AbstractVector, p) = κ(sqrt(q[1]^2 + q[2]^2), p)
Nq(q::Unitful.AbstractQuantity, p) = p.ϵ / κ(q, p)^3 * (κ(q, p)^2 + q^2)
# Nq(q::Unitful.AbstractQuantity, p) = p.ϵ / κ(q, p)^3 * (κ(q, p)^2 + q^2) + p.D / ϵ0 * p.β^2 * q^2 / (ωq_plasmons(q, p)^2 - p.β^2 * q^2)^2

Nq(q::AbstractVector, p) = Nq(sqrt(q[1]^2 + q[2]^2), p)
Efieldqdependency(q, p) = [-q[1] / sqrt(q[1]^2 + q[2]^2), -q[2] / sqrt(q[1]^2 + q[2]^2), 1im .* sqrt(q[1]^2 + q[2]^2) / κ(q, p)]
Efieldbare(q::Unitful.AbstractQuantity, p) = sqrt(
    ħ * ωq_plasmons(q, p) / (2p.a^2 * ϵ0 * Nq(q, p))
);
Efield(q::AbstractVector, p) = Efieldbare(sqrt(q[1]^2 + q[2]^2), p) .* Efieldqdependency(q, p)



##

κ(q::Unitful.AbstractQuantity, ω, p) = sqrt(q^2 - ω^2 * p.ϵ / c^2)
κ(q::AbstractVector, ω, p) = κ(sqrt(q[1]^2 + q[2]^2), p)
Nq(q::Unitful.AbstractQuantity, ω, p) = p.ϵ / κ(q, p)^3 * (κ(q, p)^2 + q^2)#+ p.D / ϵ0 * p.β^2 * q^2 / (ω^2 - p.β^2 * q^2)^2
Nq(q::AbstractVector, ω, p) = Nq(sqrt(q[1]^2 + q[2]^2), p)
Efield(q::Unitful.AbstractQuantity, ω, p) = sqrt(
    ħ * ω / (2p.a^2 * ϵ0 * Nq(q, ω, p))
);
Efield(q::AbstractVector, ω, p) = Efield(sqrt(q[1]^2 + q[2]^2), ω, p) .* Efieldqdependency(q, p)