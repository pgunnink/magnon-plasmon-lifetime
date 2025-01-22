using Parameters
using UnPack
using HCubature
using Unitful
import Unitful.ħ, Unitful.me, Unitful.ϵ0, Unitful.k, Unitful.μB, Unitful.μ0
using ProgressMeter

##
@with_kw mutable struct ParamsYIG
    a = 12.376u"Å"
    μ = 2μB
    Ms = 180u"mT" / μ0
    ρex = 5e-17u"T*m^2"
    S = Ms * a^3 / μ |> upreferred
    J = ρex * μB / a^2 / S
    Δ = μ0 * μ * Ms
    H = 700 * 1e-4u"T"
    h = μ * H
    d = 250u"nm"
    # d = 400a
    T = 1u"K"
    α = 1e-4
    Js = J * S * a^2 |> u"meV * Å^2"
    E0 = sqrt(h * (h + Δ))
end
p = ParamsYIG()
fk(k, p) = (1 - exp(-k * p.d)) / (1e-10 + k * p.d)
Dxx(k, θ, p) = μ0 * p.μ^2 / p.a^3 * (1 / 3 - fk(k, p))
Dyy(k, θ, p) = μ0 * p.μ^2 / p.a^3 * (1 / 3 - sin(θ)^2 * fk(k, p))
Dzz(k, θ, p) = μ0 * p.μ^2 / p.a^3 * (1 / 3 - cos(θ)^2 * fk(k, p))

##
function Ak(k, θ, p)
    @unpack h, J, S, a, Δ = p
    h + J * S * (4 - 2cos(k * cos(θ) * a) - 2cos(k * sin(θ) * a)) - S / 2 * (Dxx(k, θ, p) + Dyy(k, θ, p)) + Δ / 3
end
Bk(k, θ, p) = -p.S / 2 * (Dxx(k, θ, p) - Dyy(k, θ, p))
Jz(q, p) = p.J * cos(q * p.a)

Ek(k, θ, p) = sqrt(Ak(k, θ, p)^2 - Bk(k, θ, p)^2)
# Ek(k, θ, p) = sqrt(
#     (p.h + p.Js * k^2 + p.Δ * (1 - fk(k, p)) * sin(θ)^2) *
#     (p.h + p.Js * k^2 + p.Δ * fk(k, p))
# )

##
fBE(x, p) = 1 / (exp(x / (k * p.T)) - 1)
# define units:
const E_unit = u"μeV"
const χ_unit = u"μeV^-1"
const k_unit = u"μm^-1"
function integrand_full(ω, q, p)
    ωk = x -> Ek(sqrt(x[1]^2 + x[2]^2) * k_unit, atan(x[2], x[1]), p)
    x -> ustrip.(χ_unit, (2pi)^-2 * (fBE(ωk(x), p) - fBE(ωk(x .+ ustrip.(k_unit, q)), p)) / (ħ * ω * (1 + 1im * p.α) + ωk(x) - ωk(x .+ ustrip.(k_unit, q))))

end
function χ0_manual_integration(q, ω, p)
    # first estimate the thermal wavelength
    f = integrand_full(ω, q, p)
    Nk = 2000
    k_range = range(-0.1pi / p.a |> k_unit |> ustrip, stop=0.1pi / p.a |> k_unit |> ustrip, length=Nk)
    dk = (k_range[2] - k_range[1]) * k_unit
    prefactor = (2pi)^-2 * p.a^2 * dk^2

    res = 0.0 + 0.0im
    for kx in k_range
        for ky in k_range
            res += prefactor * f([kx, ky])
        end
    end
    return res
end

function integrand_isotropic_full(ω, q, p)
    ωk = x -> p.E0 + p.J * p.S * (4 - 2cos(x[1] * k_unit * p.a) - 2cos(x[2] * k_unit * p.a))
    x -> ustrip.(χ_unit, (2pi)^-2 * (fBE(ωk(x), p) - fBE(ωk(x .+ ustrip.(k_unit, q)), p)) / (ħ * ω * (1 + 1im * p.α) + ωk(x) - ωk(x .+ ustrip.(k_unit, q))))

end


function χ0_iso_manual_integration(q, ω, p)
    # first estimate the thermal wavelength
    f = integrand_isotropic_full(ω, q, p)
    Nk = 2000
    k_range = range(-0.1pi / p.a |> k_unit |> ustrip, stop=0.1pi / p.a |> k_unit |> ustrip, length=Nk)
    dk = (k_range[2] - k_range[1]) * k_unit
    prefactor = (2pi)^-2 * p.a^2 * dk^2

    res = 0.0 + 0.0im
    for kx in k_range
        for ky in k_range
            res += prefactor * f([kx, ky])
        end
    end
    return res
end

function χ0_isotropic_full(q, ω, p; rtol=1e-4)
    # first estimate the thermal wavelength
    Λ = 1 / sqrt(p.J * p.S * p.a^2 / (k * p.T))
    if Λ > pi / p.a
        Λ = pi / 2p.a
    end
    return hcubature(integrand_isotropic_full(ω, q, p), -ones(2) * 2Λ .|> x -> ustrip.(k_unit, x), ones(2) * 2Λ .|> x -> ustrip.(k_unit, x), rtol=rtol)
end


