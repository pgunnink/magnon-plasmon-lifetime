using LinearAlgebra
using Unitful
using HCubature

L(x, y, η=1e-2u"meV") = 1 / π * η / ((x - y)^2 + (η)^2)

function dos2(q, E, f, sample_grid; η=1e-2u"meV")
    # f should is a function which takes a wave vector and returns an energy
    res = 0.0u"meV^-1"
    for k in sample_grid
        res += L(E, f(k) + f(q .- k), η)
    end
    return res ./ length(sample_grid)
end


function dos2_square(q, E, f, a; Nk=100, η=1e-2u"meV")
    krange = range(-pi / a, stop=pi / a, length=Nk)
    sample_grid = Iterators.product(krange, krange)
    dos2(q, E, f, sample_grid; η=η)
end

function dos2_adaptive(q, E, f, a; rtol=1e-3, η=1e-2u"meV", maxevals=Int(1e3), N=1, d=2)
    kunit = u"nm^-1"
    resunit = u"meV^-1"
    if N == 1
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit) + f(q .- k .* kunit), η) / (2pi)^d
        )
    elseif N == 2
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit)[1] + f(q .- k .* kunit)[2], η) / (2pi)^d + L(E, f(k .* kunit)[2] + f(q .- k .* kunit)[1], η) / (2pi)^d
        )
    else
        @error "N=$N not supported"
    end
    lim = ustrip(kunit, pi / a)
    hcubature(g, -lim .* ones(d), lim .* ones(d), rtol=rtol, maxevals=maxevals) .* resunit
end

function Szplus_dos2_adaptive(q, E, f, a; rtol=1e-3, η=1e-2u"meV", maxevals=Int(1e3), N=1, d=2)
    kunit = u"nm^-1"
    resunit = u"meV^-1"
    if N == 1
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit) + f(q .- k .* kunit), η) / (2pi)^d
        )
    elseif N == 2
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit)[1] + f(q .- k .* kunit)[1], η) / (2pi)^d
        )
    else
        @error "N=$N not supported"
    end
    lim = ustrip(kunit, pi / a)
    hcubature(g, -lim .* ones(d), lim .* ones(d), rtol=rtol, maxevals=maxevals) .* resunit
end

function Szmin_dos2_adaptive(q, E, f, a; rtol=1e-3, η=1e-2u"meV", maxevals=Int(1e3), N=1, d=2)
    kunit = u"nm^-1"
    resunit = u"meV^-1"
    if N == 1
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit) + f(q .- k .* kunit), η) / (2pi)^d
        )
    elseif N == 2
        g = k -> ustrip(resunit,
            L(E, f(k .* kunit)[2] + f(q .- k .* kunit)[2], η) / (2pi)^d
        )
    else
        @error "N=$N not supported"
    end
    lim = ustrip(kunit, pi / a)
    hcubature(g, -lim .* ones(d), lim .* ones(d), rtol=rtol, maxevals=maxevals) .* resunit
end






##
function dos2_two_mode(q, E, f, sample_grid)
    # f should is a function which takes a wave vector and returns an energy
    res = 0.0u"meV^-1"
    for k in sample_grid
        res += L(E, f(k)[1] + f(q .- k)[2])
        res += L(E, f(k)[2] + f(q .- k)[1])
    end
    return res ./ length(sample_grid)
end
function dos2_square_two_mode(q, E, f, a; Nk=100)
    krange = range(-pi / a, stop=pi / a, length=Nk)
    sample_grid = Iterators.product(krange, krange)
    dos2_two_mode(q, E, f, sample_grid)
end
##

const Γη = 1e-1u"meV"
function Γ(q, E, Ek, vertex, pp, πk, sample_grid; η=Γη)
    res = 0.0im * u"meV"
    for k in sample_grid
        res += (vertex(k, q) .* Efield(q, E / ħ, pp) ⋅ πk(k, q))^2 / (E - Ek(k) - Ek(q .- k) + η * 1im)
    end
    res / length(sample_grid)
end


function Γ_adaptive(q, E, Ek, vertex, pp, πk, a; η=Γη, rtol=1e-3, atol=1e-4, maxevals=Int(1e3), N=1)
    kunit = u"nm^-1"
    resunit = u"meV"
    d = length(q)
    if N == 1
        @error "not working right now"
        g = k -> begin
            ustrip(resunit,
                1 / (E - (Ek(k .* kunit) + Ek(q .- k .* kunit)) * (1 - 1im * p.α) + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, q))^2 / (2pi)^d)
        end
    elseif N == 2
        g = k -> begin
            AB = ustrip(resunit,
                1 / (E - (Ek(k .* kunit)[1] + Ek(q .- k .* kunit)[2]) + 1im * p.η + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, q)[1])^2 / (2pi)^d)
            BA = ustrip(resunit,
                1 / (E - (Ek(k .* kunit)[2] + Ek(q .- k .* kunit)[1]) + 1im * p.η + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, q)[2])^2 / (2pi)^d)
            AB + BA
        end

    else
        @error "N= $N not supported"
    end
    lim = ustrip(kunit, pi / a)
    hcubature(g, -lim .* ones(d), lim .* ones(d), rtol=rtol, atol=atol, maxevals=maxevals) .* resunit
end


function Γ_adaptive_backwards(q, E, Ek, vertex, pp, πk, a; η=Γη, rtol=1e-3, atol=1e-4, maxevals=Int(1e3), N=1)
    kunit = u"nm^-1"
    resunit = u"meV"
    d = length(q)
    if N == 1
        @error "not working right now"
        g = k -> begin
            ustrip(resunit,
                1 / (E - (Ek(k .* kunit) + Ek(q .- k .* kunit)) * (1 - 1im * p.α) + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, q))^2 / (2pi)^d)
        end
    elseif N == 2
        g = k -> begin
            AB = ustrip(resunit,
                -1 / (E + (Ek(k .* kunit)[1] + Ek(-q .- k .* kunit)[2]) * (1 + 1im * p.α) + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, -q)[1])^2 / (2pi)^d)
            BA = ustrip(resunit,
                -1 / (E + (Ek(k .* kunit)[2] + Ek(-q .- k .* kunit)[1]) * (1 + 1im * p.α) + η * 1im) * abs(Efield(q, E / ħ, pp) ⋅ πk(k .* kunit, -q)[2])^2 / (2pi)^d)
            1 / 2 * (AB + BA)
        end

    else
        @error "N= $N not supported"
    end
    lim = ustrip(kunit, pi / a)
    hcubature(g, -lim .* ones(d), lim .* ones(d), rtol=rtol, atol=atol, maxevals=maxevals) .* resunit
end


function Γ_square(q, E, Ek, vertex, pp, πk, a; Nk=100, η=Γη)
    krange = range(-pi / a, stop=pi / a, length=Nk)
    sample_grid = Iterators.product(krange, krange)
    Γ(q, E, Ek, vertex, pp, πk, sample_grid; η=η)
end


spectral(x, y, η) = -1 / π * η / ((x - y)^2 + η^2)
spectral_plasmons(E, q, pp, Σ) = spectral(E, Eq_plasmons(q, pp) + real(Σ), -imag(Σ) + pp.γ * Eq_plasmons(q, pp))