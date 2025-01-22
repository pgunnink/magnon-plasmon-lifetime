using PyPlot
using Unitful
using Revise
using HCubature
using ProgressMeter
includet("plasmons.jl")
includet("EP_functions.jl")
include("style.jl")
includet("MnF2_3D_def.jl")


p = ParamsMnF()
pp = ParamsPlasmons(Ef=10u"meV", ϵ=30)


## iterative scheme
function iterative_Γ(p, pp; Nk=20, Nq=50, qlim=300u"μm^-1", η=1e-3u"meV")
    Nk_half = Int(Nk / 2)
    mirror_i = i -> Nk_half + (Nk_half - i) + 1
    krange = range(-pi / p.a, stop=pi / p.a, length=Nk)
    qrange = range(-qlim, stop=qlim, length=Nq)

    magnon_selfenergies = zeros(typeof(1.0u"meV"), 2, Nk, Nk, Nk)
    plasmon_selfenergies = zeros(typeof(1.0u"meV"), Nq, Nq)


    plasmon_selfenergies_new = zeros(typeof(1.0u"meV"), Nq, Nq)
    magnon_selfenergies_new = zeros(typeof(1.0u"meV"), 2, Nk, Nk, Nk)


    # first do the plasmon self-energies

    @showprogress for (i_q, qx) in enumerate(qrange)
        for (j_q, qy) in enumerate(qrange)
            q = [qx, qy]
            plasmon_energy = Eq_plasmons(q, pp) - 1im * plasmon_selfenergies[i_q, j_q]
            plasmon_field = Efield(q, plasmon_energy / ħ, pp)
            res = 0u"meV"

            for (i_k, kx) in enumerate(krange)
                for (j_k, ky) in enumerate(krange)
                    for (k_k, kz) in enumerate(krange)
                        k = [kx, ky, kz]
                        # mind the conjugation here!
                        Ek = Ek_MnF(k, p) .- 1im * magnon_selfenergies[:, i_k, j_k, k_k]
                        Ekq = Ek_MnF(.-k, p) .- 1im * magnon_selfenergies[:, i_k|>mirror_i, j_k|>mirror_i, k_k|>mirror_i]
                        res += imag(abs(plasmon_field ⋅ πk_independent(k, p)[1])^2 / (plasmon_energy - Ek[1] - Ekq[2] + η * 1im))
                        res += imag(abs(plasmon_field ⋅ πk_independent(k, p)[2])^2 / (plasmon_energy - Ek[2] - Ekq[1] + η * 1im))
                    end
                end
            end
            plasmon_selfenergies_new[i_q, j_q] = p.Nz * res / Nq^2
        end
    end

    # then use these to calculate the magnon self-energies

    @showprogress for (i_k, kx) in enumerate(krange)
        for (j_k, ky) in enumerate(krange)
            for (k_k, kz) in enumerate(krange)
                k = [kx, ky, kz]
                # mind the conjugation here!
                Ek = Ek_MnF(k, p) .- 1im * magnon_selfenergies[:, i_k, j_k, k_k]
                res_α = 0u"meV"
                res_β = 0u"meV"

                for (i_q, qx) in enumerate(qrange)
                    for (j_q, qy) in enumerate(qrange)
                        q = [qx, qy]
                        plasmon_energy = Eq_plasmons(q, pp) - 1im * plasmon_selfenergies_new[i_q, j_q]

                        plasmon_field = Efield(q, plasmon_energy / ħ, pp)


                        Ekq = Ek_MnF(k, p) .- 1im * magnon_selfenergies[:, i_k, j_k, k_k]
                        res_α += imag(
                            abs(plasmon_field ⋅ πk_independent(k, p)[2])^2 / (Ek[1] - Ekq[2] - plasmon_energy + 1im * η)
                        )
                        res_β += imag(
                            abs(plasmon_field ⋅ πk_independent(k, p)[1])^2 / (Ek[2] - Ekq[1] - plasmon_energy + 1im * η)
                        )
                    end
                end
                magnon_selfenergies_new[:, i_k, j_k, k_k] .= [res_α, res_β]
            end
        end
    end
    plasmon_selfenergies_new, magnon_selfenergies_new
end

##
pE, mE = iterative_Γ(p, pp)