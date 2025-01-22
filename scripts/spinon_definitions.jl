
Base.@kwdef mutable struct ParamsSpinons
    χ = 1u"meV"
    a = 3u"Å"
    g = 0.1^3 * 8a * Unitful.q * 1 / 2
    Δ = 0.0
end


τnn = [
    [1, 0], [1 / 2, √(3) / 2], [-1 / 2, √(3) / 2], [-1, 0], [-1 / 2, -√(3) / 2], [1 / 2, -√(3) / 2]
]


y(i, p) = k -> exp(1im * k ⋅ τnn[i] * p.a)

AA(k, p) = p.Δ
AB(k, p) = y(1, p)(k)
AC(k, p) = y(2, p)(k) - y(5, p)(k)
AD(k, p) = 0
AE(k, p) = y(4, p)(k)
AF(k, p) = 0

BA(k, p) = AB(k, p) |> conj
BB(k, p) = p.Δ
BC(k, p) = y(3, p)(k)
BD(k, p) = -y(1, p)(k)
BE(k, p) = 0
BF(k, p) = -y(6, p)(k)

CA(k, p) = AC(k, p) |> conj
CB(k, p) = BC(k, p) |> conj
CC(k, p) = p.Δ
CD(k, p) = 0
CE(k, p) = -y(3, p)(k)
CF(k, p) = 0

DA(k, p) = AD(k, p) |> conj
DB(k, p) = BD(k, p) |> conj
DC(k, p) = CD(k, p) |> conj
DD(k, p) = -p.Δ
DE(k, p) = y(1, p)(k)
DF(k, p) = y(2, p)(k) + y(5, p)(k)

EA(k, p) = AE(k, p) |> conj
EB(k, p) = BE(k, p) |> conj
EC(k, p) = CE(k, p) |> conj
ED(k, p) = DE(k, p) |> conj
EE(k, p) = -p.Δ
EF(k, p) = y(3, p)(k)


FA(k, p) = AF(k, p) |> conj
FB(k, p) = BF(k, p) |> conj
FC(k, p) = CF(k, p) |> conj
FD(k, p) = DF(k, p) |> conj
FE(k, p) = EF(k, p) |> conj
FF(k, p) = -p.Δ


Hi(k, p) = -p.χ .* [
    AA(k, p) AB(k, p) AC(k, p) AD(k, p) AE(k, p) AF(k, p)
    BA(k, p) BB(k, p) BC(k, p) BD(k, p) BE(k, p) BF(k, p)
    CA(k, p) CB(k, p) CC(k, p) CD(k, p) CE(k, p) CF(k, p)
    DA(k, p) DB(k, p) DC(k, p) DD(k, p) DE(k, p) DF(k, p)
    EA(k, p) EB(k, p) EC(k, p) ED(k, p) EE(k, p) EF(k, p)
    FA(k, p) FB(k, p) FC(k, p) FD(k, p) FE(k, p) FF(k, p)
]


##
rot_matrix(x) =
    [cos(x) -sin(x)
        sin(x) cos(x)]

prot = rot_matrix(pi / 2)
P(i, p) = (k, E) -> p.g * (prot * τnn[i] ⋅ E[1:2]) .* exp(1im * k ⋅ τnn[i] .* p.a)
pAA(k, E, p) = 0u"meV"
pAB(k, E, p) = P(1, p)(k, E) * 1
pAC(k, E, p) = -P(2, p)(k, E) * 1 - P(5, p)(k, E) * -1
pAD(k, E, p) = 0u"meV"
pAE(k, E, p) = P(4, p)(k, E) * 1
pAF(k, E, p) = 0u"meV"

pBA(k, E, p) = pAB(k, E, p) |> conj
pBB(k, E, p) = 0u"meV"
pBC(k, E, p) = P(3, p)(k, E) * 1
pBD(k, E, p) = -P(1, p)(k, E) * -1
pBE(k, E, p) = 0u"meV"
pBF(k, E, p) = P(6, p)(k, E) * -1

pCA(k, E, p) = pAC(k, E, p) |> conj
pCB(k, E, p) = pBC(k, E, p) |> conj
pCC(k, E, p) = 0u"meV"
pCD(k, E, p) = 0u"meV"
pCE(k, E, p) = -P(3, p)(k, E) * -1
pCF(k, E, p) = 0u"meV"

pDA(k, E, p) = pAD(k, E, p) |> conj
pDB(k, E, p) = pBD(k, E, p) |> conj
pDC(k, E, p) = pCD(k, E, p) |> conj
pDD(k, E, p) = 0u"meV"
pDE(k, E, p) = P(1, p)(k, E) * 1
pDF(k, E, p) = -P(2, p)(k, E) * 1 - P(5, p)(k, E) * 1

pEA(k, E, p) = pAE(k, E, p) |> conj
pEB(k, E, p) = pBE(k, E, p) |> conj
pEC(k, E, p) = pCE(k, E, p) |> conj
pED(k, E, p) = pDE(k, E, p) |> conj
pEE(k, E, p) = 0u"meV"
pEF(k, E, p) = P(3, p)(k, E) * 1

pFA(k, E, p) = pAF(k, E, p) |> conj
pFB(k, E, p) = pBF(k, E, p) |> conj
pFC(k, E, p) = pCF(k, E, p) |> conj
pFD(k, E, p) = pDF(k, E, p) |> conj
pFE(k, E, p) = pEF(k, E, p) |> conj
pFF(k, E, p) = 0u"meV"


Pij(k, E, p) = -[
    pAA(k, E, p) pAB(k, E, p) pAC(k, E, p) pAD(k, E, p) pAE(k, E, p) pAF(k, E, p)
    pBA(k, E, p) pBB(k, E, p) pBC(k, E, p) pBD(k, E, p) pBE(k, E, p) pBF(k, E, p)
    pCA(k, E, p) pCB(k, E, p) pCC(k, E, p) pCD(k, E, p) pCE(k, E, p) pCF(k, E, p)
    pDA(k, E, p) pDB(k, E, p) pDC(k, E, p) pDD(k, E, p) pDE(k, E, p) pDF(k, E, p)
    pEA(k, E, p) pEB(k, E, p) pEC(k, E, p) pED(k, E, p) pEE(k, E, p) pEF(k, E, p)
    pFA(k, E, p) pFB(k, E, p) pFC(k, E, p) pFD(k, E, p) pFE(k, E, p) pFF(k, E, p)
]
