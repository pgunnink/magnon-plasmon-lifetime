using Unitful

uP = u"Å"

# PRL 105, 087202 (2010):
@info 0.585u"μC/cm^2" * 90u"Å^3" / Unitful.q |> uP

# Maltz
@info Unitful.q * 1u"nm" * 0.1^3 / Unitful.q |> uP
@info 8Unitful.q * 3u"Å" * 0.1^3 / Unitful.q |> uP

# Curtis with SOC
@info 4.4e-3u"meV * cm / kV" / Unitful.q |> uP

# Curtis without SOC
@info 80 * 4.4e-3u"meV * cm / kV" / Unitful.q |> uP


