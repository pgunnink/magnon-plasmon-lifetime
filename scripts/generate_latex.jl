using DrWatson
@quickactivate
using Printf
include(scriptsdir("cuprate_def.jl"))
include(scriptsdir("MnF2_3D_def.jl"))

function generate_latex_command(name, cmd)
    "\\newcommand{\\$(name)}{$cmd}"
end

function export_params(p, processing_dict)
    for x in keys(processing_dict)
        println(processing_dict[x](getfield(p, x)))
    end
end

function convert_e_to_num_latex(x)
    exp_string = @sprintf "%.0e" x
    return exp_string[2:end]
end

function export_params_ParamsMnF(p=ParamsMnF())
    processing_dict = Dict(
        :J => x -> generate_latex_command("Jexchange", "\\SI{$(round(x |> u"meV" |> ustrip; digits=2))}{meV}"),
        :K => x -> generate_latex_command("Kani", "\\SI{$(round(x |> u"meV" |> ustrip; digits=2))}{meV}"),
        :bareg => x -> generate_latex_command("gpolarizationrutile", "\\num{$(@sprintf "%.0e" x / Unitful.q |> upreferred)}"))
    export_params(p, processing_dict)
end

function export_params_Cuprate(p=ParamsCuprate())
    processing_dict = Dict(
        :J => x -> generate_latex_command("Jexchangecupr", "\\SI{$(round(x |> u"meV" |> ustrip; digits=2))}{meV}"),
        :K => x -> generate_latex_command("Kanicupr", "\\SI{$(round(x |> u"meV" |> ustrip; digits=2))}{meV}"),
        :bareg => x -> generate_latex_command("gpolarizationcuprate", "\\num{$(@sprintf "%.0e" x / Unitful.q |> upreferred)}"))
    export_params(p, processing_dict)
end

function export_params_plasmons(p=ParamsPlasmons())
    processing_dict = Dict(
        :Ef => x -> generate_latex_command("Efermi", "\\SI{$(x |> u"meV" |> ustrip)}{meV}"),
        :ϵ => x -> generate_latex_command("epsilonplasmon", "\\num{$x}"),
        :γ => x -> generate_latex_command("gammaplasmon", "\\num{$(convert_e_to_num_latex(x))}"),
    )
    export_params(p, processing_dict)
end


export_params_Cuprate()
export_params_plasmons()
export_params_ParamsMnF()