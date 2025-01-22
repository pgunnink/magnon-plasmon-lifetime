using Printf

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
