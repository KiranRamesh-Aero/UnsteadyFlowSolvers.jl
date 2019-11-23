#=
s2Cmds.jl
    Written by Matthew White
    11/16/2019
=#
using DelimitedFiles

function saveMat(name,mat)
    # Save the output matrix
    exts = ["dat","txt"] # valid file extensions
    while true
        if name == ""
            name = ask("File Name: ")
        end
        if lowercase(name) == "exit" || lowercase(name) == ""
            break
        end
        parts = split(name,".")
        if occursin("/",name) # can change in the future to look for available folders
            errorHandler(name,"$name must not be inside a folder.")

        elseif length(parts) > 2 # extra periods
            errorHandler(name,"$name is an invalid file name.")

        elseif sum(parts[end] .== exts) == 0 # not a valid extention
            errorHandler(name,"$(parts[end]) is an invalid extension.")

        else
            writedlm(name,mat)
            break
        end
        name = ""
    end
end
