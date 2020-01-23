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

function runSim(commandDict,surf,curfield,nsteps,delvort,sim)
    # Run the specificed simulation, returns mat
    sim = lowercase(sim)
    simList = ["lve","ldvm"]
    while true
        if sim == "lve"
            println()
            println(" Running LVE Simulation...")
            writeInterval = floor(nsteps/commandDict.s2["wcount"])
            frames, IFR_field, mat = LVE(surf,curfield,nsteps,commandDict.s1["dt"],commandDict.s2["wcount"],"UI Window",1000,commandDict.s2["anim"],30,commandDict.s2["wflag"],writeInterval)
            println("\t Done!")
            println()
            return mat
        elseif sim == "ldvm"
            println()
            println(" Running LDVM Simulation...")
            writeInterval = commandDict.s1["runtime"]/commandDict.s2["wcount"]
            commandDict.s2["wflag"] ? wflag = 1 : wflag = 0
            commandDict.s2["sflag"] ? sflag = 1 : sflag = 0
            mat, surf, curfield = ldvm(surf, curfield, nsteps, commandDict.s1["dt"],sflag, wflag, writeInterval, delvort)
            println("\t Done!")
            println()
            return mat
        end

        sim = ask("Simulation Type: ")
        if sim == ""
            break
        elseif sum(simList .== sim) == 0 # if not on simulation list
            errorHandler(simList,"lsit")
        end
    end
end

function plotResults(mat,subCmd)
    # Plot anything from sim results
    ploty = []
    plotCnt = []
    plotx = []
    matCol = ["t","alpha","h","u","lesp","cl","cd","cm"]
    matIdx = DataStructures.SortedDict{String,Any}(matCol[i] => i for i = 1:length(matCol))
    if haskey(matIdx,lowercase(subCmd)) # subcmd is plot type (infers only one plot to be calculated)
        plotCnt = 1
        ploty = [subCmd]
    else
        plotCnt = tryparse(Int64,subCmd)
        if plotCnt == nothing # subcmd is number of plots
            plotCnt = []
        end
    end
    # Gather plot info
    while true
        if plotCnt == [] # number of plots
            plotCnt = svPrompt("i","Number of plots: ")
            if plotCnt == nothing
                return nothing
            end
        end
        if plotCnt > 7 # too many plots
            plotCnt = []
            errorHandler(plotCnt,"$plotCnt exceeds the maximum number of plots.")
            continue
        else
            break
        end
    end
    if ploty == [] # determine y axis(es)
        for i = 1:plotCnt
            while true
                ploty = [ploty;svPrompt("s","Plot $i y-axis: ")]
                if ploty[i] == nothing
                    return nothing
                end
                if !haskey(matIdx,lowercase(ploty[i])) # not a vaild y axis
                    errorHandler(matCol,"list")
                    ploty = ploty[1:end-1] # reset ploty
                    continue
                else
                    break
                end
            end
        end
    end
    while true # determine x axis
        plotx = svPrompt("s","X axis: ")
        if plotx == nothing
            return nothing
        end
        if !haskey(matIdx,lowercase(plotx)) # not a vaild x axis
            errorHandler(matCol,"list")
            continue
        else
            break
        end
    end
    plotMat(mat,plotCnt,plotx,ploty,matIdx)

    if ynAsk("Save figure?")
        fn = svPrompt("s","File name: ")
        if occursin(".",fn)
            temp = split(fn,".") # pyplot only supports .png
            fn = temp[1]
        end
        fn = "$fn.png"
        savefig(fn)
    end
end

function plotTime(tt  = "")
    # Asks for a specific time to plot step file results
    ttparse = tryparse(Float64, tt)
    if tt == "all"
        return tt
    elseif tt == "" || ttparse == nothing
        if ynAsk("Plot all step files?")
            return "all"
        else
            return svPrompt("f","Time to plot: ")
        end
    else
        return ttparse
    end
end
