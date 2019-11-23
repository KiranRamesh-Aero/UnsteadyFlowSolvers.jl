#=
UI_utils.jl
    Written by Matthew White
    11/14/2019
=#
using DataStructures

struct commands
    s1::SortedDict
    s2::SortedDict

    function commands(list1,list2)
        # Stage 1
        s1 = SortedDict{String,Any}(list1[i] => "---" for i = 1:length(list1))

        #Stage 2
        s2 = SortedDict{String,Any}(list2[i] => "---" for i = 1:length(list2))

        new(s1,s2)
    end
end

function setDefaults(commandDict)
    commandDict.s1["ndiv"] = 100
    commandDict.s1["pvt"] = 0.25
    commandDict.s1["lespcrit"] = 10000
    commandDict.s1["dt"] = .015
    commandDict.s1["veldef"] = "const(1)"
    commandDict.s1["hdef"] = "const(0)"
    commandDict.s1["rho"] = 1.225
    commandDict.s1["cambercalc"] = "radial"
    commandDict.s1["oper"] = "N/A"

    commandDict.s2["wflag"] = 1
    commandDict.s2["sflag"] = 0
    commandDict.s2["wcount"] = 20

    return commandDict
end

function ask(prompt::String = "")
    # Ask a question
    print(" $prompt")
    usrCmd = strip(chomp(readline()))
end

function cmdValues(commandDict,stage)
    # display the values of all relevant variables
    c_offset = 12

    println()
    if stage == 1 #decide paramters
        dict = commandDict.s1
        exclude = ["oper"]
    elseif stage == 2
        dict = commandDict.s2
        exclude = ["ldvm","lve","save"]
    end

    for entry in dict # print values
        if sum(entry[1] .== exclude) == 0 # not on exclude list
            print("\t$(entry[1])")
            for i = 1:c_offset - length(entry[1])
                print(" ")
            end
            print(":\t")
            println("$(entry[2])")
        end
    end
    println()
end

function errorHandler(usrCmd,type::String = "")
    # Write out error messages depending on type of input expected
    #   "" if basic error message

    print("     ERROR: ")
    if length(type) > 10 # custom error message
        println(type)
    else
        if type == ""
            println("'$usrCmd' not recognized. Type HELP for a list of available commands.")
        elseif type == "many"
            println("Too many arguments have been entered.")
        else
            print("'$usrCmd' must be ")
            if type == "i"
                println("an integer.")
            elseif type == "f"
                println("a float.")
            elseif type == "s"
                println("a string")
            elseif tpye == "def"
                println("either 'eld', 'const', 'sin', or 'cos'")
            end
        end
    end
end

function svPrompt(type, prompt, subCmd = "")
    # Ask for single output variable
    if type == "i"
        fullType = Int64
    elseif type == "f"
        fullType = Float64
    end
    while true
        if subCmd == ""
            subCmd = ask(prompt)
            if lowercase(subCmd) == "exit" || lowercase(subCmd) == ""
                break
            end
        end
        if type != "s"
            outVar = tryparse(fullType,subCmd)
        end
        if outVar != nothing
            return outVar
        else
            errorHandler(subCmd,type)
        end
        subCmd = ""
    end
end

function ynAsk(prompt)
    # ask a yes or no question. Returns a boolean
    while true
        prompt = "$prompt [Y/N]: "
        usrIn = lowercase(ask(prompt))
        if usrIn == "y"
            return true
        elseif usrIn == "n"
            return false
        else
            errorHandler(usrIn,"Must be 'Y' or 'N'")
        end
    end
end
