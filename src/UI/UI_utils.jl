#=
UI_utils.jl
    Written by Matthew White
    11/14/2019
=#
struct commands
    s1::DataStructures.SortedDict
    s2::DataStructures.SortedDict

    function commands(list1,list2)
        # Stage 1
        s1 = DataStructures.SortedDict{String,Any}(list1[i] => "---" for i = 1:length(list1))

        #Stage 2
        s2 = DataStructures.SortedDict{String,Any}(list2[i] => "---" for i = 1:length(list2))

        new(s1,s2)
    end
end

function setDefaults(commandDict)
    commandDict.s1["ndiv"] = 100
    commandDict.s1["pvt"] = 0.25
    commandDict.s1["lespcrit"] = 10000
    commandDict.s1["dt"] = .015
    commandDict.s1["veldef"] = "const(1)"
    commandDict.s1["pdef"] = "const(0)"
    commandDict.s1["rho"] = 1.225
    commandDict.s1["cambercalc"] = "radial"
    commandDict.s1["oper"] = "N/A"

    commandDict.s2["wflag"] = false
    commandDict.s2["sflag"] = false
    commandDict.s2["wcount"] = 20
    commandDict.s2["anim"] = true
    commandDict.s2["animstep"] = 50

    return commandDict
end

function helpMe(stage,help1,help2)
    # Display help data based on stage
    if stage == 1 # set relevant help results
        help = help1
    elseif stage == 2
        help = help2
    end
    c_offset = 12
    println()
    println(" Available Commands:")
    for i = 1:size(help,1) # print values
        print("\t$(help[i,1])")
        for j = 1:c_offset - length(help[i,1])
            print(" ")
        end
        print(":\t")
        println("$(help[i,2])")
    end
    println()
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
        exclude = ["save","run","plot","vortplot","infoplot"]
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
        elseif type == "valid"
            println("'$usrCmd' is not a vaild input.")
        elseif type == "list"
            prompt = join(usrCmd[1:end-1], ", ")
            println("Simulation options are $prompt and $(usrCmd[end]).")
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
        else
            outVar = subCmd
        end
        if outVar != nothing
            return outVar
        else
            errorHandler(subCmd,type)
        end
        subCmd = ""
    end
end

function ynAsk(prompt,state = "")
    # Set a value to true or false
    state = lowercase(state)

    if state == "on"
        return true
    elseif state == "off"
        return false
    elseif state != ""
        errorHandler(state,"valid")
    end

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

function endTime(kindef)
    # Returns the total time for a kinematic motion to run

    if typeof(kindef) == SinDef || typeof(kindef) == CosDef
        rt = pi / kindef.k

    elseif typeof(kindef) == EldRampReturntstartDef
        fr = kindef.K/(pi*abs(kindef.amp));
        t1 = kindef.tstart
        t2 = t1 + (1. /(2*pi*fr));
        t3 = t2 + ((1/(4*fr)) - (1/(2*pi*fr)));
        t4 = t3 + (1. /(2*pi*fr));
        rt = t4+t1;

    elseif typeof(kindef) == LinearDef
        rt = kindef.tstart + kindef.len

    else
        rt = "---"
    end
    return rt
end

function iterCalc(rt,dt)
    # Calculate iteration count and time vector
    iter = Int(ceil(rt/dt))
    trange = collect(0:dt:iter*dt)
    return iter, trange
end

function kinrange(kindef,trange)
    # calculate kinematic output at each time step
    kinem = [kindef(i) for i in trange]
    return kinem
end

function airfoilGeometry(coord_file,ndiv,camberType)
    # Extract geometry coordinates and calculate camber line

    theta = zeros(ndiv)
    x = zeros(ndiv)
    cam = zeros(ndiv)
    cam_slope = zeros(ndiv)

    if (coord_file != "FlatPlate")
        A = DelimitedFiles.readdlm(coord_file, Float64); # actual geometry
        c = maximum(A[:,1]) - minimum(A[:,1])
    else
        c = 1
    end


    if camberType == "radial"
        dtheta = pi/(ndiv-1)
        for ib = 1:ndiv
            theta[ib] = (ib-1.)*dtheta
            x[ib] = c/2. *(1-cos(theta[ib]))
        end
    elseif camberType == "linear"
        dx = c / (ndiv-1)
        for ib = 2:ndiv
            x[ib] = x[ib-1] + dx
        end
    end

    if (coord_file != "FlatPlate")
        cam, cam_slope = camber_calc(x, coord_file) # camber line
    else
        A = hcat(x,cam)
    end
    C = hcat(x,cam)
    return A,C
end

function makeInitPlot(cmds,trange,alphadef,pdef,udef)
    # Create all inputs for intitilization plot and plot results

    ## Airfoil
    if cmds["airfoil"] != "---"
        A, C = airfoilGeometry(cmds["airfoil"],cmds["ndiv"],cmds["cambercalc"])
    else
        A = 0
        C = 0
    end
    ## Kinem ranges
    arange = kinrange(alphadef,trange) .* 180 ./ pi
    hrange = kinrange(pdef,trange)
    urange = kinrange(udef,trange)

    ## Plot
    clf()
    initPlot(A,C,trange,arange,hrange,urange)
end
