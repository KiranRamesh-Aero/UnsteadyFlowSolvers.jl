#=
UI_main.jl
    Written by Matthew White
    11/14/2019
=#
function runUI()
    println("\33[2J") # clear command window
    commandHist = [] # User command history log
    hdef = ConstDef(0) # pre-define kinematics
    udef = ConstDef(1)
    alphadef = ConstDef(0)
    curfield = UnsteadyFlowSolvers.TwoDFlowField()
    delvort = UnsteadyFlowSolvers.delNone()
    nsteps = 0
    mat = []
    frames = []

    ## Dictonary of Cmdmands
    # Stage 1 commands
    list1 = ["airfoil","ndiv","lespcrit","alphadef","oper","pvt","runtime","hdef","veldef","dt","rho","cambercalc"]
    # Stage 2 commands
    list2 = ["wflag","sflag","wcount","lve","ldvm","save"]
    commandDict = commands(list1,list2)
    commandDict = setDefaults(commandDict)

    stage = 1 # set preliminary stage

    while true # Loop forever or unitl user exit
        if stage == 1 # Set stage prompts
            prompt = "UNSflow   c> "
            cmds = commandDict.s1
        elseif stage == 2
            prompt = "OPER      c> "
            cmds = commandDict.s2
        end
        ## User Input
        print(prompt)
        usrCmd = ask()
        commandHist = vcat(commandHist,usrCmd)

        ## check for multiple commands in one line
        if occursin(" ",usrCmd) && length(usrCmd) > 2 # sub command
            temp = split(usrCmd," ")
            if length(temp) <= 2 # extra command added
                usrCmd = lowercase(temp[1])
                subCmd = temp[2]
            else
                errorHandler(usrCmd,"many")
                continue
            end
        else # no sub command
            usrCmd = lowercase(usrCmd)
            subCmd = ""
        end

        ## Top level Commands
        if stage == 1 && usrCmd == "exit"
            break # exit program
        end
        if stage == 2 && (usrCmd == "" || usrCmd == "back")
            stage = 1
            continue # return to stage 1
        end
        if usrCmd == "?"
            cmdValues(commandDict,stage)
            continue # show variables
        end
        if usrCmd == "clc"
            println("\33[2J") # clear command window
            continue # clear screen
        end

        if stage == 1 && usrCmd == "load" # change from load command
            usrCmd = "airfoil"
        end
        if stage == 2 && (usrCmd == "ldvm" || usrCmd == "lve") # Calculate surf and curfield
            global surf = UnsteadyFlowSolvers.TwoDSurf(commandDict.s1["airfoil"],commandDict.s1["pvt"],kinem,[commandDict.s1["lespcrit"];]; ndiv = commandDict.s1["ndiv"], camberType = commandDict.s1["cambercalc"], rho = commandDict.s1["rho"])
            curfield = UnsteadyFlowSolvers.TwoDFlowField()
        end

        ## General Commands
        if haskey(cmds,usrCmd) # Test if command is available
            ## Evaluate Cmdmands
            if stage == 1 # Initilization Stage
                if usrCmd == "airfoil" # load airfoil
                    o = loadGeo(subCmd)
                    if o != nothing
                        commandDict.s1[usrCmd] = o
                        println(" Airfoil successfuly loaded.")
                    end

                elseif usrCmd == "ndiv" # Set panel divsions
                    o = svPrompt("i", "Panel Count: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "lespcrit" # Set panel divsions
                    o = svPrompt("f", "LESP Critical: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "pvt" # Set pivot location
                    o = svPrompt("f", "Pivot Point: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "dt" # Set time step
                    o = svPrompt("f", "Time Step: ", subCmd)
                    if o != nothing
                        commandDict.s1[usrCmd] = o
                        if commandDict.s1["runtime"] != "---"
                            println(" The simulation will take $(Int(ceil(commandDict.s1["runtime"]/o))) iterations.")
                        end
                    end

                elseif usrCmd == "runtime" # Set total time
                     o = svPrompt("f", "Total Time: ", subCmd)
                     if o != nothing
                         commandDict.s1[usrCmd] = o
                         println(" The simulation will take $(Int(ceil(o/commandDict.s1["dt"]))) iterations.")
                     end

                elseif usrCmd == "rho" # Set density
                    o = svPrompt("f", "Density [kg/m^3]: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "cambercalc" # Set camber calculation method
                    o = setCamberCalc(subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "alphadef" # Set alpha definition
                    alphavec = modef(subCmd)
                    if alphavec != nothing
                        alphadef = alphaCorr(alphavec)
                        param = join(alphavec[2:end], ", ")
                        commandDict.s1[usrCmd] = "$(alphavec[1])($param)"
                        temp, rt = alphadef(0)
                        commandDict.s1["runtime"] = rt
                    end

                elseif usrCmd == "hdef" || usrCmd == "veldef"# Set height and velocity definition
                    o = modef(subCmd)
                    if o != nothing
                        if o[1] == "const"
                            odef = ConstDef(o[2])
                        elseif o[1] == "lin"
                            odef = LinearDef(o[2],o[3],o[4],o[5])
                        elseif o[1] == "eld"
                            odef = EldRampReturntstartDef(o[2],o[3],o[4],o[5])
                        elseif o[1] == "sin"
                            odef = SinDef(o[2],o[3],o[4],o[5])
                        elseif o[1] == "cos"
                            odef = CosDef(o[2],o[3],o[4],o[5])
                        elseif o[1] == "gust"
                            odef = StepGustDef(o[2],o[3],o[4])
                        end
                        param = join(o[2:end], ", ")
                        commandDict.s1[usrCmd] = "$(o[1])($param)"

                        temp, rt = odef(0)

                        commandDict.s1["runtime"] = rt

                        if usrCmd == "hdef"
                            hdef = odef
                        else
                            udef = odef
                        end
                    end

                elseif usrCmd == "oper" # Change to stage 2
                    o = gotoOper(commandDict.s1)
                    if o
                        stage = 2
                        global kinem = KinemDef(alphadef, hdef, udef)
                        nsteps = Int(ceil(commandDict.s1["runtime"]/commandDict.s1["dt"]))
                    end
                end

            elseif stage == 2 # Operation Stage
                if usrCmd == "wcount" # Set step output file count
                    o = svPrompt("i","Step output file count: ", subCmd)
                    if o != nothing ; commandDict.s2[usrCmd] = o end

                elseif usrCmd == "wflag" # Write step output during sim
                    o = ynAsk("Write step output?")
                    o ? commandDict.s2[usrCmd] = 1 : commandDict.s2[usrCmd] = 0

                elseif usrCmd == "sflag" # Start from previous step output
                    o = ynAsk("Start from previous step output?")
                    o ? commandDict.s2[usrCmd] = 1 : commandDict.s2[usrCmd] = 0

                elseif usrCmd == "lve" # Run LVE
                    println()
                    println(" Running LVE Simulation...")
                    frames, IFR_field, mat = LVE(surf,curfield,nsteps,commandDict.s1["dt"],commandDict.s2["wcount"],"longwake")
                    println("\t Done!")
                    println()

                elseif usrCmd == "ldvm" # Run LDVM
                    println()
                    println(" Running LDVM Simulation...")
                    writeInterval = commandDict.s1["runtime"]/commandDict.s2["wcount"]
                    mat, surf, curfield = ldvm(surf, curfield, nsteps, commandDict.s1["dt"],commandDict.s2["sflag"], commandDict.s2["wflag"], writeInterval, delvort)
                    println("\t Done!")
                    println()

                elseif usrCmd == "save" # Set step output file count
                    if mat != [] # mat has been defined
                        saveMat(subCmd,mat)
                    else
                        errorHandler(mat,"A simulation must be ran before saving the results.")
                    end

                end
            end
        elseif usrCmd != "" # return error message
            errorHandler(usrCmd)
            continue
        end
    end
    #return commandDict
end
