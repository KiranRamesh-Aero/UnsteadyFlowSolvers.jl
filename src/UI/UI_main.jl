#=
UI_main.jl
    Written by Matthew White
    11/14/2019
=#
function runUI()
    println("\33[2J") # clear command window
    commandHist = [] # User command history log
    pdef = ConstDef(0) # pre-define kinematics
    udef = ConstDef(1)
    alphadef = ConstDef(0)
    alphadef = ConstDef(0)
    curfield = UnsteadyFlowSolvers.TwoDFlowField()
    delvort = UnsteadyFlowSolvers.delNone()
    nsteps = 0 # pre-define global vars
    trange = collect(1:10)
    lockTime = false
    mat = []
    frames = []

    ## Dictonary of Cmdmands
    # Stage 1 commands
    list1 = ["airfoil","ndiv","lespcrit","alphadef","oper","pvt","runtime","pdef","veldef","dt","rho","cambercalc","motionfile"]
    help1 = ["Define airfoil designation: coordinate file or FlatPlate","Define airfoil panel division count","Define critical LESP value for LEV shedding","Define pitch motion parameters","Change to operation phase","Define airfoil pivot location","Define simulation run time","Define plunge motion parameters","Define velocity profile parameters","Define time step","Define operating density","Define camber calculation method","Use a file to describe the motion"]
    help1 = hcat(list1,help1)
    help1[1,1] = "load" # overwrite to show command
    help1 = sortslices(help1,dims = 1)

    # Stage 2 commands
    list2 = ["wflag","sflag","wcount","run","save","anim","animstep","plot","vortplot","infoplot"]
    help2 = ["If enabled, writes out iteration files to default save location","If enabled, start simulation from previous output files","Define number of output files","Run a simulation","Save the polar data to a file","If enabled, shows an animation of polars and motion during LVE simulation","Define step count between animation frames","Plot simulation output data","Plot vortex location data for each write step (Must have wflag enabled)","Plot step parameters for each write step (Must have wflag enabled)"]
    help2 = sortslices(hcat(list2,help2),dims = 1)

    commandDict = commands(list1,list2)
    commandDict = setDefaults(commandDict)

    stage = 1 # set preliminary stage

    # Run 'HELP' on startup
    helpMe(stage,help1,help2)

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
            close("all")
            dirvec = readdir() # Clear step files if prompted
            if "Step Files" in dirvec
                if ynAsk("Clear step file cache?")
                    cleanWrite()
                end
            end
            break # exit program
        end
        if stage == 2 && (usrCmd == "" || usrCmd == "back")
            stage = 1
            close("all")
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
        if usrCmd == "help"
            helpMe(stage,help1,help2)
            continue
        end

        if stage == 1 && usrCmd == "load" # change from load command
            usrCmd = "airfoil"
        end
        if stage == 2 && (usrCmd == "run") # Calculate surf and curfield
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
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "ndiv" # Set panel divsions
                    o = svPrompt("i", "Panel Count: ", subCmd)
                    if o != nothing
                        commandDict.s1[usrCmd] = o
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "lespcrit" # Set panel divsions
                    o = svPrompt("f", "LESP Critical: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "pvt" # Set pivot location
                    o = svPrompt("f", "Pivot Point: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "dt" # Set time step
                    if !lockTime
                        o = svPrompt("f", "Time Step: ", subCmd)
                        if o != nothing
                            commandDict.s1[usrCmd] = o
                            if commandDict.s1["runtime"] != "---"
                                nsteps, trange = iterCalc(commandDict.s1["runtime"],o)
                                println(" The simulation will take $nsteps iterations.")
                            end
                        end
                    else
                        errorHandler(1,"The file '$(commandDict.s1["motionfile"])' prevents this from being changed.")
                    end

                elseif usrCmd == "runtime" # Set total time
                    if !lockTime
                         o = svPrompt("f", "Total Time: ", subCmd)
                         if o != nothing
                             commandDict.s1[usrCmd] = o
                             nsteps, trange = iterCalc(o,commandDict.s1["dt"])
                             println(" The simulation will take $nsteps iterations.")
                         end
                     else
                         errorHandler(1,"The file '$(commandDict.s1["motionfile"])' prevents this from being changed.")
                    end

                elseif usrCmd == "rho" # Set density
                    o = svPrompt("f", "Density [kg/m^3]: ", subCmd)
                    if o != nothing ; commandDict.s1[usrCmd] = o end

                elseif usrCmd == "cambercalc" # Set camber calculation method
                    o = setCamberCalc(subCmd)
                    if o != nothing
                        commandDict.s1[usrCmd] = o
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "motionfile" # load time, alphadef, h, and vel from file
                    t,a,h,u,fn = fileInput(subCmd)
                    if t != nothing
                        commandDict.s1["motionfile"] = fn

                        alphadef = FileDef(t,a)
                        commandDict.s1["alphadef"] = "file"

                        pdef = FileDef(t,h)
                        commandDict.s1["pdef"] = "file"

                        udef = FileDef(t,u)
                        commandDict.s1["veldef"] = "file"

                        commandDict.s1["dt"] = t[2]-t[1]
                        commandDict.s1["runtime"] = t[end]
                        lockTime = true # might need a way to toggle off if file is completely overrided

                        trange = t
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "alphadef" # Set alpha definition
                    alphavec = modef(subCmd)
                    if alphavec != nothing
                        alphadef = alphaCorr(alphavec)
                        param = join(alphavec[2:end], ", ")
                        commandDict.s1[usrCmd] = "$(alphavec[1])($param)"
                        rt = endTime(alphadef)
                        if rt != "---"
                            if  commandDict.s1["runtime"] == "---" || rt > commandDict.s1["runtime"] && !lockTime
                                commandDict.s1["runtime"] = rt
                                nsteps, trange = iterCalc(commandDict.s1["runtime"],commandDict.s1["dt"])
                            end
                        end
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "pdef" || usrCmd == "veldef"# Set height and velocity definition
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

                        rt = endTime(odef)

                        if rt != "---"
                            if commandDict.s1["runtime"] == "---" || rt > commandDict.s1["runtime"] && !lockTime
                                commandDict.s1["runtime"] = rt
                                nsteps, trange = iterCalc(commandDict.s1["runtime"],commandDict.s1["dt"])
                            end
                        end

                        if usrCmd == "pdef"
                            pdef = odef
                        else
                            udef = odef
                        end
                        makeInitPlot(commandDict.s1,trange,alphadef,pdef,udef)
                    end

                elseif usrCmd == "oper" # Change to stage 2
                    o = gotoOper(commandDict.s1)
                    if o
                        stage = 2
                        global kinem = KinemDef(alphadef, pdef, udef)
                        nsteps = Int(ceil(commandDict.s1["runtime"]/commandDict.s1["dt"]))
                    end
                    close()
                end

            elseif stage == 2 # Operation Stage
                if usrCmd == "wcount" # Set step output file count
                    o = svPrompt("i","Step output file count: ", subCmd)
                    if o != nothing ; commandDict.s2[usrCmd] = o end

                elseif usrCmd == "wflag" # Write step output during sim
                    commandDict.s2[usrCmd] = ynAsk("Write step output?",subCmd)

                elseif usrCmd == "sflag" # Start from previous step output
                    commandDict.s2[usrCmd] = ynAsk("Start from previous step output?",subCmd)

                elseif usrCmd == "anim" # Set LVE animation state
                    commandDict.s2[usrCmd] = ynAsk("Display LVE Animation?", subCmd)

                elseif usrCmd == "animstep" # set animation step
                    o = svPrompt("i", "Animation Step: ", subCmd)
                    if o != nothing ; commandDict.s2[usrCmd] = o end

                elseif usrCmd == "run" # Run simulation
                    if commandDict.s2["wflag"]
                        cleanWrite()
                    end
                    o = runSim(commandDict,surf,curfield,nsteps,delvort,subCmd)
                    if o != nothing ; mat = o end

                elseif usrCmd == "plot"
                    plotResults(mat,subCmd)

                elseif usrCmd == "save" # Set step output file count
                    if mat != [] # mat has been defined
                        saveMat(subCmd,mat)
                    else
                        errorHandler(mat,"A simulation must be ran before saving the results.")
                    end

                elseif usrCmd == "infoplot"
                    if mat != []
                        tt = plotTime(subCmd)
                        makeKinemClVortPlots2D(mat,tt)
                    else
                        errorHandler(1,"A simulation must be ran before infoplots can be generated.")
                    end

                elseif usrCmd == "vortplot"
                    tt = plotTime(subCmd)
                    makeVortPlots2D(tt)

                end
            end
        elseif usrCmd != "" # return error message
            errorHandler(usrCmd)
            continue
        end

    end
    #return commandDict
end
