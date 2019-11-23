#=
s1Cmds.jl
    Written by Matthew White
    11/14/2019
=#
function loadGeo(geo = "")
    # Assign airfoil file data
    exts = ["dat","txt","csv"] # valid file extensions
    while true
        if geo == "" # if no secondary command given
            geo = ask("Enter an airfoil geometry file: ")
            if lowercase(geo) == "exit" || lowercase(geo) == ""
                break
            end
        end
        if lowercase(geo) == "flatplate"
            return "FlatPlate"
        end
        if geo[1] == '/'
            geo = geo[2:end]
        end
        if occursin(".",geo) # check for approproite file extension
            usrExt = split(geo,".")
            usrExt = usrExt[end]
            if sum(usrExt .== exts) == 0 # if not a valid extension
                errorHandler(geo,"'$usrExt' is not a valid file extension.")
                geo = ""
                continue
            end
        else
            errorHandler(geo,"File must have an extension.")
            geo = ""
            continue
        end
        if occursin("/",geo) # Check if file exsists
            path = split(geo,"/")
            fn = path[end]
            path =join(path[1:end-1],"/")
            files = cd(readdir, path)
        else
            files = readdir()
            fn = geo
        end
        if sum(files .== fn) == 0 # file does not exist in directory
            errorHandler(fn,"'$fn' does not exist.")
            geo = ""
            continue
        end
        return geo
    end
end

function modef(motype = "")
    # Create alphadef/hdef/udef motion vectors

    modefs = ["const","lin","eld","sin","cos"]#,"gust"]
    outVar = []

    while true # Motion type
        if motype == "" # if no secondary command given
            motype = lowercase(ask("Motion Type: "))
            if motype == "exit" || motype == ""
                return nothing
            end
        end
        if sum(motype .== modefs) > 0 # motion type matches
            outVar = [outVar;motype]
            break # next question
        else
            # Error message
            prompt = join(modefs[1:end-1], ", ")
            prompt = "Motion options are $prompt, and $(modefs[end])"
            errorHandler(motype,prompt)
        end
        motype = ""
    end

    if motype == "const"
        amp = svPrompt("f","Magnitude: ")
        if amp == nothing ; return nothing end
        outVar = [outVar;amp]
    elseif motype == "lin"
        sm = svPrompt("f","Start Magnitude: ")
        if sm == nothing ; return nothing end
        fm = svPrompt("f","Final Magnitude: ")
        if fm == nothing ; return nothing end
        l = svPrompt("f","Line Length: ")
        if l == nothing ; return nothing end
        ynAsk("Change start time?") ? tstart = svPrompt("f","Start Time: ") : tstart = 1
        if tstart == nothing ; return nothing end
        outVar = [outVar;tstart;sm;fm;l]
    elseif motype == "eld"
        amp = svPrompt("f","Amplitude: ")
        if amp == nothing ; return nothing end
        K = svPrompt("f","K: ")
        if K == nothing ; return nothing end
        ynAsk("Change a?") ? a = svPrompt("f","a: ") : a = (pi^2 * K*180 )/(2*amp*pi *(1-0.1))
        if a == nothing ; return nothing end
        ynAsk("Change start time?") ? tstart = svPrompt("f","Start Time: ") : tstart = 1
        if tstart == nothing ; return nothing end
        outVar = [outVar;amp;K;a;tstart]
    elseif motype == "sin" || motype == "cos"
        mean = svPrompt("f","Mean: ")
        if mean == nothing ; return nothing end
        amp = svPrompt("f","Amplitude: ")
        if amp == nothing ; return nothing end
        K = svPrompt("f","K: ")
        if K == nothing ; return nothing end
        phi = svPrompt("f","Phi: ")
        if phi == nothing ; return nothing end
        outVar = [outVar;mean;amp;K;phi]
    elseif motype == "gust"
        amp = svPrompt("f","Magnitude: ")
        if amp == nothing ; return nothing end
        ynAsk("Change start time?") ? tstart = svPrompt("f","Start Time: ") : tstart = 1
        if tstart == nothing ; return nothing end
        tgust = svPrompt("f","Gust Start Time: ")
        if tgust == nothing ; return nothing end
        outVar = [outVar;amp;tstart;tgust]
    end
    return outVar
end

function setCamberCalc(meth)
    methods = ["radial","linear"]
    eStr = join(methods[1:end-1], ",")
    eStr = "Accepts either $eStr, or $(methods[end])"
    while true
        if meth == ""
            prompt = "Camber calculation method: "
            meth = lowercase(ask(prompt))
        end
        if sum(meth .== methods) == 0 # not a method
            errorHandler(meth,eStr)
        else
            return meth
        end
        meth = ""
    end
end

function alphaCorr(alphavec)
    # Turn alpha into radians
    alphavec[2] = alphavec[2]*pi/180

    if alphavec[1] == "const"
        alphadef = ConstDef(alphavec[2])
    elseif alphavec[1] == "lin"
        alphavec[3:4] = alphavec[3:4].*pi./180
        alphadef = LinearDef(alphavec[2],alphavec[3],alphavec[4],alphavec[5])
    elseif alphavec[1] == "eld"
        alphadef = EldRampReturntstartDef(alphavec[2],alphavec[3],alphavec[4],alphavec[5])
    elseif alphavec[1] == "sin"
        alphavec[3] = alphavec[3]*pi/180
        alphavec[5] = alphavec[5]*pi/180
        alphadef = SinDef(alphavec[2],alphavec[3],alphavec[4],alphavec[5])
    elseif alphavec[1] == "cos"
        alphavec[3] = alphavec[3]*pi/180
        alphavec[5] = alphavec[5]*pi/180
        alphadef = CosDef(alphavec[2],alphavec[3],alphavec[4],alphavec[5])
    elseif alphavec[1] == "gust"
        alphadef = StepGustDef(alphavec[2],alphavec[3],alphavec[4])
    end
    return alphadef
end

function gotoOper(cmds)
    # Cheack whether there are enough inputs to move to operation stage (2)
    o = false
    notdef = [k for (k,v) in cmds if v=="---"]
    if length(notdef) > 1 # incomplete intitilization
        prompt = join(notdef[1:end-1], ", ")
        prompt = "$prompt, and $(notdef[end]) have not been defined."
        errorHandler(o,prompt)
    elseif length(notdef) == 1
        prompt = "$(notdef[1]) has not been defined."
        errorHandler(o,prompt)
    else
        println(" The simulation will take $(Int(ceil(cmds["runtime"]/cmds["dt"]))) iterations.")
        o = true
    end
    return o
end
