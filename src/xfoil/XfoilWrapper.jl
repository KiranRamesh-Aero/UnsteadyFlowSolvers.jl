#=
XfoilWrapper.jl
    Written by Matthew White
    4/2/2019

    Input:  config file by name of Config.txt
    Output: Cn,Cm and polar file
=#
function xfoilWrapper(cnfg::String = "bin/config.txt")
    # import config file
    cnfg = DelimitedFiles.readdlm(cnfg, String)
    # airfoil
    uppercase(cnfg[1,2]) == "NACA" ? airfoil = cnfg[1,3] : airfoil = cnfg[1,2]
    re = cnfg[2,2] # Reynolds
    mach = cnfg[3,2] # Mach #
    changePanels = uppercase(cnfg[4,2]) .== "TRUE" # Manually change panels (otherwise runs PANE)
        panels = cnfg[5,2] # Panels to change to
    changeIter = uppercase(cnfg[6,2]) .== "TRUE" # change the number of iterations
        iter = cnfg[7,2] # iteration count to change to
    # AoA range
    aStart = parse(Float32, cnfg[8,2])
    aEnd = parse(Float32, cnfg[9,2])
    aStep = parse(Float32, cnfg[10,2])
    # XFOIL path
    xPath = cnfg[11,2]
    # Polar file path
    polarPath = cnfg[12,2]

    #   --- Open XFOIL Pipe ---
    p=open(pipeline(`$xPath`),"r+")
    run(`rm -rf $polarPath`) # remove any old polar file
    # Write to pipe
    length(airfoil) <= 5 ? write(p,"naca $airfoil\n") : write(p,"load $airfoil\n")
    changePanels ? write(p,"ppar $panels\n\n") : write(p,"pane\n")
    write(p,"oper\n")
    write(p,"visc $re\n")
    if changeIter write(p,"iter $iter\n") end
    write(p,"mach $mach\n")
    write(p,"pacc\n")
    write(p,"$polarPath\n\n")
    write(p,"aseq $aStart $aEnd $aStep\n")
    write(p,"\nquit")

    close(p)

    while !isfile(polarPath) end # wait unitl polar file has been created
    while filesize(polarPath) < 1000 end # wait until at least one value has been written

    #   --- Reading the Polar File ---
    prevSize = 0
    sec = 0
    n = 10
    while true # wait unit the polar file is done changing
        global polar = DelimitedFiles.readdlm(polarPath,Float64,skipstart = 12)
        if polar[end,1] > aEnd - aStep
            println("Done writing polar file")
            break #File is fully loaded
        else #check every n seconds for stable file size
            sleep(1)
            sec += 1
            if prevSize == filesize(polarPath) && sec == n
                println("Method could not converge for the final alphas")
                println("   Maximum angle of attack at $(polar[end,1]) degrees")
                break # after n seconds if filesize is constant then use current polar
            elseif sec >= n
                sec = 0 # reset sec every n seconds
            end
            prevSize = filesize(polarPath)
        end
    end

    # polar: alpha | Cl | CD | CDp | Cm | Top_Xtr | Bot_Xtr
    # Cn = Cl*cos(a) + Cd*sin(a)
    Cn = polar[:,2].*cosd.(polar[:,1]) + polar[:,3].*sind.(polar[:,1])
    Cm = polar[:,5]

    Cn,Cm

end
