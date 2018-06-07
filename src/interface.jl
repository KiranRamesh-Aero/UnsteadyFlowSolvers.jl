#using Tk: TopLevel, pack_stop_propagate, frame, pack, Label, Radio, Button

# function callback3(path)
#     if get_value(cb2) == "Constant"
#         l1 = Label(fInput, "Value")
#         e1 = Entry(fInput)
#     end
#     map(u -> pack(u, anchor="w"), (l1, e1))    
# end


# function callpack2(path)
#     probtype = get_value(rb1)
#     if probtype == "Prescribed kinematics"
#         l1 = Label(fInput, "Surface definition")
#         l2 = Label(fInput, "Geometry")
#         geomTypes = ["Flat plate", "NACA", "SD7003", "Read from file"]
#         cb1 = Combobox(fInput, geomTypes)
#         l3 = Label(fInput, "Pitch axis location (0-1)")
#         e1 = Entry(fInput)
#         l4 = Label(fInput, "Kinematics")
#         kinemTypes = ["Constant", "Sine", "Cosine", "Eldredge ramp", "Eldredge ramp-hold-return"]
#         l5 = Label(fInput, "Pitch")
#         cb2 = Combobox(fInput, kinemTypes)
#         bind(cb2, "command", callback3)
#         l6 = Label(fInput, "Plunge/Heave")
#         cb3 = Combobox(fInput, kinemTypes)
#         bind(cb3, "command", callback4)

#         map(u -> pack(u, anchor="w"), (l1, l2, cb1, l3, e1, l4, l5, cb2, l6, cb3))
#     end
# end

# function callback1(path)
#     simtype = get_value(rb)
#     if simtype == "2D"
#         l1 = Label(fInput, "Problem type")
#         rb1 = Radio(fInput, ["Prescribed kinematics"])
                
#         bind(rb1, "command", callback2)
        
#         map(u1 -> pack(u1, anchor="w"), (l1, rb1))
#     end
# end


w = Toplevel("UNSflow - unsteady flow solver", 600, 600)
pack_stop_propagate(w)
fInput = Frame(w)
pack(fInput, expand=true, fill="both")

fD = Frame(fInput, padding = [3, 3, 2, 2])
grid(fD, 1, 1)

fProb = Frame(fInput, padding = [3, 3, 2, 2])
grid(fProb, 2, 1)

fSurf = Frame(fInput, padding = [3, 3, 2, 2])
grid(fSurf, 3, 1)

fKinem = Frame(fInput, padding = [3, 3, 2, 2])
grid(fKinem, 4, 1)

l1 = Label(fD, "Simulation type")
rb1 = Radio(fD, ["2D", "3D"])
b1 = Button(fD, "ok")
grid(l1, 1, 1)
grid(rb1, 2, 1)
grid(b1, 3, 1)

bind(b1, "command", path -> giveProbs(get_value(rb1)))

function giveProbs(simtype::String)
    if simtype == "2D"
        l2 = Label(fProb, "Problem type")
        rb2 = Radio(fProb, ["Prescribed kinematics"])
        b2 = Button(fProb, "ok")
        grid(l2, 1, 1)
        grid(rb2, 2, 1)
        grid(b2, 3, 1)
        Tk.update()

        bind(b2, "command", path -> giveSurf(get_value(rb2)))
    end
end

function giveSurf(probtype::String)

    if probtype == "Prescribed kinematics"
        l3 = Label(fSurf, "Surface definition")
        l4 = Label(fSurf, "Geometry")
        geomTypes = ["Flat plate", "NACA", "SD7003", "Read from file"]
        cb1 = Combobox(fSurf, geomTypes)
        l5 = Label(fSurf, "Pitch axis location (0-1)")
        e1 = Entry(fSurf)
        l6 = Label(fSurf, "Kinematics")
        kinemTypes = ["Constant", "Sine", "Cosine", "Eldredge ramp", "Eldredge ramp-hold-return"]
        l7 = Label(fSurf, "Pitch")
        cb2 = Combobox(fSurf, kinemTypes)
        l8 = Label(fSurf, "Plunge/Heave")
        cb3 = Combobox(fSurf, kinemTypes)
        b3 = Button(fSurf, "ok")
        grid(l3, 1, 1)
        grid(l4, 2, 1)
        grid(cb1, 3, 1)
        grid(l5, 4, 1)
        grid(e1, 5, 1)
        grid(l6, 6, 1)
        grid(l7, 7, 1)
        grid(cb2, 8, 1)
        grid(l8, 9, 1)
        grid(cb3, 10, 1)
        grid(b3, 11, 1)
        Tk.update()

        bind(b3, "command", path -> giveKinem(get_value(cb2), get_value(cb3)))
    end
end

function giveKinem(pitchtype::String, plungetype::String)
    if pitchtype == "Constant"
        l1 = Label(fKinem, "Pitch $pitchtype properties") 
        l2 = Label(fKinem, "Value")
        e1 = Entry(fKinem)
        grid(l1, 1, 1)
        grid(l2, 2, 1)
        grid(e1, 3, 1)

        Tk.update()
    end
    if pitchtype == "Sine"
        l1 = Label(fKinem, "Pitch $pitchtype properties") 
        l2 = Label(fKinem, "Mean")
        e1 = Entry(fKinem)
        l3 = Label(fKinem, "Amplitude (deg)")
        e2 = Entry(fKinem)
        l4 = Label(fKinem, "Reduced frequency (omega*c/(2*u))")
        e3 = Entry(fKinem)
        l5 = Label(fKinem, "Phase angle (deg)")
        e4 = Entry(fKinem)
        grid(l1, 1, 1)
        grid(l2, 2, 1)
        grid(e1, 3, 1)
        grid(l3, 4, 1)
        grid(e2, 5, 1)
        grid(l4, 6, 1)
        grid(e3, 7, 1)
        grid(l5, 8, 1)
        grid(e4, 9, 1)

        Tk.update()
    end
end

# bind(rb, "command", callback)

# map(u -> pack(u, anchor="fD"), (l, rb))

# fSurf = Frame(fInput)


# surfTypes = [""]


# #surfTypes = ["Prescribed", "Constrained-2DOF"]
# 


# l1 = Label(fSurf, "Surface properties:")
# cb2 = Combobox(fSurf, surfTypes)
# b = Button(fSurf, "ok")
# map(u -> pack(u, anchor="w"), (l1, l2, cb1, l3, e1, cb2, b))

#

# fKinem = Frame(w)
# pack(fKinem, expand=true, fill="both")

# l = Label(fKinem, "Kinematics definition:")
# l = Label(fKinem, "Pitch")

# cb = Combobox(fKinem, kinemTypes)
# b = Button(fKinem, "ok")
# map(u -> pack(u, anchor="w"), (l, cb, b))



