import matplotlib
matplotlib.rcParams["mathtext.fontset"] = "cm"
import schemdraw
import schemdraw.elements as elm
import users.danielsank.schemdraw.elements as pylelements


def fig_unmatched_resistor():
    d = schemdraw.Drawing()
    noise_source = d.add(pylelements.NoiseSource())
    noise_source_b = d.add(pylelements.NoiseSource().at(xy=noise_source.center, dx=6))
    d.add(elm.GroundSignal().at(noise_source.input))
    d.add(elm.Dot().at(noise_source.output).label("$S_{V_s}=4k_b T R$", "left"))
    d.add(elm.Resistor().label("$R_s$").at(noise_source.output).length(2).up())

    d.add(elm.GroundSignal().at(noise_source_b.input))
    d.add(elm.Dot().at(noise_source_b.output).label(r"$V_{s,rms}=\sqrt{4 k_b T R B}$", "left"))
    resistor_b = d.add(elm.Resistor().label("$R_s$").at(noise_source_b.output).length(2).up())
    wire = d.add(elm.Wire().at(resistor_b.end).to(resistor_b.end, dx=2))
    d.add(elm.Dot().label(r"$V_{l,rms} = V_{s,rms} \, \frac{R_l}{R + R_l}$", "right"))
    load = d.add(elm.Resistor().label("$R_l$", "bottom").at(wire.end).down().length(2))
    d.add(elm.Wire().at(load.end).to(load.end, dx=2).down())
    d.add(elm.GroundSignal())
    return d


def fig_resistor_with_tline():
    separation = 6
    d = schemdraw.Drawing()
    noise_source_a = d.add(pylelements.NoiseSource())
    noise_source_b = d.add(pylelements.NoiseSource().at(xy=noise_source_a.center, dx=separation))
    d.add(elm.GroundSignal().at(noise_source_a.input))
    d.add(elm.Dot().at(noise_source_a.output).label(r"$V_{s,rms}$", "left"))
    resistor_a = d.add(elm.Resistor().label("$R_s$").at(noise_source_a.output).length(2).up())
    d.add(elm.Dot().at(resistor_a.end).label("$V_{s,rms}/2$", ofst=(-1, 0)))
    tline_a = d.add(elm.Coax().label("$Z=R_s$", ofst=-0.55).at(resistor_a.end).right())
    #d.add(elm.ZLabel().label("$Z=R_s$").at((tline_a.start[0]+2, tline_a.start[1])))
    #d.add(elm.ZLabel(theta=180).label("$Z=R_s$", "top").at((tline_a.end[0]-2, tline_a.end[1])).flip())
    d.add(elm.CurrentLabel(ofst=0.5).at(tline_a).label("$V^+ = V_{s,rms} / 2$"))
    load_a = d.add(elm.Resistor().down().length(2).label("$R_s$", "top"))
    d.add(elm.Wire().at(load_a.end).to(load_a.end, dx=2).down())
    d.add(elm.GroundSignal())

    d.add(elm.GroundSignal().at(noise_source_b.input))
    resistor_b = d.add(elm.Resistor().up().at(noise_source_b.output).length(2))
    tline_b = d.add(elm.Coax().label("$Z=R_s$", ofst=-0.55).at(resistor_b.end).right())
    tline_c = d.add(elm.Coax().label("$Z=R_s$", ofst=-0.55).at(xy=tline_b.center, dx=separation-3).right())
    d.add(elm.CurrentLabel(ofst=0.5).at(tline_b).label("$V^+$"))
    d.add(elm.CurrentLabel(ofst=-0.5).at(tline_b).reverse().label("$V^-$", "bottom"))
    load_b = d.add(elm.RBox().at(tline_b.end).down().label("$Z_l$", "bottom"))
    d.add(elm.Wire().at(load_b.end).to(load_b.end, dy=-1))
    d.add(elm.GroundSignal())

    d.add(elm.CurrentLabel(ofst=0.5).at(tline_c).label(r"$V^+ = \sqrt{P_a R_s}$"))
    d.add(elm.CurrentLabel(ofst=-0.5).at(tline_c).reverse().label("$V^-$", "bottom"))
    load_c = d.add(elm.RBox().at(tline_c.end).down().label("$Z_l$", "bottom"))
    d.add(elm.Wire().at(load_c.end).to(load_c.end, dy=-1))
    d.add(elm.GroundSignal())
    return d


def fig_attenuator():
    d = schemdraw.Drawing()
    coax = d.add(elm.Coax())
    att = d.add(elm.RBox().at(coax.end).label("$G$", "bottom"))
    d.add(elm.CurrentLabel().at((coax.center[0], coax.center[1])).label("noise in"))
    d.add(elm.ZLabel().at((att.start[0]+1.2, att.start[1]-0.7)).label("noise added"))
    return d


def fig_worked_example():
    d = schemdraw.Drawing()
    noise_source = d.add(pylelements.NoiseSource())
    signal = d.add(elm.SourceSin().up().at(noise_source.output).length(2))
    d.add(elm.GroundSignal().at(noise_source.input))
    resistor_source = d.add(elm.Resistor().at(signal.end).length(2).label("$R_s$"))
    coax_a = d.add(elm.Coax().at(resistor_source.end).length(2).right())
    d.add(elm.CurrentLabel(ofst=0.5).at(coax_a).label("$S_{P_a}^{source}$"))
    d.add(elm.CurrentLabel(ofst=-0.5).at(coax_a).label("$P_a^{source}$", "bottom"))
    attenuator_a = d.add(elm.RBox().right().length(2).label(r"$G_1, \, T_1$"))
    d.add(elm.ZLabel().at((attenuator_a.start[0]+1.2, attenuator_a.start[1]-1)).label("$S_{P_a}^{att. 1}$"))
    d.add(elm.Coax().at(attenuator_a.end).length(2).right())
    attenuator_b = d.add(elm.RBox().right().length(2).label(r"$G_2, \, T_2$"))
    d.add(elm.ZLabel().at((attenuator_b.start[0]+1.2, attenuator_a.start[1]-1)).label("$S_{P_a}^{att. 2}$"))
    d.add(elm.Coax().at(attenuator_b.end).length(2).right())
    capacitor_drive = d.add(elm.Capacitor().label("$C_d$").right().length(1))
    wire = d.add(elm.Wire().at(capacitor_drive.end).to(capacitor_drive.end, dx=2))
    junction = d.add(elm.Josephson().at(wire.end).down())
    capacitor = d.add(elm.Capacitor().at(wire.start).down())
    wire_bottom = d.add(elm.Wire().at(junction.end).to(capacitor.end))
    d.add(elm.GroundSignal().at(wire_bottom.end, dx=1))

    return d
