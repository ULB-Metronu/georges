import matplotlib.pyplot as plt
import georges
PATH = "/Users/chernals/Dropbox/IBA/Work/IBA-Optics/beamlines"

brho = georges.physics.momentum_to_brho(georges.physics.energy_to_momentum(177.5))
context = {
    'PARTICLE': 'PROTON',
    'PC': georges.physics.energy_to_momentum(177.5),
    'BETAX': 0.0846155,
    'BETAY': 0.0846155,
    'ALPHAX': 0.0,
    'ALPHAY': 0.0,
    'DELTAP': 0.0,
    'DPP': 0.5e-2,
    'IQ1E': 2.741 / 10 / brho / 0.0325 * 0.297 / 0.295,
    'IQ2E': -5.745 / 10 / brho / 0.0325 * 0.297 / 0.295,
    'IQ3E': 4.079 / 10 / brho / 0.0325 * 0.297 / 0.295,
    'IQ47E': -3.062 / 10 / brho / 0.0325,
    'IQ56E': 2.66 / 10 / brho / 0.0325,
    'IQ8E': 3.5484 / 10 / brho / 0.0325,
    'IQ9E': -4.0276 / 10 / brho / 0.0325,
    'IQ10E': 4.0352 / 10 / brho / 0.0325,
    'MOMENTUM_SLITS_OPENING': 1.015,
    'DIVERGENCE_SLITS_OPENING_X': 1.015,
    'DIVERGENCE_SLITS_OPENING_Y': 1.015,
    'N_TRACKING': 5000,
    'EMITX': 14.3e-6,
    'EMITY': 14.3e-6,
}
m = georges.madx.Madx(madx='/usr/local/bin/madx', path=PATH, context=context)
bl = georges.Beamline('ess', path=PATH, prefix='generic', elements='elements')
bl_twiss = georges.madx.twiss(line=bl, madx=m)
b = georges.Beam(energy=200).from_5d_multigaussian_distribution(5,
                                                                XRMS=0.01,
                                                                PXRMS=0.001,
                                                                YRMS=0.01,
                                                                PYRMS=0.001,
                                                                DPPRMS=0.00)
from georges.plotting import *

with plt.style.context('word'):
    fig = plt.figure(1)

    ax1 = fig.add_subplot(111)
    prepare(ax1, bl)
    aperture(ax1, bl)
    # twiss(ax1, bl_twiss)
    # track(ax1, bl, b)
bl_track = georges.madx.track(line=bl, madx=m, beam=b)