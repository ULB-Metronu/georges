import matplotlib.pyplot as plt
from georges.manzoni import Input
from georges.manzoni.integrators import *
from georges.manzoni import observers
from georges import vis

_ureg = georges.ureg


def plot_mean_observer(sequence, beam_observer_mean):
    # Using Plotly
    manzoni_plot = vis.ManzoniPlotlyArtist(width=800, height=600)
    manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
    manzoni_plot.tracking(beam_observer_mean, plane='X')
    manzoni_plot.tracking(beam_observer_mean, plane='Y')

    # Using matplotlib
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df)
    manzoni_plot.plot_beamline(sequence.df, print_label=True, with_aperture=True)
    manzoni_plot.tracking(beam_observer_mean, plane='X')
    manzoni_plot.tracking(beam_observer_mean, plane='Y')


def plot_std_observer(sequence, beam_observer_std):
    # Using Plotly
    manzoni_plot = vis.ManzoniPlotlyArtist(width=800, height=600)
    manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
    manzoni_plot.tracking(beam_observer_std, plane='both', fill_between=False)

    # Using matplotlib
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df)
    manzoni_plot.plot_beamline(sequence.df, print_label=True, with_aperture=True, plane='both')
    manzoni_plot.tracking(beam_observer_std, plane='both', fill_between=False)


def plot_beam_observer(sequence, beam_observer_beam):
    # Using Plotly
    manzoni_plot = vis.ManzoniPlotlyArtist(width=800, height=600)
    manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
    manzoni_plot.tracking(beam_observer_beam, fill_between=False, plane='X', mean=False, std=False, halo=True)

    # Using matplotlib
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df)
    manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=True, plane='X')  # Preparation of the plot
    manzoni_plot.tracking(beam_observer_beam, fill_between=False, plane='X', mean=False, std=False, halo=True)


def plot_losses_observer(sequence, beam_observer_losses):
    # Using Plotly
    manzoni_plot = vis.ManzoniPlotlyArtist(width=800, height=600)
    manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
    manzoni_plot.losses(beam_observer_losses, log_scale=True)

    # # Using matplotlib
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df, print_label=True)  # Preparation of the plot
    manzoni_plot.losses(beam_observer_losses, log_scale=True)


def plot_phase_space(sequence, beam_observer_beam):
    # Using matplotlib
    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df)  # Preparation of the plot
    manzoni_plot.phase_space(observer=beam_observer_beam,
                             element='Q1',
                             dim=['X', 'PX'],
                             location='OUT',
                             nbins=[51, 51])


def test_plotting():
    pass
    # integrator = MadXIntegrator
    # aper1 = 10
    # aper2 = 15
    #
    # d1 = georges.Element.Drift(NAME="D1",
    #                            integrator=integrator,
    #                            L=0.3 * _ureg.m,
    #                            APERTYPE="RECTANGULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper2 * _ureg.cm])
    #
    # qf = georges.Element.Quadrupole(NAME="Q1",
    #                                 integrator=integrator,
    #                                 L=0.3 * _ureg.m,
    #                                 K1=2 * _ureg.m ** -2,
    #                                 APERTYPE="RECTANGULAR",
    #                                 APERTURE=[aper1 * _ureg.cm, aper2 * _ureg.cm])
    #
    # d2 = georges.Element.Drift(NAME="D2",
    #                            integrator=integrator,
    #                            L=0.3 * _ureg.m,
    #                            APERTYPE="CIRCULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper1 * _ureg.cm])
    #
    # b1 = georges.Element.SBend(NAME="B1",
    #                            integrator=integrator,
    #                            L=1 * _ureg.m,
    #                            ANGLE=30 * _ureg.degrees,
    #                            K1=0 * _ureg.m ** -2,
    #                            APERTYPE="CIRCULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper1 * _ureg.cm])
    #
    # d3 = georges.Element.Drift(NAME="D3",
    #                            integrator=integrator,
    #                            L=0.3 * _ureg.m,
    #                            APERTYPE="CIRCULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper1 * _ureg.cm])
    #
    # qd = georges.Element.Quadrupole(NAME="Q2",
    #                                 integrator=integrator,
    #                                 L=0.3 * _ureg.m,
    #                                 K1=-2 * _ureg.m ** -2,
    #                                 APERTYPE="RECTANGULAR",
    #                                 APERTURE=[aper1 * _ureg.cm, aper2 * _ureg.cm])
    #
    # d4 = georges.Element.Drift(NAME="D4",
    #                            integrator=integrator,
    #                            L=0.3 * _ureg.m,
    #                            APERTYPE="CIRCULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper1 * _ureg.cm])
    #
    # b2 = georges.Element.SBend(NAME="B2",
    #                            integrator=integrator,
    #                            L=1 * _ureg.m,
    #                            ANGLE=-30 * _ureg.degrees,
    #                            K1=0 * _ureg.m ** -2,
    #                            APERTYPE="RECTANGULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper2 * _ureg.cm])
    #
    # d5 = georges.Element.Drift(NAME="D5",
    #                            integrator=integrator,
    #                            L=0.3 * _ureg.m,
    #                            APERTYPE="CIRCULAR",
    #                            APERTURE=[aper1 * _ureg.cm, aper1 * _ureg.cm])
    #
    # sequence = georges.PlacementSequence(name="Sequence")
    #
    # sequence.place(d1, at_entry=0 * _ureg.m)
    # sequence.place_after_last(qf)
    # sequence.place_after_last(d2)
    # sequence.place_after_last(b1)
    # sequence.place_after_last(d3)
    # sequence.place_after_last(qd)
    # sequence.place_after_last(d4)
    # sequence.place_after_last(b2)
    # sequence.place_after_last(d5)
    #
    # kin = georges.Kinematics(230 * _ureg.MeV,
    #                          particle=georges.particles.Proton,
    #                          kinetic=True)
    #
    # # Add kinematics to the sequence
    # sequence.metadata.kinematics = kin
    # input_beam = georges.manzoni.Beam(kinematics=kin,
    #                                   distribution=georges.generate_from_5d_sigma_matrix(n=100000,
    #                                                                                      s11=0.0001,
    #                                                                                      s22=0.0001,
    #                                                                                      s33=0.0005)
    #                                   )
    # mi = georges.manzoni.Input.from_sequence(sequence=sequence)
    # mi.freeze()
    # beam_observer_std = mi.track(beam=input_beam, observers=observers.SigmaObserver())
    # beam_observer_mean = mi.track(beam=input_beam, observers=observers.MeanObserver())
    # beam_observer_beam = mi.track(beam=input_beam, observers=observers.BeamObserver(with_input_beams=True))
    # beam_observer_losses = mi.track(beam=input_beam, observers=observers.LossesObserver())
    # plot_mean_observer(sequence, beam_observer_mean)
    # plot_std_observer(sequence, beam_observer_std)
    # plot_beam_observer(sequence, beam_observer_beam)
    # plot_losses_observer(sequence, beam_observer_losses)
    # plot_phase_space(sequence, beam_observer_beam)
    #
    # assert True
