*******************************************
Georges’s visualization and plotting module
*******************************************

This module uses the `georges-core's plotting module <https://ulb-metronu.github.io/georges-core/modules/plotting.html>`_ . It can uses the *Matplotlib* ot the *Plotly* library as backend.
After defining a *ManzoniArtist()*, the beamline can be displayed with different options:

+----------------+------------------+----------------+
| argument       |  type            |  default value |
+================+==================+================+
| beamline       | pandas.DataFrame | None           |
+----------------+------------------+----------------+
| print_label    | bool             | False          |
+----------------+------------------+----------------+
| with_aperture  | bool             | False          |
+----------------+------------------+----------------+
| plane          | str              | None           |
+----------------+------------------+----------------+

Example::

    import matplotlib.pyplot as plt
    plt.rc('text', usetex=False)
    fig = plt.figure(figsize=(20,8))
    ax = fig.add_subplot(111)
    manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
    manzoni_plot.plot_cartouche(sequence.df) # Preparation of the plot
    manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=False, plane='X')

If an observer has been used in the simulation, the results can be added to the previous figure depending the
instance of the observer:

.. jupyter-execute::
   :hide-output:
   :hide-code:

    import georges
    from georges.manzoni import Input
    from georges.manzoni.beam import MadXBeam
    from georges.manzoni.integrators import *
    from georges.manzoni import observers
    from georges import vis
    _ureg = georges.ureg
    aper1 = 10
    aper2 = 5
    d1 = georges.Element.Drift(NAME="D1",
                               L=0.3* _ureg.m,
                               APERTYPE="RECTANGULAR",
                               APERTURE=[aper1*_ureg.cm, aper2*_ureg.cm])

    qf = georges.Element.Quadrupole(NAME="Q1",
                                    L=0.3*_ureg.m,
                                    K1=2*_ureg.m**-2,
                                    APERTYPE="RECTANGULAR",
                                    APERTURE=[aper1*_ureg.cm, aper2*_ureg.cm])

    d2 = georges.Element.Drift(NAME="D2",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[aper1*_ureg.cm, aper1*_ureg.cm])

    b1 = georges.Element.SBend(NAME="B1",
                               L=1*_ureg.m,
                               ANGLE=30*_ureg.degrees,
                               K1=0*_ureg.m**-2,
                               APERTYPE="CIRCULAR",
                               APERTURE=[aper1*_ureg.cm, aper1*_ureg.cm])

    d3 = georges.Element.Drift(NAME="D3",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[aper1*_ureg.cm, aper1*_ureg.cm])

    qd = georges.Element.Quadrupole(NAME="Q2",
                                    L=0.3*_ureg.m,
                                    K1=-2*_ureg.m**-2,
                                    APERTYPE="RECTANGULAR",
                                    APERTURE=[aper1*_ureg.cm, aper2*_ureg.cm])

    d4 = georges.Element.Drift(NAME="D4",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[aper1*_ureg.cm, aper1*_ureg.cm])

    b2 = georges.Element.SBend(NAME="B2",
                               L=1*_ureg.m,
                               ANGLE=-30*_ureg.degrees,
                               K1=0*_ureg.m**-2,
                               APERTYPE="RECTANGULAR",
                               APERTURE=[aper1*_ureg.cm, aper2*_ureg.cm])

    d5 = georges.Element.Drift(NAME="D5",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[aper1*_ureg.cm, aper1*_ureg.cm])

    sequence = georges.PlacementSequence(name="Sequence")

    sequence.place(d1,at_entry=0*_ureg.m)
    sequence.place_after_last(qf)
    sequence.place_after_last(d2)
    sequence.place_after_last(b1)
    sequence.place_after_last(d3)
    sequence.place_after_last(qd)
    sequence.place_after_last(d4)
    sequence.place_after_last(b2)
    sequence.place_after_last(d5);

    kin = georges.Kinematics(230 * _ureg.MeV,
                             particle=georges.particles.Proton,
                             kinetic=True)

    # Add kinematics to the sequence
    sequence.metadata.kinematics=kin

    beam = MadXBeam(kinematics=kin,
                distribution=georges.generate_from_5d_sigma_matrix(n=10000,
                                                               s11=0.001,
                                                               s22=0.001,
                                                               s33=0.005,
                                                               s44=0.005)
           )
    mi = Input.from_sequence(sequence=sequence)

    beam_observer_std = mi.track(beam=beam, observers=observers.SigmaObserver())
    beam_observer_mean = mi.track(beam=beam, observers=observers.MeanObserver())
    beam_observer_beam = mi.track(beam=beam, observers=observers.BeamObserver(with_input_beams=True))
    beam_observer_losses = mi.track(beam=beam, observers=observers.LossesObserver())
    beam_observer_tw = mi.track(beam=beam, observers=observers.TwissObserver())

.. jupyter-execute::
   :hide-output:
   :hide-code:

    import matplotlib.pyplot as plt
    plt.rc('text', usetex=False)

Mean Observer
#############

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df)
        manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=True, plane='X')
        manzoni_plot.tracking(beam_observer_mean, plane='X')

   .. tab:: Plotly

      .. jupyter-execute::

        manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
        manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
        manzoni_plot.fig['layout']['legend'] = dict(yanchor="top",
                                                    y=0.99,
                                                    xanchor="left",
                                                    x=0.01)
        manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
        manzoni_plot.tracking(beam_observer_mean, plane='X')
        manzoni_plot.fig['data'][0]['showlegend'] = True
        manzoni_plot.render()

Std Observer
############

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df)
        manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=True, plane='X')
        manzoni_plot.tracking(beam_observer_std, plane='both')

   .. tab:: Plotly

      .. jupyter-execute::

        manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
        manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
        manzoni_plot.fig['layout']['legend'] = dict(yanchor="top",
                                                    y=0.99,
                                                    xanchor="left",
                                                    x=0.01)
        manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
        manzoni_plot.tracking(beam_observer_std, plane='both')
        manzoni_plot.fig['data'][0]['showlegend'] = True
        manzoni_plot.fig['data'][1]['showlegend'] = True
        manzoni_plot.render()

Beam Observer
#############

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df)
        manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=True, plane='X')
        manzoni_plot.tracking(beam_observer_beam, fill_between=False, plane='X', mean=False, std=False, halo=True)

   .. tab:: Plotly

      .. jupyter-execute::

        manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
        manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
        manzoni_plot.fig['layout']['legend'] = dict(yanchor="top",
                                                    y=0.99,
                                                    xanchor="left",
                                                    x=0.01)
        manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
        manzoni_plot.tracking(beam_observer_beam, fill_between=False, plane='X', mean=False, std=False, halo=True)
        manzoni_plot.fig['data'][0]['showlegend'] = True
        manzoni_plot.render()


Phase Space Observer
####################

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df)
        manzoni_plot.plot_beamline(sequence.df, print_label=False, with_aperture=True, plane='X')
        manzoni_plot.phase_space(observer=beam_observer_beam,
                                 element='Q1',
                                 dim=['X', 'PX'],
                                 location='OUT',
                                 nbins=[51, 51])

   .. tab:: Plotly

      .. warning::

        This method is not yet implemented.

Losses Observer
###############

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df) # Preparation of the plot
        manzoni_plot.losses(beam_observer_losses, log_scale=True)

   .. tab:: Plotly

      .. jupyter-execute::

        manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
        manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
        manzoni_plot.fig['layout']['legend'] = dict(yanchor="top",
                                                    y=0.99,
                                                    xanchor="left",
                                                    x=0.01)
        manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
        manzoni_plot.losses(beam_observer_losses, log_scale=True)
        manzoni_plot.render()


Twiss Observer
##############

.. tabs::

   .. tab:: Matplotlib

      .. jupyter-execute::

        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot(111)
        manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
        manzoni_plot.plot_cartouche(sequence.df)
        manzoni_plot.twiss(beam_observer_tw, with_beta=True, with_alpha=False, with_dispersion=False)

   .. tab:: Plotly

      .. jupyter-execute::

        manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
        manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
        manzoni_plot.fig['layout']['legend'] = dict(yanchor="top",
                                                    y=0.99,
                                                    xanchor="left",
                                                    x=0.01)
        manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.15)
        manzoni_plot.twiss(beam_observer_tw, with_beta=True, with_alpha=False, with_dispersion=False)
        manzoni_plot.fig['data'][0]['showlegend'] = True
        manzoni_plot.render()
