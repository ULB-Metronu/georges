*******
Manzoni
*******

Tracking
########
* Python3
* numba
* Elements

Integrators
###########

* MadX
* Transport
* Mad8

Observers
#########

Different observers are implemented in `Manzoni`. The syntax is as follow ::

    observer = Observer(element=['eleA', 'eleB'])

The results are stored in a `pandas.DataFrame` and are available using::

    result = observer.to_df()

In the sections below, we detail each observer implemented in `Manzoni`

MeanObserver
************

This observer store the mean of each quantity in the beam
(:math:`\bar{x}`, :math:`\bar{px}`, :math:`\bar{y}`, :math:`\bar{py}`, :math:`\bar{dpp}`)

::

    observer = MeanObserver(element=['eleA', 'eleB'])


StdObserver
***********

This observer store the standard deviation of each quantity in the beam
(:math:`\sigma_x`, :math:`\sigma_{px}`, :math:`\sigma_y`, :math:`\sigma_{p}y`, :math:`\sigma_{dpp}`)

::

    observer = StdObserver(element=['eleA', 'eleB'])

BeamObserver
************

This observer store the entire beam distribution at the exit of an element.
If the flag input `with_input_beams` is set to `True`, the beam at the entrance of the
element is also stored.

::

    observer = BeamObserver(element=['eleA', 'eleB'], with_input_beams=True)

LossesObserver
**************

This observer store the number of particles, the losses and the transmission of an element.

::

    observer = LossesObserver(element=['eleA', 'eleB'])

SymmetryObserver
****************

This observer compute the symmetry of the beam at the entrance
and at the exit of an element. The symmetry is given by the following relation:

.. math::

     sym = \left| \frac{(\sigma_x - \sigma_y)}{(\sigma_x + \sigma_y)} \right|

::

    observer = SymmetryObserver(element=['eleA', 'eleB'])


IbaBPMObserver
**************

This observer is used to validate the model with the IBA' Beam Profile Monitor.
A Gaussian fit is performed directly on the beam position data. The element to store the
data must be a `Marker`.

::

    observer = IbaBpmObserver(element=['eleA', 'eleB'])

UserObserver
************

It is also possible to define his own observer. First, you must create a class that inherits
from the main Observer. The method `__call__(self, element, b1, b2)` receives the element
where to store data and the beam at the entry and exit of an element. The example below
illustrates the implementation of a mean observer.

::

    class MyMeanObserver(georges.manzoni.Observer):
        def __init__(self, elements = None):
            super().__init__(elements)
            self.headers = ('NAME',
                            'AT_ENTRY',
                            'AT_CENTER',
                            'AT_EXIT',
                            'BEAM_IN_X',
                            'BEAM_OUT_X',
                            'BEAM_IN_Y',
                            'BEAM_OUT_Y',
                            'BEAM_IN_XP',
                            'BEAM_OUT_XP',
                            'BEAM_IN_YP',
                            'BEAM_OUT_YP',
                            'BEAM_IN_DPP',
                            'BEAM_OUT_DPP',
                            )

        def __call__(self, element, b1, b2):
            if super().__call__(element, b1, b2):
                self.data.append((element.NAME,
                                  element.AT_ENTRY,
                                  element.AT_CENTER,
                                  element.AT_EXIT,
                                  b1[:, 0].mean(),
                                  b2[:, 0].mean(),
                                  b1[:, 2].mean(),
                                  b2[:, 2].mean(),
                                  b1[:, 1].mean(),
                                  b2[:, 1].mean(),
                                  b1[:, 3].mean(),
                                  b2[:, 3].mean(),
                                  b1[:, 4].mean(),
                                  b2[:, 4].mean(),
                                  ))

Now this observer can be used by `Manzoni`::

    beam_Myobserver = mi.track(beam=beam, observers=MyMeanObserver())


Validation
##########

The Manzoni tracking code can be validated against `MAD-X` using the converters described in
`georges-core`. The example here below illustrates how to define a line, then use `MAD-X`
and finally compare the results with `Manzoni`.

.. jupyter-execute::
   :hide-output:

    import matplotlib.pyplot as plt
    import georges
    from georges.manzoni import Input
    from georges.manzoni.beam import MadXBeam
    from georges.manzoni import observers
    from georges import vis
    _ureg = georges.ureg

We define the kinematics and the beam properties for the line
*************************************************************

.. jupyter-execute::
   :hide-output:

    kin = georges.Kinematics(230 *_ureg.MeV,
                         particle=georges.particles.Proton,
                         kinetic=True)
    betax = 3.81481846*_ureg.m
    betay = 2.30182336*_ureg.m

    alfax = -1
    alfay = 0.75

Let's define a line with drifts, quadrupoles and dipoles
********************************************************

.. jupyter-execute::
   :hide-output:

    d1 = georges.Element.Drift(NAME="D1",
                           L=0.3* _ureg.m,
                           APERTYPE="RECTANGULAR",
                           APERTURE=[10*_ureg.cm, 5*_ureg.cm])

    qf = georges.Element.Quadrupole(NAME="Q1",
                                    L=0.3*_ureg.m,
                                    K1=2*_ureg.m**-2,
                                    APERTYPE="RECTANGULAR",
                                    APERTURE=[10*_ureg.cm, 5*_ureg.cm])

    d2 = georges.Element.Drift(NAME="D2",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[10*_ureg.cm, 10*_ureg.cm])

    b1 = georges.Element.SBend(NAME="B1",
                               L=1*_ureg.m,
                               ANGLE=30*_ureg.degrees,
                               K1=0*_ureg.m**-2,
                               APERTYPE="CIRCULAR",
                               APERTURE=[10*_ureg.cm, 10*_ureg.cm])

    d3 = georges.Element.Drift(NAME="D3",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[10*_ureg.cm, 10*_ureg.cm])

    qd = georges.Element.Quadrupole(NAME="Q2",
                                    L=0.3*_ureg.m,
                                    K1=-2*_ureg.m**-2,
                                    APERTYPE="RECTANGULAR",
                                    APERTURE=[10*_ureg.cm, 5*_ureg.cm])

    d4 = georges.Element.Drift(NAME="D4",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[10*_ureg.cm, 10*_ureg.cm])

    b2 = georges.Element.SBend(NAME="B2",
                               L=1*_ureg.m,
                               ANGLE=-30*_ureg.degrees,
                               K1=0*_ureg.m**-2,
                               APERTYPE="RECTANGULAR",
                               APERTURE=[10*_ureg.cm, 5*_ureg.cm])

    d5 = georges.Element.Drift(NAME="D5",
                               L=0.3*_ureg.m,
                               APERTYPE="CIRCULAR",
                               APERTURE=[10*_ureg.cm, 10*_ureg.cm])

    sequence = georges.sequence.PlacementSequence(name="fodo")

    sequence.place(d1,at_entry=0*_ureg.m)
    sequence.place_after_last(qf)
    sequence.place_after_last(d2)
    sequence.place_after_last(b1)
    sequence.place_after_last(d3)
    sequence.place_after_last(qd)
    sequence.place_after_last(d4)
    sequence.place_after_last(b2)
    sequence.place_after_last(d5)

    sequence.metadata.kinematics = kin


.. jupyter-execute::
   :hide-output:

    mad_input = georges.madx.MadX(sequence=sequence);
    tfs_data = mad_input.twiss(sequence='fodo',
                               betx=betax.m_as('m'),
                               bety=betay.m_as('m'),
                               alfx=alfax,
                               alfy=alfay);

.. jupyter-execute::
   :hide-output:

    beam = MadXBeam(kinematics=kin,
            distribution=georges.Distribution.from_twiss_parameters(n=100000,
                                           betax=betax,
                                           betay=betay,
                                           alphax=alfax,
                                           alphay=alfay).distribution[["X", "PX", "Y", "PY", "DPP"]].values
           )

.. todo::

    Correct when georges-core distribution will be updated

.. jupyter-execute::
   :hide-output:

    mi = Input.from_sequence(sequence=sequence)
    beam_observer_tw = mi.track(beam=beam, observers=observers.TwissObserver())

Compare results between MAD-X and Manzoni
*****************************************

.. tabs::

   .. tab:: Matplotlib

      .. tabs::

         .. tab:: Beta

            .. jupyter-execute::

                fig = plt.figure(figsize=(10,4))
                ax = fig.add_subplot(111)
                manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
                manzoni_plot.plot_cartouche(sequence.df)
                manzoni_plot.twiss(beam_observer_tw, tfs_data=tfs_data)
                ax.legend(loc='upper left')

         .. tab:: Alpha

            .. jupyter-execute::

                fig = plt.figure(figsize=(10,4))
                ax = fig.add_subplot(111)
                manzoni_plot = vis.ManzoniMatplotlibArtist(ax=ax)
                manzoni_plot.plot_cartouche(sequence.df)
                manzoni_plot.twiss(beam_observer_tw, with_beta = False, with_alpha = True, tfs_data=tfs_data)
                ax.legend(loc='upper left')

         .. tab:: Dispersion

            .. todo::

                To implement

   .. tab:: Plotly

      .. tabs::

         .. tab:: Beta

            .. jupyter-execute::

                manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
                manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
                manzoni_plot.fig['layout']['legend'] =dict(
                    yanchor="top",
                    y=0.99,
                    xanchor="left",
                    x=0.01
                )
                manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.12)
                manzoni_plot.twiss(beam_observer_tw, with_beta=True, tfs_data=tfs_data)
                manzoni_plot.fig['data'][0]['showlegend'] = True
                manzoni_plot.fig['data'][1]['showlegend'] = True
                manzoni_plot.fig['data'][1]['showlegend'] = True
                manzoni_plot.render()

         .. tab:: Alpha

            .. jupyter-execute::

                manzoni_plot = vis.ManzoniPlotlyArtist(width=600, height=400)
                manzoni_plot.fig["layout"]["margin"] = dict(l=0, r=0, b=0)
                manzoni_plot.fig['layout']['legend'] =dict(
                    yanchor="top",
                    y=0.99,
                    xanchor="left",
                    x=0.01
                )
                manzoni_plot.plot_cartouche(sequence.df, unsplit_bends=False, vertical_position=1.12)
                manzoni_plot.twiss(beam_observer_tw, with_beta=False, with_alpha=True, tfs_data=tfs_data)
                manzoni_plot.fig['data'][0]['showlegend'] = True
                manzoni_plot.fig['data'][1]['showlegend'] = True
                manzoni_plot.fig['data'][1]['showlegend'] = True
                manzoni_plot.render()

         .. tab:: Dispersion

            .. todo::

                To implement
