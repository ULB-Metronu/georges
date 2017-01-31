def aperture_plot(ax, beamline, context, planes):
    global palette

    if planes == 'both':
        up = 1
        down = 0
    elif planes == 'X':
        up = 0
        down = 0
    elif planes == 'Y':
        up = 1
        down = 1
    else:
        up = 1
        down = 0

    def draw_quad(ax, e):
        # Poles
        ax.add_patch(
            matplotlib.patches.Rectangle(
                (e['S'] - e['L'] / 2.0, 1000*e['POLES'][0]),  # (x,y)
                e['L'],  # width
                100,  # height
                hatch='.', facecolor=palette['quad']
            )
        )
        ax.add_patch(
            matplotlib.patches.Rectangle(
                (e['S'] - e['L'] / 2.0, -1000*e['POLES'][0]),  # (x,y)
                e['L'],  # width
                -100,  # height
                hatch='.', facecolor=palette['quad']
            )
        )


    def draw_bend(ax, e):
        ax.add_patch(
            matplotlib.patches.Rectangle(
                (e['S'] - e['L'] / 2.0, 1000*e['POLES'][up]),  # (x,y)
                e['L'],  # width
                100,  # height
                hatch='/', facecolor=palette['bend']
            )
        )
        ax.add_patch(
            matplotlib.patches.Rectangle(
                (e['S'] - e['L'] / 2.0, -1000*e['POLES'][down]),  # (x,y)
                e['L'],  # width
                -100,  # height
                hatch='/', facecolor=palette['bend']
            )
        )


    beamline.sequence[beamline.sequence['CLASS'] == 'QUADRUPOLE'].apply(lambda e: draw_quad(ax,e), axis=1);
    beamline.sequence[beamline.sequence['CLASS'] == 'RBEND'].apply(lambda e: draw_bend(ax,e), axis=1);
    beamline.sequence[beamline.sequence['CLASS'] == 'SBEND'].apply(lambda e: draw_bend(ax,e), axis=1);
    # Beam pipe
    if planes == 'X':
        ax.plot(beamline.sequence['S'],
                 1000*beamline.sequence['APERTURE_X'].fillna(method='ffill'), 'k')
        ax.plot(beamline.sequence['S'],
                 -1000*beamline.sequence['APERTURE_X'].fillna(method='ffill'), 'k')
    elif planes == 'Y':
        ax.plot(beamline.sequence['S'],
                 1000*beamline.sequence['APERTURE_Y'].fillna(method='ffill'), 'k')
        ax.plot(beamline.sequence['S'],
                 -1000*beamline.sequence['APERTURE_Y'].fillna(method='ffill'), 'k')
    elif planes == 'both':
        ax.plot(beamline.sequence['S'],
                 1000*beamline.sequence['APERTURE_Y'].fillna(method='ffill'), 'k')
        ax.plot(beamline.sequence['S'],
                 -1000*beamline.sequence['APERTURE_X'].fillna(method='ffill'), 'k')