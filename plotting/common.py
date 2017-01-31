PALETTE = {}
PALETTE['solarized'] = {'base03': '#002b36',
                        'base02': '#073642',
                        'base01': '#586e75',
                        'base00': '#657b83',
                        'base0': '#839496',
                        'base1': '#93a1a1',
                        'base2': '#eee8d5',
                        'base3': '#fdf6e3',
                        'yellow': '#b58900',
                        'orange': '#cb4b16',
                        'red': '#dc322f',
                        'magenta': '#d33682',
                        'violet': '#6c71c4',
                        'blue': '#268bd2',
                        'cyan': '#2aa198',
                        'green': '#859900'
                        }

# Define default color palette
palette = PALETTE['solarized']
# Define "logical" colors
palette['quad'] = palette['blue']
palette['bend'] = palette['red']
palette['X'] = palette['cyan']
palette['Y'] = palette['orange']

def style_boxplot(bp, color):
    """Apply fancy styles to a matplotlib boxplot."""
    for box in bp['boxes']:
        box.set(color=color, linewidth=1)
        box.set(facecolor=color, alpha=0.4)
    for whisker in bp['whiskers']:
        whisker.set(color=color, linewidth=1, alpha=0.5)
    for cap in bp['caps']:
        cap.set(color=color, linewidth=1)
    for median in bp['medians']:
        median.set(color=color, linewidth=1)
    for flier in bp['fliers']:
        flier.set(marker='o', markersize=6, color=color, markeredgecolor='none', alpha=0.5)