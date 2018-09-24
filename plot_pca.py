import matplotlib.pyplot as plt
import seaborn as sns

def plot_pca(outfile, proj, pcs=(0, 1), labels=None, label_points=True,
             levels=None, colors=None, legend=True, s=100):
    """
    Plots a PCA projection along two selected principal components.

    Parameters
    ----------
    outfile : str
        The file to write the plot to.
    proj : np.ndarray
        The matrix of PCA-projected replicates.
    pcs : Tuple[int]
        Which two (zero-indexed) principle components should be plotted.
    labels : Optional[List[str]]
        String names identifying the replicates (the rows of ``proj``). Pass
        None to simply label them with their row index within ``proj``.
    label_points : bool
        Pass True to annotate each point with its label.
    levels : Optional[List[str]]
        The "level" for each replicate (each row of ``proj``). Each "level" gets
        one color and one entry in the legend. If None is passed each replicate
        gets its own level (``levels = labels``).
    colors : Optional[Dict[str, str]]
        Mapping from levels as strings to the color to use for that level. Pass
        None to use randomly assigned colors.
    legend : bool
        Pass True to include a legend.
    s : float
        The area of the points to plot on the scatterplot.
    """
    # clear plot area
    plt.clf()

    # resolve labels
    if labels is None:
        labels = [str(i) for i in range(len(proj))]

    # resolve levels
    if levels is None:
        levels = labels

    # resolve colors
    if colors is None:
        palette = sns.color_palette('husl', len(levels))
        colors = {levels[i]: palette[i] for i in range(len(levels))}

    # plot, keeping track of the handles we encouter
    handles = [plt.scatter(proj[i, pcs[0]], proj[i, pcs[1]],
                           c=colors[levels[i]], s=s)
               for i in range(len(proj))]

    # add legend
    if legend:
        legend_labels = list(set(levels))
        legend_handles = [handles[levels.index(l)] for l in legend_labels]
        lgd = plt.legend(legend_handles, legend_labels, scatterpoints=1,
                         loc='upper left', bbox_to_anchor=(1, 1.05))

    # add labels to points
    if label_points:
        for i in range(len(proj)):
            plt.annotate(labels[i], (proj[i, pcs[0]], proj[i, pcs[1]]))

    # add labels to axes
    plt.xlabel('PC %i' % (pcs[0] + 1))
    plt.ylabel('PC %i' % (pcs[1] + 1))

    # save figure
    if legend:
        plt.savefig(outfile, bbox_inches='tight')
    else:
        plt.savefig(outfile, bbox_inches='tight')

