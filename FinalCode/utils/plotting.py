import matplotlib.pyplot as plt


def pretty_ax(ax: plt.Axes, linew: int = 1):
    # Hide the right and top spines
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(linew)


def transform_plot_ax(
    ax: plt.Axes,
    legend_title: str,
    ftsize: int = 15,
    leg_ftsize: int = 12,
    linew: int = 4,
    remove_ticks: bool = False,
):
    pretty_ax(ax)
    ax.legend(
        frameon=False,
        title=legend_title,
        fontsize=leg_ftsize,
        title_fontsize=leg_ftsize,
    )
    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(linew)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=ftsize)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=ftsize)
    ax.set_xlabel(ax.get_xlabel(), fontsize=ftsize)
    ax.set_ylabel(ax.get_ylabel(), fontsize=ftsize)
