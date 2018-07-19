"""
    manhattan_generator
    ~~~~~~~~~~~~~~~~~~~

    Help creating beautiful graphs of linkage results

    Version: 1.7.4

    Author: Louis-Philippe Lemieux Perreault

    Email:  louis-philippe.lemieux.perreault@statgen.org

"""


from __future__ import print_function
from __future__ import division

import os
import sys
import logging
import argparse

import numpy as np
import pandas as pd


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault"]
__license__ = "CC BY-NC 4.0"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"
__version__ = "1.7.4"


# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("manhattan-generator")


class DraggableAnnotation:
    """Creates draggable annotations for markers."""
    lock = None  # only one can be animated at a time

    def __init__(self, annot):
        """Creates an annotation which is draggable."""
        self.annot = annot
        self.press = None
        self.background = None

    def connect(self):
        """Connect to all the events we need."""
        self.cidpress = self.annot.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.annot.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.annot.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        "On button press, we will see if the mouse is over and store data."""
        if event.inaxes != self.annot.axes:
            return
        if DraggableAnnotation.lock is not None:
            return
        contains, attrd = self.annot.contains(event)
        if not contains:
            return
        x0, y0 = None, None
        if not hasattr(self.annot, "xyann"):
            # Quick fix for a deprecation in annotation...
            x0, y0 = self.annot.xytext
        else:
            x0, y0 = self.annot.xyann
        self.press = x0, y0, event.xdata, event.ydata
        DraggableAnnotation.lock = self

        # draw everything but the selected annotation and store the pixel
        # buffer
        canvas = self.annot.figure.canvas
        axes = self.annot.axes
        self.annot.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.annot.axes.bbox)

        # now redraw just the annotation
        axes.draw_artist(self.annot)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        """On motion we will move the annot if the mouse is over us."""
        if DraggableAnnotation.lock is not self:
            return
        if event.inaxes != self.annot.axes:
            return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        if not hasattr(self.annot, "xyann"):
            # Quick fix for a deprecation in annotation...
            self.annot.xytext = (x0+dx, y0+dy)
        else:
            self.annot.xyann = (x0+dx, y0+dy)

        canvas = self.annot.figure.canvas
        axes = self.annot.axes
        # restore the background region
        canvas.restore_region(self.background)

        # redraw just the current annotation
        axes.draw_artist(self.annot)

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """On release we reset the press data."""
        if DraggableAnnotation.lock is not self:
            return

        self.press = None
        DraggableAnnotation.lock = None

        # turn off the annot animation property and reset the background
        self.annot.set_animated(False)
        self.background = None

        # redraw the full figure
        self.annot.figure.canvas.draw()

    def disconnect(self):
        "Disconnect all the stored connection ids."""
        self.annot.figure.canvas.mpl_disconnect(self.cidpress)
        self.annot.figure.canvas.mpl_disconnect(self.cidrelease)
        self.annot.figure.canvas.mpl_disconnect(self.cidmotion)


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        """Creates a string representation of the message."""
        return self.message


def main():
    """The main method of the program."""
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Reading the input file for two point linkage
    two_point = None
    if args.twopoint is not None:
        two_point = read_input_file(args.twopoint, args.phys_pos_flag,
                                    args.use_pvalues_flag, args)

    # Reading the input file for multipoint linkage
    multi_point = None
    if args.multipoint is not None:
        multi_point = read_input_file(args.multipoint, args.phys_pos_flag,
                                      args.use_pvalues_flag, args)

    # Creating the plots
    create_manhattan_plot(two_point, multi_point, args)


def read_input_file(i_fn, use_bp, use_p, options):
    """Reads input file.

    Args:
        i_fn (str): the name of the input file.
        use_bp (bool): use physical position (bp) rather than genetic position?
        use_p (bool): use *p values* instead of *lod score*?
        options (argparse.Namespace): the options.


    Returns:
        pandas.DataFrame:  The array will contain the following names: ``chr``,
                           ``pos``, ``snp`` and ``conf``.

    This function reads any kind of input file, as long as the file is
    tab-separated and that it contains columns with the following headers:

    ======================  ===============================================
            Header                           Description
    ======================  ===============================================
    ``chr``                 The name of the chromosome
    ``snp``                 The name of the marker
    ``pos`` or ``cm``       The physical or genetic position (respectively)
    ``lod`` or ``p_value``  The confidence value (either *lod score* or *p
                            value*, respectively)
    ======================  ===============================================

    Note
    ----

        If there is a problem while reading the input file(s),
        a :py:class:`ProgramError` will be raised, and the program will be
        terminated.

    """
    # Reading the data
    csv_iterator = pd.read_csv(i_fn, sep="\t", chunksize=1e6, low_memory=False)
    data = next(csv_iterator).dropna()
    for chunk in csv_iterator:
        data = data.append(chunk.dropna(), ignore_index=True)

    # Checking we have the required column
    required_cols = {options.col_chr, options.col_name,
                     options.col_pos if use_bp else options.col_cm,
                     options.col_pvalue if use_p else options.col_lod}
    same_col = required_cols & set(data.columns)
    if same_col != required_cols:
        raise ProgramError("{}: missing columns {}".format(
            i_fn,
            ", ".join(required_cols - same_col),
        ))

    # Renaming the columns
    data = data.rename(columns={
        options.col_chr: "chrom",
        options.col_pos if use_bp else options.col_cm: "pos",
        options.col_name: "snp",
        options.col_pvalue if use_p else options.col_lod: "conf",
    })
    data = data[["chrom", "pos", "snp", "conf"]]

    # Encoding the chromosomes
    data["chrom"] = [encode_chr(chrom) for chrom in data.chrom]

    # Keeping only the required chromosome
    if options.only_chr is not None:
        data = data.loc[data.chrom == options.only_chr, :]
    else:
        data = data[~data.chrom.isin(options.exclude_chr)]

    # If p values, we modify
    if use_p:
        data["conf"] = -1 * np.log10(data.conf)

    # Ordering and returning
    return data.sort_values(by=["chrom", "pos"])


def create_manhattan_plot(twopoint, multipoint, args):
    """Creates the manhattan plot from marker data.

    Args:
        twopoint (pandas.DataFrame): the two point data
                                     (``None`` if not available).
        multipoint (pandas.DataFrame): the multipoint data
                                       (``None`` if not available).
        args (argparse.Namespace): the options and arguments of the program.

    Creates manhattan plots from two point or multipoint data. Two point
    results are shown in a manhattan plot using points (different color for
    each of the chromosomes). Multi point results are shown using lines.

    If both two and mutli point data are available, multi point results are
    shown above two point data.

    """
    import matplotlib as mpl
    if args.no_annotation:
        mpl.use("Agg")
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ProgramError("Could not import matplotlib. The most common "
                           "cause is that there is no available display, but "
                           "annotation has been asked for... Try using the "
                           "--no_annotation option.")

    if args.no_annotation:
        plt.ioff()

    # The available chromosomes
    available_chrom = []
    if args.twopoint is not None:
        available_chrom.append(sorted(twopoint.chrom.unique()))
    if args.multipoint is not None:
        available_chrom.append(sorted(multipoint.chrom.unique()))
    if len(available_chrom) == 1:
        available_chrom = available_chrom[0]
    else:
        if available_chrom[0] != available_chrom[1]:
            raise ProgramError("chromosomes are not the same for twopoint and "
                               "multipoint data")
        available_chrom = available_chrom[0]

    # Creating the figure
    figure = None
    figure = plt.figure(figsize=(args.graph_width, args.graph_height),
                        frameon=True)

    # Getting the maximum and minimum of the confidence value
    conf_min = [0.0]
    conf_max = []
    if args.twopoint is not None:
        conf_min.append(twopoint.conf.min())
        conf_max.append(twopoint.conf.max())
    if args.multipoint is not None:
        conf_min.append(multipoint.conf.min())
        conf_max.append(multipoint.conf.max())
    conf_min = min(conf_min)
    conf_max = max(conf_max)
    if args.max_ylim is not None:
        conf_max = args.max_ylim
    if args.min_ylim is not None:
        conf_min = args.min_ylim
    if args.no_negative_values or args.use_pvalues_flag:
        conf_min = 0.0

    # The chromosome spacing
    chrom_spacing = 25.0
    if args.phys_pos_flag:
        chrom_spacing = 25000000

    # Creating the ax and modify it
    ax = figure.add_subplot(111)
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("left")
    ax.set_ylabel("LOD")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    if args.only_chr is None:
        ax.set_xticks([])
        ax.set_xticklabels([])
    if args.use_pvalues_flag:
        ax.set_ylabel(r'$-\log_{10}$ (p value)', fontsize=args.label_text_size)
    else:
        ax.set_ylabel(args.graph_y_label, fontsize=args.label_text_size)
    ax.set_xlabel(args.graph_x_label, fontsize=args.label_text_size)
    ax.set_title(args.graph_title, fontsize=16, weight="bold")

    # Now plotting for each of the chromosome
    starting_pos = 0
    annots = []
    ticks = []
    for i, chrom in enumerate(available_chrom):
        chrom_twopoint = None
        chr_multipoint = None
        max_pos = []
        if args.twopoint is not None:
            chrom_twopoint = twopoint[twopoint.chrom == chrom]
            max_pos.append(chrom_twopoint.pos.max())
        if args.multipoint is not None:
            chr_multipoint = multipoint[multipoint.chrom == chrom]
            max_pos.append(chr_multipoint.pos.max())
        max_pos = max(max_pos)

        # The color of the points
        color = args.even_chromosome_color
        if i % 2 == 0:
            color = args.odd_chromosome_color
        multipoint_color = color

        # The box
        xmin = starting_pos - (chrom_spacing / 2)
        xmax = max_pos + starting_pos + (chrom_spacing / 2)
        if i % 2 == 1:
            ax.axvspan(xmin=xmin, xmax=xmax, color=args.chromosome_box_color)

        # The chromosome label
        ticks.append((xmin + xmax) / 2)

        # Plotting the twopoint
        if args.twopoint is not None:
            ax.plot(chrom_twopoint.pos + starting_pos, chrom_twopoint.conf,
                    marker="o", ms=args.point_size, mfc=color, mec=color,
                    ls="None")
            multipoint_color = args.multipoint_color

        # Plotting the multipoint
        if args.multipoint is not None:
            ax.plot(chr_multipoint.pos + starting_pos, chr_multipoint.conf,
                    ls="-", color=multipoint_color, lw=1.2)

        # Plotting the abline
        for abline_position in args.abline:
            ax.axhline(y=abline_position, color="black", ls="--", lw=1.2)
        if conf_min < 0:
            ax.axhline(y=0, color="black", ls="-", lw=1.2)

        # Plotting the significant markers
        if args.twopoint is not None:
            sig_mask = chrom_twopoint.conf >= args.significant_threshold
            ax.plot(chrom_twopoint.pos[sig_mask] + starting_pos,
                    chrom_twopoint.conf[sig_mask], marker="o", ls="None",
                    ms=args.significant_point_size, mfc=args.significant_color,
                    mec=args.significant_color)

            # If we want annotation
            if not args.no_annotation:
                for m_index, m in chrom_twopoint[sig_mask].iterrows():
                    # The confidence to write
                    the_conf = "{:.3f}".format(m.conf)
                    if args.use_pvalues_flag:
                        the_conf = str(10 ** (-1 * m.conf))

                    # The label of the annotation
                    label = "\n".join([m.snp, the_conf])

                    annot = ax.annotate(
                        label,
                        xy=(m.pos + starting_pos, m.conf),
                        xycoords="data",
                        size=10,
                        xytext=(m.pos + starting_pos, conf_max),
                        va="center",
                        bbox=dict(boxstyle="round", fc="white", ec="black"),
                        textcoords="data",
                        arrowprops=dict(arrowstyle="->", shrinkA=6, shrinkB=5),
                    )
                    annots.append(annot)

        # Changing the starting point for the next chromosome
        starting_pos = max_pos + starting_pos + chrom_spacing

    # Make the annotation draggable
    drs = []
    for annot in annots:
        dr = DraggableAnnotation(annot)
        dr.connect()
        drs.append(dr)

    # Setting the limits
    padding = 0.39
    if args.no_y_padding:
        padding = 0
    ax.set_ylim(conf_min - padding, conf_max + padding)
    ax.set_xlim(0 - chrom_spacing, starting_pos + chrom_spacing)

    # Putting the xticklabels
    if args.only_chr is None:
        ax.set_xticks(ticks)
        ax.set_xticklabels(available_chrom)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(args.axis_text_size)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(args.chr_text_size)

    # Saving or plotting the figure
    mpl.rcParams['savefig.dpi'] = args.dpi
    mpl.rcParams['ps.papersize'] = "auto"
    mpl.rcParams['savefig.orientation'] = "landscape"

    if args.no_annotation or (args.twopoint is None):
        # Annotation is for two-point only, se we save the figure
        plt.savefig(args.outFile_name + "." + args.graph_format,
                    bbox_inches="tight")
        if args.graph_format != "png":
            plt.savefig(args.outFile_name + ".png", bbox_inches="tight")
        if args.web:
            print(args.outFile_name + ".png")

    else:
        # There is some two-point data and annotation is asked, se we show
        # the figure
        plt.show()


def encode_chr(chromosome):
    """Encode a chromosome in integer format.

    Args:
        chromosome (str): the chromosome to encode in integer.

    Returns:
        int: the chromosome encoded in integer instead of string.

    This function encodes sex chromosomes, pseudo-autosomal regions and
    mitochondrial chromosomes in 23, 24, 25 and 26, respectively. If the
    chromosome is none of the above, the function returns the integer
    representation of the chromosome, if possible.

    Note
    ----

        If the chromosome is invalid, a :py:class:`ProgramError` will be
        raised, and the program terminated.

    Warning
    -------

        No check is done whether the chromosome is higher than 26 and below 1.
        As long as the chromosome is an integer or equal to ``X``, ``Y``,
        ``XY`` or ``MT``, no :py:class:`ProgramError` is raised.

    """
    try:
        return int(chromosome)

    except ValueError:
        chromosome = chromosome.upper()
        if chromosome == 'X':
            return 23
        elif chromosome == 'Y':
            return 24
        elif chromosome == 'XY':
            return 25
        elif chromosome == 'MT':
            return 26

        msg = "%(chromosome)s: not a valid chromosome" % locals()
        raise ProgramError(msg)


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): a :py:class:`Namespace` object containing
                                   the options of the program.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    if (args.max_ylim is not None) and (args.min_ylim is not None):
        if args.max_ylim <= args.min_ylim:
            msg = "Y max limit (%f) is <= Y min limit " \
                  "(%f)" % (args.max_ylim, args.min_ylim)
            raise ProgramError(msg)

    # The type of graph
    if (args.twopoint is None) and (args.multipoint is None):
        msg = "Meed to specify at least one graph type (option -t or -m)"
        raise ProgramError(msg)

    # Check for input file (two-point)
    if args.twopoint is not None:
        if not os.path.isfile(args.twopoint):
            msg = "%s: no such file or directory" % args.twopoint
            raise ProgramError(msg)

    # Check for input file (multipoint)
    if args.multipoint is not None:
        if not os.path.isfile(args.multipoint):
            msg = "%s: no such file or directory" % args.multipoint
            raise ProgramError(msg)

    try:
        args.abline = [float(i) for i in args.abline.split(',')]
    except ValueError:
        msg = "%s: not a valid LOD score (must be float)" % args.abline
        raise ProgramError(msg)

    # Checking that only one of --exclude-chr or --only-chr is used
    if args.exclude_chr is not None and args.only_chr is not None:
        msg = "Use only one of '--exclude-chr' or '--only-chr'"
        raise ProgramError(msg)

    # Checking if there are some chromosome to exclude
    if args.exclude_chr is None:
        args.exclude_chr = set()
    else:
        args.exclude_chr = {encode_chr(i) for i in args.exclude_chr.split(",")}

    if args.only_chr is not None:
        args.only_chr = encode_chr(args.only_chr)

    # Checking the graph title for unicode (python2)
    try:
        args.graph_title = unicode(args.graph_title, "utf-8")
        args.graph_x_label = unicode(args.graph_x_label, "utf-8")
        args.graph_y_label = unicode(args.graph_y_label, "utf-8")
    except NameError:
        pass


def parse_args():
    """Parses the command line options and arguments.

    Returns:
    argparse.Namespace: An object created by the :py:mod:`argparse` module. It
                        contains the values of the different options.

    ============================  =======  ====================================
         Options                   Type                   Description
    ============================  =======  ====================================
    ``--twopoint``                File     The input *file* for two-point
                                           linkage
    ``--multipoint``              File     The input *file* for multipoint
                                           linkage
    ``--output``                  String   The name of the ouput *file*
    ``--format``                  String   The format of the plot (ps, pdf
                                           png)
    ``--dpi``                     Int      The quality of the output (in dpi)
    ``--bp``                      Boolean  Use physical positions (bp) instead
                                           of genetic positions (cM).
    ``--use-pvalues``             Boolean  Use pvalues instead of LOD score
                                           requires to compute
                                           :math:`-log_{10}(pvalue)`
    ``--no-negative-values``      Boolean  Do not plot negative values
    ``--max-ylim``                Float    The maximal Y *value* to plot
    ``--min-ylim``                Float    The minimal Y *value* to plot
    ``--graph-title``             String   The *title* of the graph
    ``--graph-xlabel``            String   The *text* for the x label
    ``--graph-ylabel``            String   The *text* for the y label
    ``--graph-width``             Int      The *width* of the graph, in
                                           inches
    ``--graph-height``            Int      The *height* of the graph, in
                                           inches
    ``--point-size``              Float    The *size* of each points.
    ``--significant-point-size``  Float    The *size* of each significant
                                           points
    ``--abline``                  String   The y *value* where to create a
                                           horizontal line, separated by a
                                           comma
    ``--significant-threshold``   Float    The significant threshold for
                                           linkage
    ``--no-annotation``           Boolean  Do not draw annotation (SNP names)
                                           for the significant results
    ``--chromosome-box-color``    String   The *color* for the box surrounding
                                           even chromosome numbers
    ``--even-chromosome-color``   String   The *color* for the box surrounding
                                           even chromosome numbers
    ``--odd-chromosome-color``    String   The *color* for the box surrounding
                                           odd chromosome numbers
    ``--multipoint-color``        String   The *color* for the multipoint plot
    ``--significant-color``       String   The *color* for points representing
                                           significant linkage
    ============================  =======  ====================================

    Note
    ----

        No option check is done here (except for the one automatically done
        by :py:mod:`argparse`. Those need to be done elsewhere
        (see :py:func:`check_args`).

    """
    # Creating the parser object
    parser = argparse.ArgumentParser(
        description="This script produces nice Manhattan plots for either "
                    "linkage or GWAS results.",
    )

    # Adding the version option
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s {}".format(__version__),
    )

    # The input options
    group = parser.add_argument_group(
        "Input Options",
        "Options for the input file(s) (name of the file, type of graph, "
        "etc.). Note that for GWAS results, only the '--twopoint' option "
        "should be used.",
    )

    # The input file (for two point)
    group.add_argument(
        "--twopoint", type=str, metavar="FILE",
        help="The input FILE for two-point linkage.",
    )

    # The input file (for multipoint)
    group.add_argument(
        "--multipoint", type=str, metavar="FILE",
        help="The input FILE for multipoint linkage.",
    )

    # The column options
    group = parser.add_argument_group(
        "Column Options",
        "The name of the different columns in the input file(s).",
    )

    # The chromosome column
    group.add_argument(
        "--col-chr", type=str, metavar="COL", default="chr",
        help="The name of the column containing the chromosomes "
             "[Default: %(default)s].",
    )

    # The marker name column
    group.add_argument(
        "--col-name", type=str, metavar="COL", default="name",
        help="The name of the column containing the marker names "
             "[Default: %(default)s].",
    )

    # The marker position column
    group.add_argument(
        "--col-pos", type=str, metavar="COL", default="pos",
        help="The name of the column containing the marker positions "
             "[Default: %(default)s].",
    )

    # The marker cM column
    group.add_argument(
        "--col-cm", type=str, metavar="COL", default="cm",
        help="The name of the column containing the marker cM "
             "[Default: %(default)s].",
    )

    # The marker p value
    group.add_argument(
        "--col-pvalue", type=str, metavar="COL", default="p_value",
        help="The name of the column containing the marker p values "
             "[Default: %(default)s]",
    )

    # The marker lod score
    group.add_argument(
        "--col-lod", type=str, metavar="COL", default="lod",
        help="The name of the column containing the marker LOD score "
             "[Default: %(default)s]",
    )

    # The output options
    group = parser.add_argument_group(
        "Graph Output Options",
        "Options for the ouput file (name of the file, type of graph, etc.).",
    )

    # The output file name
    group.add_argument(
        "-o", "--output", dest="outFile_name", type=str, default="manhattan",
        metavar="NAME",
        help="The NAME of the ouput file [Default: %(default)s].",
    )

    # The type of the graph (png, ps or pdf)
    format_choices = ["ps", "pdf", "png", "eps"]
    group.add_argument(
        "-f", "--format", dest="graph_format", type=str, default="png",
        metavar="FORMAT", choices=format_choices,
        help="The FORMAT of the plot ({}) "
             "[Default: %(default)s].".format(", ".join(format_choices)),
    )

    group.add_argument(
        "--web", action="store_true",
        help="Always write a PNG file for web display, and return the path of "
             "the PNG file.",
    )

    group.add_argument(
        "--dpi", type=int, default=600, metavar="INT",
        help="The quality of the output (in dpi) [Default: %(default)d].",
    )

    # The graph type options
    group = parser.add_argument_group(
        "Graph Options",
        "Options for the graph type (two-point, multipoint, etc.).",
    )

    # Use physical position instead of genetic posiition
    group.add_argument(
        "--bp", dest="phys_pos_flag", action="store_true",
        help="Use physical positions (bp) instead of genetic positions (cM).",
    )

    # Using p values instead of LOD score
    group.add_argument(
        "--use-pvalues", dest='use_pvalues_flag', action="store_true",
        help="Use pvalues instead of LOD score. Requires to compute "
             "-log10(pvalue).",
    )

    # Exclude some chromosomes
    group.add_argument(
        "--exclude-chr", metavar="STRING",
        help="Exclude those chromosomes (list of chromosomes, separated by a "
             "coma) [Default: None].",
    )
    group.add_argument(
        "--only-chr", metavar="CHR",
        help="Print only the results for a single chromosome. The xaxis will "
             "hence show the positions/cm instead of the chromosome number.",
    )

    # The graph presentation options
    group = parser.add_argument_group(
        "Graph Presentation Options",
        "Options for the graph presentation (title, axis label, etc.).",
    )

    # print negative values
    group.add_argument(
        "--no-negative-values", action="store_true",
        help="Do not plot negative values.",
    )

    # The maximal y limit of the graph
    group.add_argument(
        "--max-ylim", type=float, metavar="FLOAT",
        help="The maximal Y value to plot [Default: maximum of max(LOD) "
             "and 1+significant-threshold].",
    )

    # The minimal y limit of the graph
    group.add_argument(
        "--min-ylim", type=float, default=-2.0, metavar="FLOAT",
        help="The minimal Y value to plot [Default: %(default).1f].",
    )

    # Do we want padding?
    group.add_argument(
        "--no-y-padding", action="store_true",
        help="Do not add Y padding to the Y limit",
    )

    # The graph's title
    group.add_argument(
        "--graph-title", type=str, dest='graph_title', default="",
        metavar="TITLE",
        help="The TITLE of the graph [Default: empty].",
    )

    # The graph's x label
    group.add_argument(
        "--graph-xlabel", dest='graph_x_label', type=str, default="Chromosome",
        metavar="TEXT",
        help="The TEXT for the x label. [Default: %(default)s].",
    )

    # The graph's y label
    group.add_argument(
        "--graph-ylabel", dest='graph_y_label', type=str, default="LOD",
        metavar="TEXT",
        help="The TEXT for the y label. [Default: %(default)s].",
    )

    # The graph width
    group.add_argument(
        "--graph-width", type=int, default=14, metavar="WIDTH",
        help="The WIDTH of the graph, in inches [Default: %(default)d].",
    )

    # The graph height
    group.add_argument(
        "--graph-height", type=int, default=7, metavar="HEIGHT",
        help="The HEIGHT of the graph, in inches [Default: %(default)d].",
    )

    # The size of each point
    group.add_argument(
        "--point-size", type=float, default=2.1, metavar="SIZE",
        help="The SIZE of each points [Default: %(default).1f].",
    )

    # The size of each significant point
    group.add_argument(
        "--significant-point-size", type=float, default=4.5, metavar="SIZE",
        help="The SIZE of each significant points [Default: %(default).1f].",
    )

    # The ablines positions
    group.add_argument(
        "--abline", type=str, default="3,-2", metavar="POS1,POS2,...",
        help="The y value where to create a horizontal line, separated by a "
             "comma [Default: %(default)s].",
    )

    # The significant threshold
    group.add_argument(
        "--significant-threshold", type=float, default=3.0, metavar="FLOAT",
        help="The significant threshold for linkage or association "
             "[Default: %(default).1f]",
    )

    # The annotation flag
    group.add_argument(
        "--no-annotation", action="store_true",
        help="Do not draw annotation (SNP names) for the significant results.",
    )

    # The size of the text
    group.add_argument(
        "--axis-text-size", type=int, default=12, metavar="INT",
        help="The axis font size [Default: %(default)d]",
    )

    group.add_argument(
        "--chr-text-size", type=int, default=12, metavar="INT",
        help="The axis font size [Default: %(default)d]",
    )

    group.add_argument(
        "--label-text-size", type=int, default=12, metavar="INT",
        help="The axis font size [Default: %(default)d]",
    )

    # The graph color options
    group = parser.add_argument_group(
        "Graph Colors Options",
        "Options for the graph colors.",
    )

    group.add_argument(
        "--chromosome-box-color", type=str, default="#E5E5E5", metavar="COLOR",
        help="The COLOR for the box surrounding even chromosome numbers "
             "[Default: %(default)s].",
    )

    group.add_argument(
        "--even-chromosome-color", type=str, default="#1874CD",
        metavar="COLOR",
        help="The COLOR for the box surrounding even chromosome numbers "
             "[Default: %(default)s].",
    )

    group.add_argument(
        "--odd-chromosome-color", type=str, default="#4D4D4D", metavar="COLOR",
        help="The COLOR for the box surrounding odd chromosome numbers "
             "[Default: %(default)s].",
    )

    group.add_argument(
        "--multipoint-color", type=str, default="#FF8C00", metavar="COLOR",
        help="The COLOR for the multipoint plot [Default: %(default)s].",
    )

    group.add_argument(
        "--significant-color", type=str, default="#FF0000", metavar="COLOR",
        help="The COLOR for points representing significant linkage "
             "[Default: %(default)s].",
    )

    return parser.parse_args()


def safe_main():
    """A main function that catches errors."""
    try:
        main()

    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    except ProgramError as e:
        logger.critical(e.message)


if __name__ == "__main__":
    safe_main()
