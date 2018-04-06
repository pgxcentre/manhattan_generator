# Manhattan Plot Generator

`manhattan_generator` is a python tool to create beautiful Manhattan plots.


## Dependencies

Here are the dependencies of the tool:

- [Python](http://python.org/) version 2.7 or 3.4 or latest
- [numpy](http://www.numpy.org/) version 1.8.0 or latest
- [matplotlib](http://matplotlib.org/) version 1.3.1 or latest
- [pandas](http://pandas.pydata.org/) version 0.17 or latest


## Installation

We recommend installing the tool in a Python virtual environment.

`manhattan_generator` should work on Windows and MacOS, even though it hasn't
been fully tested for full compatibility.


## Basic usage

```console
$ manhattan_generator --help
usage: manhattan_generator [-h] [-v] [--twopoint FILE] [--multipoint FILE]
                           [--col-chr COL] [--col-name COL] [--col-pos COL]
                           [--col-cm COL] [--col-pvalue COL] [--col-lod COL]
                           [-o NAME] [-f FORMAT] [--web] [--dpi INT] [--bp]
                           [--use-pvalues] [--exclude-chr STRING]
                           [--no-negative-values] [--max-ylim FLOAT]
                           [--min-ylim FLOAT] [--no-y-padding]
                           [--graph-title TITLE] [--graph-xlabel TEXT]
                           [--graph-ylabel TEXT] [--graph-width WIDTH]
                           [--graph-height HEIGHT] [--point-size SIZE]
                           [--significant-point-size SIZE]
                           [--abline POS1,POS2,...]
                           [--significant-threshold FLOAT] [--no-annotation]
                           [--axis-text-size INT] [--chr-text-size INT]
                           [--label-text-size INT]
                           [--chromosome-box-color COLOR]
                           [--even-chromosome-color COLOR]
                           [--odd-chromosome-color COLOR]
                           [--multipoint-color COLOR]
                           [--significant-color COLOR]

This script produces nice Manhattan plots for either linkage or GWAS results.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Input Options:
  Options for the input file(s) (name of the file, type of graph, etc.).
  Note that for GWAS results, only the '--twopoint' option should be used.

  --twopoint FILE       The input FILE for two-point linkage.
  --multipoint FILE     The input FILE for multipoint linkage.

Column Options:
  The name of the different columns in the input file(s).

  --col-chr COL         The name of the column containing the chromosomes
                        [Default: chr].
  --col-name COL        The name of the column containing the marker names
                        [Default: name].
  --col-pos COL         The name of the column containing the marker positions
                        [Default: pos].
  --col-cm COL          The name of the column containing the marker cM
                        [Default: cm].
  --col-pvalue COL      The name of the column containing the marker p values
                        [Default: p_value]
  --col-lod COL         The name of the column containing the marker LOD score
                        [Default: lod]

Graph Output Options:
  Options for the ouput file (name of the file, type of graph, etc.).

  -o NAME, --output NAME
                        The NAME of the ouput file [Default: manhattan].
  -f FORMAT, --format FORMAT
                        The FORMAT of the plot (ps, pdf, png, eps) [Default:
                        png].
  --web                 Always write a PNG file for web display, and return
                        the path of the PNG file.
  --dpi INT             The quality of the output (in dpi) [Default: 600].

Graph Options:
  Options for the graph type (two-point, multipoint, etc.).

  --bp                  Use physical positions (bp) instead of genetic
                        positions (cM).
  --use-pvalues         Use pvalues instead of LOD score. Requires to compute
                        -log10(pvalue).
  --exclude-chr STRING  Exclude those chromosomes (list of chromosomes,
                        separated by a coma) [Default: None].

Graph Presentation Options:
  Options for the graph presentation (title, axis label, etc.).

  --no-negative-values  Do not plot negative values.
  --max-ylim FLOAT      The maximal Y value to plot [Default: maximum of
                        max(LOD) and 1+significant-threshold].
  --min-ylim FLOAT      The minimal Y value to plot [Default: -2.0].
  --no-y-padding        Do not add Y padding to the Y limit
  --graph-title TITLE   The TITLE of the graph [Default: empty].
  --graph-xlabel TEXT   The TEXT for the x label. [Default: Chromosome].
  --graph-ylabel TEXT   The TEXT for the y label. [Default: LOD].
  --graph-width WIDTH   The WIDTH of the graph, in inches [Default: 14].
  --graph-height HEIGHT
                        The HEIGHT of the graph, in inches [Default: 7].
  --point-size SIZE     The SIZE of each points [Default: 2.1].
  --significant-point-size SIZE
                        The SIZE of each significant points [Default: 4.5].
  --abline POS1,POS2,...
                        The y value where to create a horizontal line,
                        separated by a comma [Default: 3,-2].
  --significant-threshold FLOAT
                        The significant threshold for linkage or association
                        [Default: 3.0]
  --no-annotation       Do not draw annotation (SNP names) for the significant
                        results.
  --axis-text-size INT  The axis font size [Default: 12]
  --chr-text-size INT   The axis font size [Default: 12]
  --label-text-size INT
                        The axis font size [Default: 12]

Graph Colors Options:
  Options for the graph colors.

  --chromosome-box-color COLOR
                        The COLOR for the box surrounding even chromosome
                        numbers [Default: #E5E5E5].
  --even-chromosome-color COLOR
                        The COLOR for the box surrounding even chromosome
                        numbers [Default: #1874CD].
  --odd-chromosome-color COLOR
                        The COLOR for the box surrounding odd chromosome
                        numbers [Default: #4D4D4D].
  --multipoint-color COLOR
                        The COLOR for the multipoint plot [Default: #FF8C00].
  --significant-color COLOR
                        The COLOR for points representing significant linkage
                        [Default: #FF0000].
```


## Example

As an example, we used the dataset publicly provided by Wood *et al.* 2014 (doi: [10.1038/ng.3097](http://dx.doi.org/10.1038/ng.3097)], as
part of the GIANT consortium for height (available
[here](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2014_Height)).

```bash
manhattan_generator \
    --twopoint ../data/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.grch37.txt \
    --col-chr Chromosome \
    --col-name MarkerName \
    --col-pos Position \
    --col-pvalue p \
    --bp \
    --use-pvalues \
    --abline 95 \
    --significant-threshold 95 \
    --no-annotation \
    --significant-point-size 2 \
    --point-size 1 \
    --graph-title "GIANT height (Wood et al. 2014, public release)" \
    --chr-text-size 10 \
    --exclude-chr 23,24
```

<img src=https://raw.github.com/pgxcentre/manhattan_generator/master/example_giant.png width=728 />
