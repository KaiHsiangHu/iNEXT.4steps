<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.4steps (R package)

<h5 align="right">
Latest version: 2024-02-24
</h5>
<font color="394CAE">
<h3 color="394CAE" style="font-weight: bold">
Introduction to iNEXT.4steps (R package): Excerpt from iNEXT.4steps
User’s Guide
</h3>
</font> <br>
<h5>
<b>Anne Chao, Kai-Hsiang Hu</b> <br><br> <i>Institute of Statistics,
National Tsing Hua University, Hsin-Chu, Taiwan 30043</i>
</h5>

<br> `iNEXT.4steps` (iNterpolation and EXTrapolation for four steps of
biodiversity) is an original R package which provide an easy complete
biological analysis computation. In Chao et al. (2020) paper, they
proposed a complete biological analysis process:

-   `STEP1`. Sample completeness profiles.

-   `STEP2`. Size-based rarefaction and extrapolation analysis and the
    asymptotic diversity profile.

-   `STEP3`. Non-asymptotic coverage-based rarefaction and extrapolation
    analysis.

-   `STEP4`. Evenness among species abundances.

These are the compositions of `iNEXT.4steps`. Here we will also
introduce functions for `STEP1` and `STEP4`, particularly. If you want
to grasp the functions of `STEP2` or `STEP3` only, then please refer to
the package
[iNEXT.3D](https://cran.r-project.org/web/packages/iNEXT.3D/index.html)
from CRAN for more details. `iNEXT.3D` features two statistical analyses
(non-asymptotic and asymptotic):

1.  A non-asymptotic approach based on interpolation and extrapolation
    for 3D diversity (i.e., Hill-Chao numbers)

computes the estimated 3D diversity for standardized samples with a
common sample size or sample completeness. This approach aims to compare
diversity estimates for equally-large (with a common sample size) or
equally-complete (with a common sample coverage) samples; it is based on
the seamless rarefaction and extrapolation (R/E) sampling curves of
Hill-Chao numbers for q = 0, 1 and 2.

1.  An asymptotic approach to infer asymptotic 3D diversity (i.e.,
    Hill-Chao numbers)

`iNEXT.4steps` provides four functions for users to compute and
visualize sample completeness curves and evenness curves: `Completeness`
and `Evenness` for computation, and `ggCompleteness` and `ggEvenness`
for visualization. The most comprehensive function `iNEXT4steps` gathers
`iNEXT3D`, `AO3D`, `Completeness` and `Evenness` into a summary table
and a figure. In this document, we will give a quick introduction
demonstrating how to run these functions. Detailed information about
these function settings is provided in the `iNEXT.4steps` manual, which
will be submitted to \[CRAN\]. See Chao et al. (2020) to get more
details.

## How to cite

If you publish your work based on results from iNEXT.4steps package, you
should make references to the following methodology paper and the
package:

-   Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B.
    Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello,
    and R. K. Colwell (2020). Quantifying sample completeness and
    comparing diversities among assemblages. Ecological Research, 35,
    292-314.

## SOFTWARE NEEDED TO RUN INEXT IN R

-   Required: [R](https://cran.r-project.org/)
-   Suggested: [RStudio
    IDE](https://posit.co/products/open-source/rstudio/#Desktop)

## HOW TO RUN INEXT.4STEPS:

The `iNEXT.4steps` package will be submitted to \[CRAN\] or can be
downloaded from Anne Chao’s
[iNEXT.4steps_github](https://github.com/AnneChao/iNEXT.4steps). For a
first-time installation, some additional packages must be installed and
loaded; see package manual.

``` r
## install iNEXT.4steps package from CRAN
# install.packages("iNEXT.4steps")  # Coming soon

## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('AnneChao/iNEXT.4steps')

## import packages
library(iNEXT.4steps)
```

An online version of iNEXT.4steps
(<https://chao.shinyapps.io/iNEXT_4steps/>) is also available for users
without R background.

## SIMPLE STEPS BRIEF

-   **Step 1: Sample completeness profile**

`Sample Completeness` represent the proportion of observed species in
the population (Chao et al., 2020). Usually, the sampling data
represents the abundant species in the population so that we will ignore
the rare species. Here we will use Turing’s sample coverage theory to
reconstruct the population proportion. Besides, sample completeness can
correspond to order q, which is an weight index. When order q tends to
zero, then we will give more weight to rare species. If order q tends to
unity, then we will equally treat each species, which is also called
sample coverage at unity. In contrast, if order q tends to larger than
unity, we will give more weights to abundant species. By sample
completeness, we can easily plot the estimated curve with respect to
order q and associated 95% confidence interval.

-   **Step 2.1 and step 3: Size-based and coverage-based Interpolation
    and Extrapolation**

`Interpolation and Extrapolation (iNEXT)` focuses on three measures of
order q: species richness (q = 0), Shannon diversity (q = 1, the
exponential of Shannon entropy), and Simpson diversity (q = 2, the
inverse of Simpson concentration) (Chao and Jost, 2012; Chao et
al. 2014). For each diversity measures, `iNEXT` uses observed sample to
compute expected diversity estimates and associated 95% confidence
intervals according two different unit types of rarefaction and
extrapolation (R/E):

1.  Sample-size-based R/E sampling curves versus diversity in each order
    q.

2.  Coverage-based R/E sampling curves versus diversity in each order q.

A unified framework based on Hill numbers (for TD) and their
generalizations (Hill-Chao numbers, for PD and FD) is adopted to
quantify 3D. In this framework, TD quantifies the effective number of
species, PD quantifies the effective total branch length, mean-PD (PD
divided by tree depth) quantifies the effective number of lineages, and
FD quantifies the effective number of virtual functional groups (or
functional “species”). Thus, TD, mean-PD, and FD are all in the same
units of species/lineage equivalents and can be meaningfully compared.
For more particular usage about `iNEXT.3D`, please refer to Chao et
al. (2021).

-   **Step 2.2: Asymptotic diversity profile**

`Asymptotic Diversity (or called Hill numbers)` computes the estimated
asymptotic 3D diversity and also plots 3D diversity profiles
(q-profiles) for q between 0 and 2, in comparison with the observed
diversity. Typically, the asymptotic estimates for q ≥ 1 are reliable,
but for q \< 1 (especially for q = 0, species richness), the asymptotic
estimates represent only lower bounds. `Asymptotic Diversity` also
features a time-profile (which depicts the observed and asymptotic
estimate of PD or mean PD with respect to reference times), and a
tau-profile (which depicts the observed and asymptotic estimate of FD
with respect to threshold level tau).

-   **Step 4: Evenness profile**

`Evenness` is an function to compute whether a assemblage is uniform or
not. There are five main classes according to different transformation
by species and diversity (Chao and Ricotta, 2019). In these five
classes, they all have range from zero to one. When the value reaches to
minimum (zero), it means that the assemblage tends to uneven. On the
contrary, when the value reaches to maximum (one), it means that the
assemblage tends to even. Evenness considers different order q under
each classes. When order q tends to zero, we will focus on rare species.
In other sides, when order q tends to far from zero, then we will give
more weights on abundant species. For each sampling data, it usually has
biase because of unobserved rare species. Here, we propose a
“standardized sample coverage” (default is `Cmax`, which means the
minimum sample coverage among all samples extrapolated to double
reference sizes) as a judged criterion. According to this criterion, we
can also plot the Evenness curves with order q and associated 95%
confidence interval.

## <span style="color:red;">DATA INPUT FORMAT</span>

### Species abundance/incidence data format

For `iNEXT.3D` package, information on species identity (or any unique
identification code) and assemblage affiliation is required for PD and
FD. Two types of species abundance/incidence data are supported:

1.  Individual-based abundance data (`datatype = "abundance"`): When
    there are multiple assemblages, in addition to the assemblage/site
    names (as column names) and the species names (as row names),
    species abundance data (reference sample) can be input as a species
    (in rows) by assemblage (in columns) matrix/data.frame or a list of
    species abundance vectors. In the special case that there is only
    one assemblage, all data should be read in one column.

2.  Sampling-unit-based incidence data: There are two kinds of input
    data.  

<!-- -->

1.  Sampling-unit-based incidence data: Incidence-raw data
    (`datatype = "incidence_raw"`): for each assemblage, input data for
    a reference sample consist of a species-by-sampling-unit matrix, in
    addition to the sampling-unit names (as column names) and the
    species names (as row names). When there are N assemblages, input
    data consist of N lists of matrices, and each matrix is a
    species-by-sampling-unit matrix. Each element in the incidence raw
    matrix is 1 for a detection, and 0 for a non-detection. Input a
    matrix which combines data for all assemblages is allowed, but the
    argument `nT` in the function `iNEXT3D` must be specified so that
    the number of sampling units in each assemblage is specified.

2.  Sampling-unit-based incidence-frequency data
    (`datatype = "incidence_freq"`): input data for each assemblage
    consist of species sample incidence frequencies (i.e., row sums of
    the corresponding incidence raw matrix). When there are N
    assemblages, input data consist of an (S+1) by N matrix, or N lists
    of species incidence frequencies. The first entry of each
    column/list must be the total number of sampling units, followed by
    the species incidence frequencies.

For example, the dataset `Brazil_abun_data` included in the `iNEXT.3D`
package consists of species sample abundances of two
assemblages/habitats: “Edge” and “Interior”. Run the following code to
view the first 15 rows of the abundance data.

``` r
data("Brazil_abun_data")
Brazil_abun_data
```

                             Edge Interior
    Carpotroche_brasiliensis   11       21
    Astronium_concinnum       110       11
    Astronium_graveolens       36        7
    Spondias_macrocarpa        12        1
    Spondias_venulosa           2        0
    Tapirira_guianensis         7        1
    Thyrsodium_spruceanum      11       11
    Anaxagorea_silvatica        1       13
    Annona_acutiflora           1        1
    Annona_cacans               0        2
    Annona_dolabripetala        3        3
    Annona_sp                   0        1
    Duguetia_chrysocarpa        1        1
    Ephedranthus_sp1            1        0
    Ephedranthus_sp2            0        1

We use incidence frequency data (`Woody_plants`) collected from four
sites, namely `"Upper_cloud"`, `"Lower_cloud"`, `"Lowland"`, and
`"Monsoon"`, as an example. Note that the first row of incidence
frequency data should be the total sampling units for each assemblage.
Run the following code to view the first 16 rows and first 2 columns of
the incidence frequency data.

``` r
data("Woody_plants")
Woody_plants
```

                                                                          Upper_cloud Lower_cloud
    plots                                                                         153         203
    Abelia_chinensis_R._Br._var._ionandra_(Hayata)_Masam.                           1           0
    Abies_kawakamii_(Hayata)_Ito                                                    2           1
    Acacia_confusa_Merr.                                                            0           0
    Acer_albopurpurascens_Hayata                                                    1           2
    Acer_kawakamii_Koidzumi                                                         7          36
    Acer_morrisonense_Hayata                                                       25          25
    Acer_palmatum_Thunb._var._pubescens_Li                                         12           1
    Acer_serrulatum_Hayata                                                          6          26
    Actinidia_chinensis_Planch._var._setosa_Li                                      0           0
    Adinandra_formosana_Hayata                                                     12          16
    Adinandra_lasiostyla_Hayata                                                    34          15
    Aeschynanthus_acuminatus_Wall._ex_A._DC.                                        0           0
    Ailanthus_altissima_(Miller)_Swingle_var._tanakai_(Hayata)_Sasaki               0           1
    Akebia_trifoliata_(Thunb.)_Koidz._subsp._australis_(Diels)_T._Shimizu           0           1
    Alangium_chinense_(Lour.)_Rehder                                                0           1

### Phylogenetic tree format for PD

To perform PD analysis, the phylogenetic tree (in Newick format) spanned
by species observed in the pooled data is required. For the data
`Brazil_abun_data`, the phylogenetic tree for all observed species
(including species in both Edge and Interior habitats) is stored in the
file `Brazil_phylo_tree`. A partial list of the tip labels and node
labels are shown below.

``` r
data("Brazil_phylo_tree")
Brazil_phylo_tree

Phylogenetic tree with 425 tips and 205 internal nodes.

Tip labels:
  Carpotroche_brasiliensis, Casearia_ulmifolia, Casearia_sp4, Casearia_sylvestris, Casearia_sp2, Casearia_sp3, ...
Node labels:
  magnoliales_to_asterales, poales_to_asterales, , , , celastrales_to_malpighiales, ...

Rooted; includes branch lengths.
```

### Species pairwise distance matrix format for FD

To perform FD analysis, the species-pairwise distance matrix (Gower
distance computed from species traits) for species observed in the
pooled data is required in a matrix/data.frame format. For the data
`Brazil_abun_data`, the distance matrix for all species (including
species in both Edge and Interior habitats) is stored in the file
`Brazil_rainforest_dist_matrix`. The distance matrix for the first 3
Brazil rainforest tree species is shown below.

``` r
data("Brazil_distance_matrix")
Brazil_distance_matrix
```

                             Carpotroche_brasiliensis Astronium_concinnum Astronium_graveolens
    Carpotroche_brasiliensis                    0.000               0.522                0.522
    Astronium_concinnum                         0.522               0.000                0.000
    Astronium_graveolens                        0.522               0.000                0.000

## <span style="color:red;">MAIN FUNCTION iNEXT4steps()</span>

We first describe the main function `iNEXT4steps()` with default
arguments:

``` r
iNEXT4steps(data, diversity = "TD", q = seq(0, 2, 0.2), datatype = "abundance", 
            nboot = 50, conf = 0.95, nT = NULL, PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD', 
            FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50, details = FALSE)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style='text-align: left;'>

1.  For `datatype = 'abundance'`, data can be input as a vector of
    species abundances (for a single assemblage), matrix/data.frame
    (species by assemblages), or a list of species abundance
    vectors. (2) For `datatype = 'incidence_freq'`, data can be input as
    a vector of incidence frequencies (for a single assemblage),
    matrix/data.frame (species by assemblages), or a list of incidence
    frequencies; the first entry in all types of input must be the
    number of sampling units in each assemblage. (3) For
    `datatype = 'incidence_raw'`, data can be input as a list of
    matrix/data.frame (species by sampling units); data can also be
    input as a matrix/data.frame by merging all sampling units across
    assemblages based on species identity; in this case, the number of
    sampling units (nT, see below) must be input.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    diversity
    </td>
    <td style="text-align: left;">
    selection of diversity type: `TD` = Taxonomic diversity, `PD` =
    Phylogenetic diversity, and `FD` = Functional diversity.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    q
    </td>
    <td style="text-align: left;">
    a numerical vector specifying the diversity orders for q-profile
    output. Default is <code>seq(0, 2, by = 0.2)</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    datatype
    </td>
    <td style="text-align: left;">
    data type of input data: individual-based abundance data (datatype =
    ‘abundance’), sampling-unit-based incidence frequencies data
    (datatype = ‘incidence_freq’), or species by sampling-units
    incidence matrix (datatype = ‘incidence_raw’) with all entries being
    0 (non-detection) or 1 (detection).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nboot
    </td>
    <td style="text-align: left;">
    a positive integer specifying the number of bootstrap replications
    when assessing sampling uncertainty and constructing confidence
    intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    conf
    </td>
    <td style="text-align: left;">
    a positive number \< 1 specifying the level of confidence interval.
    Default is 0.95.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nT
    </td>
    <td style="text-align: left;">
    (required only when datatype = ‘incidence_raw’ and input data is
    matrix/data.frame) a vector of nonnegative integers specifying the
    number of sampling units in each assemblage. If assemblage names are
    not specified, then assemblages are automatically named as
    ‘assemblage1’, ‘assemblage2’,…, etc.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDtree
    </td>
    <td style="text-align: left;">
    (required argument for <code>diversity = ‘PD’</code>), a
    phylogenetic tree in Newick format for all observed species in the
    pooled assemblage.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDreftime
    </td>
    <td style="text-align: left;">
    (argument only for <code>diversity = ‘PD’</code>), a vector of
    numerical values specifying reference times for PD. Default is
    <code>NULL</code> (i.e., the age of the root of
    <code>PDtree</code>).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    PDtype
    </td>
    <td style="text-align: left;">
    (argument only for <code>diversity = ‘PD’</code>), select PD type:
    <code>PDtype = ‘PD’</code> (effective total branch length) or
    <code>PDtype = ‘meanPD’</code> (effective number of equally
    divergent lineages). Default is <code>‘meanPD’</code>, where
    <code>meanPD = PD/tree depth</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDdistM
    </td>
    <td style="text-align: left;">
    (required argument for <code>diversity = ‘PD’</code>), a species
    pairwise distance matrix for all species in the pooled assemblage.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDtype
    </td>
    <td style="text-align: left;">
    (argument only for <code>diversity = ‘FD’</code>), select FD type:
    <code>FDtype = ‘tau_values’</code> for FD under specified threshold
    values, or <code>FDtype = ‘AUC’</code> (area under the curve of
    tau-profile) for an overall FD which integrates all threshold values
    between zero and one. Default is <code>‘AUC’</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDtau
    </td>
    <td style="text-align: left;">
    (argument only for <code>diversity = ‘FD’</code> and <code>FDtype =
    ‘tau_values’</code>), a numerical vector between 0 and 1 specifying
    tau values (threshold levels). If <code>NULL</code> (default), then
    threshold is set to be the mean distance between any two individuals
    randomly selected from the pooled assemblage (i.e., quadratic
    entropy).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    FDcut_number
    </td>
    <td style="text-align: left;">
    (argument only for <code>diversity = ‘FD’</code> and <code>FDtype =
    ‘AUC’</code>), a numeric number to cut \[0, 1\] interval into
    equal-spaced sub-intervals to obtain the AUC value by integrating
    the tau-profile. Equivalently, the number of tau values that will be
    considered to compute the integrated AUC value. Default is
    <code>FDcut_number = 50</code>. A larger value can be set to obtain
    more accurate AUC value.
    </td>
    </tr>
    <tr>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    details
    </td>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    a logical variable to decide whether do you want to print out the
    detailed value for each plots, default is `FALSE`.
    </td>
    </tr>
    </tbody>
    </table>

Data type has following formats, such as a vector gathered by factors
(`abundance` and `incidence_freq`), a matrix/data frame with species
versus a assemblage (`abundance` `incidence_freq`, and `incidence_raw`),
a list of several vectors (`abundance` and `incidence_freq`), or a list
correspond to a assemblage (`incidence_raw`). `data` should comform the
format of each datatype. When `datatype = "incidence_raw"` and class of
data is matrix/data frame, user should input `nT` for each assemblage
which represents sampling units.

`diversity` contains three attributes diversity dimensions:
`Taxonomic diversity`, `Phylogenetic diversity`, and
`Functional diversity`. User should choose one `diversity`: `TD` means
`Taxonomic diversity`, `PD` means `Phylogenetic diversity` under a
specified reference time (default is root height), and `FD` means
`Functional diversity`. For the `"Functional diversity"`,
`FDtype = "AUC"` consider overall functional diversity which integrates
all threshold values between 0 to 1. And `FDtype = "tau_values"`
computes functional diversity under specified thresholds.

When `diversity = "PD"`, user should input `PDtree` Newick format data
for all observed species. When `diversity = "FD"`, user should input
`FDdistM` data matrix. Each element of the matrix is the pairwise
distance between any two observed species. And the species
identification names should be listed as row names and column names of
the matrix.

`nboot` is applied to get confidence interval, which is estimated by
bootstrap method. `details` means a logical setting whether print out
the computation value of all figures.

The output of `iNEXT4steps` will have three parts (if `details = TRUE`):
`$summary`, `$figure`, and `$details`. It may take some time to compute
when data size is large or `nboot` is large.

## <span style="color:blue;">TAXONOMIC DIVERSITY (TD): 4 STEPS VIA EXAMPLES</span>

`"abundance Data"` is used for a random sampling scheme. If the species
has aggregation effect, such as trees or plants, then set
`datatype = "incidence_freq"` or `datatype = "incidence_raw"`. First, we
use data `Spider` to compute taxonomic diversity.

### EXAMPLE 1: TD 4 steps for abundance data

Datasets `Spider` were sampled in a mountain forest ecosystem in the
Bavarian Forest National Park, Germany (Thorn et al. 2016, 2017). A
total of 12 experimental plots were established in “closed forest”
stands (6 plots) and “open forest” stands with naturally occurring gaps
and edges (6 plots) to assess the effects of microclimate on communities
of epigeal (ground-dwelling) spiders. Epigeal spiders were sampled over
three years with four pitfall traps in each plot, yielding a total of
3171 individuals belonging to 85 species recorded in the pooled habitat.
More details refer to data Source : A mountain forest ecosystem in the
Bavarian Forest National Park, Germany (Thorn et al. 2016, 2017).

Run the following code to get the output from iNEXT4steps.

``` r
data(Spider)
out1 <- iNEXT4steps(data = Spider, diversity = "TD", datatype = "abundance")
out1$summary
$`STEP1. Sample completeness profiles`
  Assemblage q = 0 q = 1 q = 2
1     Closed  0.61  0.99     1
2       Open  0.77  0.99     1

$`STEP2. Asymptotic analysis`
  Assemblage               qTD TD_obs TD_asy  s.e. qTD.LCL qTD.UCL
1     Closed  Species richness  44.00  72.11 18.27   36.30  107.91
2     Closed Shannon diversity  10.04  10.30  0.47    9.37   11.23
3     Closed Simpson diversity   5.71   5.73  0.29    5.16    6.30
4       Open  Species richness  74.00  96.31 11.26   74.25  118.37
5       Open Shannon diversity  16.34  16.84  0.56   15.74   17.93
6       Open Simpson diversity   9.41   9.46  0.36    8.75   10.17

$`STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.994 q = 0 q = 1 q = 2
1       Closed 55.62 10.18  5.72
2         Open 86.51 16.59  9.43

$`STEP4. Evenness among species abundances`
       Pielou J' q = 1 q = 2
Closed      0.58  0.17  0.09
Open        0.63  0.18  0.10
out1$figure[[6]]
```

<img src="README/README-unnamed-chunk-11-1.png" width="100%" style="display: block; margin: auto;" />

`$summary` lists all biological summaries according to Chao et
al. (2020). There are four parts corresponding to each step in the
paper. They analysis and explain biological data from different and
superimposed side. User can easily compare difference between each
assemblages.

`$figure` visualize the statistics by continuous curves. From the above
five figures, `iNEXT4stpes` provides a standard analysis process from
figure (a) to figure (e). User can analyze the process of biodiversity
through these figures.

`$details` contains four parts: `Sample Completeness`, `iNEXT`,
`Asymptotic Diversity`, `Evenness`. They are the computing values which
are used to plot each figure in `$figure`.

### EXAMPLE 2: TD 4 steps for incidence data

Incidence raw data is a species-by-incidence-sampling-units
matrix/data.frame. We split a space into several quadrats and only
record whether the species is detected or undetected in each quadrat.
According to this sampling scheme, `incidence_raw` data has only value
“zero” (undetected) or “one” (detected) in matrix/data frame (species by
assemblages). `incidence_freq` data is the total incidence frequency for
each species (i.e., row sums of the corresponding incidence raw matrix).
`incidence_freq` data should contain total sampling units (number of
quadrats) in the first row/entry.

*Note:* The phylogenetic diversity can only accept
`datatype = "incidence_raw"` for incidence-based data.

Datasets `Woody plants` are a subset of The National Vegetation Database
of Taiwan (AS-TW-001), sampled between 2003 and 2007 within the first
national vegetation inventory project (Chiou et al. 2009). Over 3600
vegetation plots, each 20x20-m in area, were set up in various locations
in Taiwan, and all woody plant individuals taller than 2 meters were
recorded in each plot. For illustration here, we selected only plots
belonging to two vegetation types (according to Li et al. 2013):
Pyrenaria-Machilus subtropical winter monsoon forest and Chamaecyparis
montane mixed cloud forest, sampled in the northern part of Taiwan (in
ecoregions 7 and 8 according to Su 1985).

Run the following code to get the output from iNEXT4steps.

``` r
data(Woody_plants)
out2 <- iNEXT4steps(data = Woody_plants[,c(1,4)], diversity = "TD", datatype = "incidence_freq")
out2$summary
$`STEP1. Sample completeness profiles`
   Assemblage q = 0 q = 1 q = 2
1     Monsoon  0.78  0.99     1
2 Upper_cloud  0.78  0.98     1

$`STEP2. Asymptotic analysis`
   Assemblage               qTD TD_obs TD_asy  s.e. qTD.LCL qTD.UCL
1     Monsoon  Species richness 329.00 421.67 27.03  368.69  474.65
2     Monsoon Shannon diversity 145.65 150.15  1.72  146.78  153.52
3     Monsoon Simpson diversity 102.33 103.35  1.37  100.67  106.03
4 Upper_cloud  Species richness 239.00 307.78 19.63  269.30  346.25
5 Upper_cloud Shannon diversity 105.53 110.50  1.99  106.60  114.39
6 Upper_cloud Simpson diversity  71.17  72.23  1.24   69.80   74.65

$`STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.993  q = 0  q = 1  q = 2
1      Monsoon 359.80 147.29 102.67
2  Upper_cloud 278.96 108.52  71.69

$`STEP4. Evenness among species abundances`
            Pielou J' q = 1 q = 2
Monsoon          0.85  0.41  0.28
Upper_cloud      0.83  0.39  0.25
out2$figure[[6]]
```

<img src="README/README-unnamed-chunk-12-1.png" width="100%" style="display: block; margin: auto;" />

## <span style="color:blue;">PHYLOGENETIC DIVERSITY (PD): 4 STEPS VIA EXAMPLES</span>

Here use abundance data: `"Brazil_abun_data"` as example.
`"Brazil_abun_data"` data has two main communities: “Edge”, “Interior”.
Here we also provide phylogenetic tree (`Brazil_phylo_tree`) data for
`"Brazil_abun_data"` data to compute phylogenetic diversity.

Run the following code to get the output from iNEXT4steps.

``` r
data(Brazil_abun_data)
data(Brazil_phylo_tree)
out3 <- iNEXT4steps(data = Brazil_abun_data, diversity = "PD", datatype = "abundance", nboot = 20, 
                    PDtree = Brazil_phylo_tree)
out3$summary
$`STEP1. Sample completeness profiles`
  Assemblage q = 0 q = 1 q = 2
1       Edge  0.72  0.94     1
2   Interior  0.69  0.94     1

$`STEP2. Asymptotic analysis`
  Assemblage      qPD PD_obs PD_asy s.e. qPD.LCL qPD.UCL Reftime   Type
1       Edge q = 0 PD  61.29  80.03 4.72   70.79   89.27     400 meanPD
2       Edge q = 1 PD   5.25   5.37 0.13    5.13    5.62     400 meanPD
3       Edge q = 2 PD   1.80   1.80 0.03    1.75    1.85     400 meanPD
4   Interior q = 0 PD  69.32  86.38 3.41   79.69   93.06     400 meanPD
5   Interior q = 1 PD   5.72   5.85 0.11    5.64    6.07     400 meanPD
6   Interior q = 2 PD   1.91   1.91 0.03    1.86    1.97     400 meanPD

$`STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.973 q = 0 q = 1 q = 2
1         Edge 71.76  5.32  1.80
2     Interior 80.32  5.80  1.91

$`STEP4. Evenness among species abundances`
         Pielou J' q = 1 q = 2
Edge          0.86  0.43  0.21
Interior      0.85  0.40  0.16
out3$figure[[6]]
```

<img src="README/README-unnamed-chunk-13-1.png" width="100%" style="display: block; margin: auto;" />

## <span style="color:blue;">FUNCTIONAL DIVERSITY (FD): 4 STEPS VIA EXAMPLES</span>

Here use abundance data `Brazil_abun_data` and its pairwise distance
matrix (`Brazil_distance_matrix`) to compute functional diversity. Under
`FDtype = "tau_values"`, user can key in `FDtau` as thresholds (default
is `Cmax`). Under `FDtype = "AUC"`, it considers overall functional
diversity which integrates all threshold values between 0 to 1. And
`FDtype = "tau_values"` computes functional diversity under specified
thresholds.

Run the following code to get the output from iNEXT4steps.

``` r
data(Brazil_abun_data)
data(Brazil_distance_matrix)
out4 <- iNEXT4steps(data = Brazil_abun_data, diversity = "FD", datatype = "abundance", nboot = 20, 
                    FDdistM = Brazil_distance_matrix, FDtype = 'tau_values')
out4$summary
$`STEP1. Sample completeness profiles`
  Assemblage q = 0 q = 1 q = 2
1       Edge  0.72  0.94     1
2   Interior  0.69  0.94     1

$`STEP2. Asymptotic analysis`
  Assemblage                  qFD FD_obs FD_asy s.e. qFD.LCL qFD.UCL  Tau
1       Edge q = 0 FD(single tau)   6.86   6.86 0.17    6.52    7.20 0.35
2       Edge q = 1 FD(single tau)   6.52   6.52 0.12    6.28    6.76 0.35
3       Edge q = 2 FD(single tau)   6.26   6.26 0.11    6.05    6.47 0.35
4   Interior q = 0 FD(single tau)   5.91   5.91 0.07    5.76    6.05 0.35
5   Interior q = 1 FD(single tau)   5.19   5.19 0.09    5.02    5.36 0.35
6   Interior q = 2 FD(single tau)   4.72   4.72 0.09    4.54    4.90 0.35

$`STEP3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.973 q = 0 q = 1 q = 2
1         Edge  6.86  6.53  6.26
2     Interior  5.91  5.19  4.72

$`STEP4. Evenness among species abundances`
         Pielou J' q = 1 q = 2
Edge          0.86  0.43  0.21
Interior      0.85  0.40  0.16
out4$figure[[6]]
```

<img src="README/README-unnamed-chunk-14-1.png" width="100%" style="display: block; margin: auto;" />

## <span style="color:red;">FUNCTION Completeness: SAMPLE COMPLETENESS PROFILES</span>

`iNEXT.4steps` provides function `Completeness()` to compute estimated
sample completeness with order q. The arguments is below:

``` r
Completeness(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 50, conf = 0.95, nT = NULL)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style='text-align: left;'>

1.  For `datatype = 'abundance'`, data can be input as a vector of
    species abundances (for a single assemblage), matrix/data.frame
    (species by assemblages), or a list of species abundance
    vectors. (2) For `datatype = 'incidence_freq'`, data can be input as
    a vector of incidence frequencies (for a single assemblage),
    matrix/data.frame (species by assemblages), or a list of incidence
    frequencies; the first entry in all types of input must be the
    number of sampling units in each assemblage. (3) For
    `datatype = 'incidence_raw'`, data can be input as a list of
    matrix/data.frame (species by sampling units); data can also be
    input as a matrix/data.frame by merging all sampling units across
    assemblages based on species identity; in this case, the number of
    sampling units (nT, see below) must be input.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    q
    </td>
    <td style="text-align: left;">
    a numerical vector specifying the diversity orders for q-profile
    output. Default is <code>seq(0, 2, by = 0.2)</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    datatype
    </td>
    <td style="text-align: left;">
    data type of input data: individual-based abundance data (datatype =
    ‘abundance’), sampling-unit-based incidence frequencies data
    (datatype = ‘incidence_freq’), or species by sampling-units
    incidence matrix (datatype = ‘incidence_raw’) with all entries being
    0 (non-detection) or 1 (detection).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nboot
    </td>
    <td style="text-align: left;">
    a positive integer specifying the number of bootstrap replications
    when assessing sampling uncertainty and constructing confidence
    intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    conf
    </td>
    <td style="text-align: left;">
    a positive number \< 1 specifying the level of confidence interval.
    Default is 0.95.
    </td>
    </tr>
    <tr>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    nT
    </td>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    (required only when datatype = ‘incidence_raw’ and input data is
    matrix/data.frame) a vector of nonnegative integers specifying the
    number of sampling units in each assemblage. If assemblage names are
    not specified, then assemblages are automatically named as
    ‘assemblage1’, ‘assemblage2’,…, etc.
    </td>
    </tr>
    </tbody>
    </table>

## <span style="color:red;">FUNCTION ggCompleteness(): GRAPHIC DISPLAYS OF SAMPLE COMPLETENESS PROFILES</span>

`iNEXT.4steps` also provides a visualized function `ggCompleteness` to
plot the output from `Completeness()`:

``` r
ggCompleteness(output)
```

Following are two simple examples for functions `Completeness` and
`ggCompleteness`.

### Sample completeness profiles for abundance data

Use abundance data `Spider` to compute estimated sample completeness.

``` r
data(Spider)
out1 <- Completeness(data = Spider, datatype = "abundance")
out1
   Order.q Estimate.SC         s.e.    SC.LCL    SC.UCL Assemblage
1      0.0   0.7683622 9.481347e-02 0.5825312 0.9541932       Open
2      0.2   0.8181065 7.064529e-02 0.6796443 0.9565687       Open
3      0.4   0.8768761 4.269957e-02 0.7931865 0.9605658       Open
4      0.6   0.9302044 2.004604e-02 0.8909149 0.9694940       Open
5      0.8   0.9664828 7.473631e-03 0.9518347 0.9811308       Open
6      1.0   0.9858045 2.378624e-03 0.9811425 0.9904665       Open
7      1.2   0.9944729 7.210652e-04 0.9930596 0.9958861       Open
8      1.4   0.9979624 2.443664e-04 0.9974834 0.9984413       Open
9      1.6   0.9992757 9.868063e-05 0.9990823 0.9994691       Open
10     1.8   0.9997489 4.254166e-05 0.9996656 0.9998323       Open
11     2.0   0.9999146 1.804458e-05 0.9998792 0.9999499       Open
12     0.0   0.6102206 1.290200e-01 0.3573462 0.8630951     Closed
13     0.2   0.7182647 9.851865e-02 0.5251717 0.9113577     Closed
14     0.4   0.8338575 5.765765e-02 0.7208506 0.9468644     Closed
15     0.6   0.9215916 2.330820e-02 0.8759084 0.9672748     Closed
16     0.8   0.9692426 6.962278e-03 0.9555968 0.9828884     Closed
17     1.0   0.9893733 1.873850e-03 0.9857006 0.9930459     Closed
18     1.2   0.9966174 5.822561e-04 0.9954762 0.9977586     Closed
19     1.4   0.9989804 2.203640e-04 0.9985485 0.9994123     Closed
20     1.6   0.9997042 8.700597e-05 0.9995337 0.9998748     Closed
21     1.8   0.9999166 3.326519e-05 0.9998514 0.9999818     Closed
22     2.0   0.9999770 1.223136e-05 0.9999530 1.0000000     Closed
```

The following commands plot sample completeness curves, along with its
confidence interval for q between 0 to 2.

``` r
ggCompleteness(out1)
```

<img src="README/README-unnamed-chunk-19-1.png" width="70%" style="display: block; margin: auto;" />

### Sample completeness profiles for incidence data

Use incidence frequency data `woody plants` to compute sample
completeness.

``` r
data(Woody_plants)
out2 <- Completeness(data = Woody_plants[,c(1,4)], datatype = "incidence_freq")
out2
   Order.q Estimate.SC         s.e.    SC.LCL    SC.UCL  Assemblage
1      0.0   0.7765330 3.945082e-02 0.6992108 0.8538552 Upper_cloud
2      0.2   0.8358241 2.730851e-02 0.7823004 0.8893478 Upper_cloud
3      0.4   0.8915125 1.626604e-02 0.8596316 0.9233933 Upper_cloud
4      0.6   0.9354738 8.335205e-03 0.9191371 0.9518105 Upper_cloud
5      0.8   0.9649731 3.755301e-03 0.9576129 0.9723334 Upper_cloud
6      1.0   0.9823020 1.554049e-03 0.9792562 0.9853479 Upper_cloud
7      1.2   0.9915212 6.445519e-04 0.9902579 0.9927845 Upper_cloud
8      1.4   0.9960932 3.032244e-04 0.9954989 0.9966875 Upper_cloud
9      1.6   0.9982514 1.641998e-04 0.9979296 0.9985733 Upper_cloud
10     1.8   0.9992348 9.264975e-05 0.9990532 0.9994164 Upper_cloud
11     2.0   0.9996711 5.117167e-05 0.9995708 0.9997714 Upper_cloud
12     0.0   0.7802244 4.946689e-02 0.6832711 0.8771778     Monsoon
13     0.2   0.8490462 3.254597e-02 0.7852573 0.9128352     Monsoon
14     0.4   0.9086970 1.762230e-02 0.8741580 0.9432361     Monsoon
15     0.6   0.9508262 7.923499e-03 0.9352965 0.9663560     Monsoon
16     0.8   0.9758389 3.049352e-03 0.9698623 0.9818155     Monsoon
17     1.0   0.9888942 1.039192e-03 0.9868574 0.9909310     Monsoon
18     1.2   0.9951304 3.404490e-04 0.9944631 0.9957977     Monsoon
19     1.4   0.9979362 1.379567e-04 0.9976658 0.9982066     Monsoon
20     1.6   0.9991474 7.495364e-05 0.9990005 0.9992943     Monsoon
21     1.8   0.9996547 4.242513e-05 0.9995716 0.9997379     Monsoon
22     2.0   0.9998625 2.270661e-05 0.9998179 0.9999070     Monsoon
```

The following commands plot sample completeness curves, along with its
confidence interval for q between 0 to 2.

``` r
ggCompleteness(out2)
```

<img src="README/README-unnamed-chunk-21-1.png" width="70%" style="display: block; margin: auto;" />

## <span style="color:red;">FUNCTION Evenness: EVENNESS PROFILES</span>

`iNEXT.4steps` provides the function `Evenness()` to compute observed
eveness or estimated evenness under specified sample coverage. The
argument is below:

``` r
Evenness(data, q = seq(0, 2, 0.2), datatype = "abundance", method = "Estimated",
         nboot = 50, conf = 0.95, nT = NULL, E.class = 1:5, SC = NULL)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Argument
</th>
<th style="font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;">
Description
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;">
data
</td>
<td style='text-align: left;'>

1.  For `datatype = 'abundance'`, data can be input as a vector of
    species abundances (for a single assemblage), matrix/data.frame
    (species by assemblages), or a list of species abundance
    vectors. (2) For `datatype = 'incidence_freq'`, data can be input as
    a vector of incidence frequencies (for a single assemblage),
    matrix/data.frame (species by assemblages), or a list of incidence
    frequencies; the first entry in all types of input must be the
    number of sampling units in each assemblage. (3) For
    `datatype = 'incidence_raw'`, data can be input as a list of
    matrix/data.frame (species by sampling units); data can also be
    input as a matrix/data.frame by merging all sampling units across
    assemblages based on species identity; in this case, the number of
    sampling units (nT, see below) must be input.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    q
    </td>
    <td style="text-align: left;">
    a numerical vector specifying the diversity orders for q-profile
    output. Default is <code>seq(0, 2, by = 0.2)</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    datatype
    </td>
    <td style="text-align: left;">
    data type of input data: individual-based abundance data (datatype =
    ‘abundance’), sampling-unit-based incidence frequencies data
    (datatype = ‘incidence_freq’), or species by sampling-units
    incidence matrix (datatype = ‘incidence_raw’) with all entries being
    0 (non-detection) or 1 (detection).
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    method
    </td>
    <td style="text-align: left;">
    a binary calculation method with <code>‘Estimated’</code> or
    <code>‘Observed’</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nboot
    </td>
    <td style="text-align: left;">
    a positive integer specifying the number of bootstrap replications
    when assessing sampling uncertainty and constructing confidence
    intervals. Enter 0 to skip the bootstrap procedures. Default is 50.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    conf
    </td>
    <td style="text-align: left;">
    a positive number \< 1 specifying the level of confidence interval.
    Default is 0.95.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nT
    </td>
    <td style="text-align: left;">
    (required only when datatype = ‘incidence_raw’ and input data is
    matrix/data.frame) a vector of nonnegative integers specifying the
    number of sampling units in each assemblage. If assemblage names are
    not specified, then assemblages are automatically named as
    ‘assemblage1’, ‘assemblage2’,…, etc.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    E.class
    </td>
    <td style="text-align: left;">
    an integer vector between 1 to 5.
    </td>
    </tr>
    <tr>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    SC
    </td>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    (required only when method = ‘Estimated’) a standardized coverage
    for calculating estimated evenness. If <code>NULL</code>, then this
    function computes the diversity estimates for the minimum sample
    coverage among all samples extrapolated to double reference sizes
    (Cmax).
    </td>
    </tr>
    </tbody>
    </table>

## <span style="color:red;">FUNCTION ggEvenness(): GRAPHIC DISPLAYS OF EVENNESS PROFILES</span>

`iNEXT.4steps` provide a function `ggEvenness()` to plot the output from
`Evenness()`.

``` r
ggEvenness(output)
```

Following are two simple examples for functions `Evenness` and
`ggEvenness`.

### Evenness profiles for abundance data

Use abundance data `Spider` to compute estimated evenness under
`C = Cmax`. Here only show the output for first class of evenness.

``` r
data(Spider)
out1 <- Evenness(data = Spider, datatype = "abundance", 
                 method = "Estimated", SC = NULL, E.class = 1:5)
out1$E1
   Order.q  Evenness        s.e.  Even.LCL  Even.UCL Assemblage    Method        SC
1      0.0 1.0000000 0.000000000 1.0000000 1.0000000       Open Estimated 0.9937676
2      0.2 0.7276194 0.017474610 0.6933698 0.7618690       Open Estimated 0.9937676
3      0.4 0.6008010 0.022861236 0.5559938 0.6456082       Open Estimated 0.9937676
4      0.6 0.5628529 0.022811667 0.5181429 0.6075630       Open Estimated 0.9937676
5      0.8 0.5801167 0.020159768 0.5406043 0.6196292       Open Estimated 0.9937676
6      1.0 0.6297250 0.016344703 0.5976899 0.6617600       Open Estimated 0.9937676
7      1.2 0.6941362 0.012378998 0.6698738 0.7183985       Open Estimated 0.9937676
8      1.4 0.7600773 0.008997610 0.7424423 0.7777123       Open Estimated 0.9937676
9      1.6 0.8190763 0.006528693 0.8062803 0.8318723       Open Estimated 0.9937676
10     1.8 0.8673381 0.004889999 0.8577538 0.8769223       Open Estimated 0.9937676
11     2.0 0.9044427 0.003809103 0.8969770 0.9119084       Open Estimated 0.9937676
12     0.0 1.0000000 0.000000000 1.0000000 1.0000000     Closed Estimated 0.9937676
13     0.2 0.6996114 0.034852544 0.6313017 0.7679211     Closed Estimated 0.9937676
14     0.4 0.5665175 0.045837792 0.4766771 0.6563579     Closed Estimated 0.9937676
15     0.6 0.5253117 0.046456519 0.4342586 0.6163648     Closed Estimated 0.9937676
16     0.8 0.5365988 0.041914965 0.4544469 0.6187506     Closed Estimated 0.9937676
17     1.0 0.5775431 0.034823871 0.5092895 0.6457966     Closed Estimated 0.9937676
18     1.2 0.6330777 0.027063319 0.5800346 0.6861209     Closed Estimated 0.9937676
19     1.4 0.6926130 0.020057779 0.6533004 0.7319255     Closed Estimated 0.9937676
20     1.6 0.7491390 0.014616559 0.7204911 0.7777870     Closed Estimated 0.9937676
21     1.8 0.7988183 0.010883244 0.7774875 0.8201490     Closed Estimated 0.9937676
22     2.0 0.8402207 0.008545150 0.8234725 0.8569689     Closed Estimated 0.9937676
```

The following commands plot the evenness curves, along with its
confidence interval for q between 0 to 2.

``` r
ggEvenness(out1)
```

<img src="README/README-unnamed-chunk-26-1.png" width="100%" style="display: block; margin: auto;" />

### Evenness profiles for incidence data

Use incidence frequency data `Woody plants` to calculate estimated
evenness under `C = NULL` (Cmax). Here only show the output for first
class of evenness.

``` r
data(Woody_plants)
out2 <- Evenness(data = Woody_plants[,c(1,4)], datatype = "incidence_freq", 
                 method = "Estimated", SC = NULL, E.class = 1:5)
out2$E1
   Order.q  Evenness         s.e.  Even.LCL  Even.UCL  Assemblage    Method        SC
1      0.0 1.0000000 0.0000000000 1.0000000 1.0000000 Upper_cloud Estimated 0.9925846
2      0.2 0.8302829 0.0076490280 0.8152911 0.8452747 Upper_cloud Estimated 0.9925846
3      0.4 0.7565487 0.0100862627 0.7367800 0.7763174 Upper_cloud Estimated 0.9925846
4      0.6 0.7473532 0.0098008078 0.7281440 0.7665625 Upper_cloud Estimated 0.9925846
5      0.8 0.7795871 0.0080657630 0.7637785 0.7953957 Upper_cloud Estimated 0.9925846
6      1.0 0.8323270 0.0057570559 0.8210434 0.8436106 Upper_cloud Estimated 0.9925846
7      1.2 0.8869413 0.0035944817 0.8798962 0.8939863 Upper_cloud Estimated 0.9925846
8      1.4 0.9312253 0.0020038459 0.9272979 0.9351528 Upper_cloud Estimated 0.9925846
9      1.6 0.9614213 0.0010339467 0.9593948 0.9634478 Upper_cloud Estimated 0.9925846
10     1.8 0.9795937 0.0005149306 0.9785845 0.9806030 Upper_cloud Estimated 0.9925846
11     2.0 0.9895988 0.0002596511 0.9890899 0.9901077 Upper_cloud Estimated 0.9925846
12     0.0 1.0000000 0.0000000000 1.0000000 1.0000000     Monsoon Estimated 0.9925846
13     0.2 0.8270279 0.0083714125 0.8106202 0.8434355     Monsoon Estimated 0.9925846
14     0.4 0.7581724 0.0108983810 0.7368120 0.7795328     Monsoon Estimated 0.9925846
15     0.6 0.7554881 0.0104049651 0.7350947 0.7758814     Monsoon Estimated 0.9925846
16     0.8 0.7929638 0.0083519563 0.7765943 0.8093333     Monsoon Estimated 0.9925846
17     1.0 0.8482450 0.0057547836 0.8369658 0.8595241     Monsoon Estimated 0.9925846
18     1.2 0.9022994 0.0034232908 0.8955899 0.9090089     Monsoon Estimated 0.9925846
19     1.4 0.9437403 0.0017898253 0.9402323 0.9472483     Monsoon Estimated 0.9925846
20     1.6 0.9703214 0.0008497857 0.9686558 0.9719869     Monsoon Estimated 0.9925846
21     1.8 0.9853020 0.0003809243 0.9845554 0.9860486     Monsoon Estimated 0.9925846
22     2.0 0.9930196 0.0001685398 0.9926893 0.9933500     Monsoon Estimated 0.9925846
```

The following commands plot the evenness curves, along with its
confidence interval for q between 0 to 2.

``` r
ggEvenness(out2)
```

<img src="README/README-unnamed-chunk-28-1.png" width="100%" style="display: block; margin: auto;" />

## License

The iNEXT.4steps package is licensed under the GPLv3. To help refine
`iNEXT.4steps`, your comments or feedback would be welcome (please send
them to Anne Chao or report an issue on the iNEXT.4steps github
[iNEXT.4steps_github](https://github.com/AnneChao/iNEXT.4steps).

## References

-   Chao, A., Gotelli, N. G., Hsieh, T. C., Sander, E. L., Ma, K. H.,
    Colwell, R. K. and Ellison, A. M. (2014). Rarefaction and
    extrapolation with Hill numbers: a framework for sampling and
    estimation in species biodiversity studies. Ecological Monographs
    84, 45-67.

-   Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H.,
    Dornelas, M and Magurran, A. E. (2021). Measuring temporal change in
    alpha diversity: a framework integrating taxonomic, phylogenetic and
    functional diversity and the iNEXT.3D standardization. Methods in
    Ecology and Evolution, 12, 1926-1940.

-   Chao, A. and Jost. L. (2012) Coverage-based rarefaction and
    extrapolation: standardizing samples by completeness rather than
    size. Ecology, 93, 2533-2547.

-   Chao, A. and Jost, L. (2015). Estimating diversity and entropy
    profiles via discovery rates of new species. Methods in Ecology and
    Evolution, 6, 873-882.

-   Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it
    to diversity, beta diversity, and similarity. Ecology, 100(12),
    e02852.

-   Chao, A., Y. Kubota, D. Zeleny, C.-H. Chiu, C.-F. Li, B.
    Kusumoto, M. Yasuhara, S. Thorn, C.-L. Wei, M. J. Costello,
    and R. K. Colwell (2020). Quantifying sample completeness and
    comparing diversities among assemblages. Ecological Research, 35,
    292-314.

-   Chiou, C.-R., Hsieh, C.-F., Wang, J.-C., Chen, M.-Y., Liu, H.-Y.,
    Yeh, C.-L., … Song, M. G.-Z. (2009). The first national vegetation
    inventory in Taiwan. Taiwan Journal of Forest Science, 24, 295–302.

-   Li, C.-F., Chytry, M., Zeleny, D., Chen, M. -Y., Chen, T.-Y., Chiou,
    C.-R., … Hsieh, C.-F. (2013). Classification of Taiwan forest
    vegetation. Applied Vegetation Science, 16, 698–719.
    <https://doi.org/10.1111/avsc.12025>

-   Su, H. -J. (1985). Studies on the climate and vegetation types of
    the natural forests in Taiwan (III) A Scheme of Geographical
    Climatic Regions. Quarterly Journal of Chinese Forestry, 18, 33–44.

-   Thorn, S., Bassler, C., Svoboda, M., & Müller, J. (2017). Effects of
    natural disturbances and salvage logging on biodiversity - lessons
    from the bohemian Forest. Forest Ecology and Management, 388,
    113–119. <https://doi.org/10.1016/j.foreco.2016.06.006>

-   Thorn, S., BuBler, H., Fritze, M. -A., Goeder, P., Muller, J., WeiB,
    I., & Seibold, S. (2016). Canopy closure determines arthropod
    assemblages in microhabitats created by windstorms and salvage
    logging. Forest Ecology and Managemen
