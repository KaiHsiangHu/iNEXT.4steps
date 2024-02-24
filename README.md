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

-   **Step 1**: `Sample completeness profile`

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

-   **Step 2.1 and step 3**:
    `Size-based and coverage-based Interpolation and Extrapolation`

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

-   **Step 2.2**: `Asymptotic diversity profile`

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

-   **Step 4**: `Evenness profile`

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
Run the following code to view the first 16 rows and first 3 columns of
the incidence frequency data.

``` r
data("Woody_plants")
Woody_plants
```

     [1] 129   0   0   4   2   1   0   0  16   0   4   0   0   0   0   0

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
1     Closed  Species richness  44.00  72.11 17.27   38.25  105.96
2     Closed Shannon diversity  10.04  10.30  0.36    9.60   11.01
3     Closed Simpson diversity   5.71   5.73  0.22    5.30    6.15
4       Open  Species richness  74.00  96.31 15.39   66.14  126.47
5       Open Shannon diversity  16.34  16.84  0.60   15.66   18.01
6       Open Simpson diversity   9.41   9.46  0.33    8.81   10.11

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
1     Monsoon  Species richness 329.00 421.67 27.01  368.73  474.62
2     Monsoon Shannon diversity 145.65 150.15  1.71  146.79  153.51
3     Monsoon Simpson diversity 102.33 103.35  1.24  100.92  105.78
4 Upper_cloud  Species richness 239.00 307.78 21.03  266.55  349.01
5 Upper_cloud Shannon diversity 105.53 110.50  1.56  107.43  113.57
6 Upper_cloud Simpson diversity  71.17  72.23  1.02   70.22   74.23

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
1       Edge q = 0 PD  61.29  80.03 7.20   65.92   94.14     400 meanPD
2       Edge q = 1 PD   5.25   5.37 0.11    5.16    5.58     400 meanPD
3       Edge q = 2 PD   1.80   1.80 0.03    1.74    1.85     400 meanPD
4   Interior q = 0 PD  69.32  86.38 4.87   76.83   95.92     400 meanPD
5   Interior q = 1 PD   5.72   5.85 0.10    5.66    6.05     400 meanPD
6   Interior q = 2 PD   1.91   1.91 0.02    1.87    1.96     400 meanPD

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
1       Edge q = 0 FD(single tau)   6.86   6.86 0.22    6.42    7.30 0.35
2       Edge q = 1 FD(single tau)   6.52   6.52 0.14    6.26    6.79 0.35
3       Edge q = 2 FD(single tau)   6.26   6.26 0.11    6.05    6.47 0.35
4   Interior q = 0 FD(single tau)   5.91   5.91 0.06    5.79    6.03 0.35
5   Interior q = 1 FD(single tau)   5.19   5.19 0.08    5.04    5.34 0.35
6   Interior q = 2 FD(single tau)   4.72   4.72 0.08    4.55    4.88 0.35

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
1      0.0   0.7683622 7.881909e-02 0.6138796 0.9228448       Open
2      0.2   0.8181065 5.830653e-02 0.7038278 0.9323852       Open
3      0.4   0.8768761 3.520988e-02 0.8078660 0.9458862       Open
4      0.6   0.9302044 1.661199e-02 0.8976455 0.9627633       Open
5      0.8   0.9664828 6.246960e-03 0.9542389 0.9787266       Open
6      1.0   0.9858045 2.001376e-03 0.9818819 0.9897271       Open
7      1.2   0.9944729 6.077648e-04 0.9932817 0.9956641       Open
8      1.4   0.9979624 2.085246e-04 0.9975537 0.9983711       Open
9      1.6   0.9992757 8.695740e-05 0.9991052 0.9994461       Open
10     1.8   0.9997489 3.871000e-05 0.9996731 0.9998248       Open
11     2.0   0.9999146 1.684662e-05 0.9998815 0.9999476       Open
12     0.0   0.6102206 1.096237e-01 0.3953622 0.8250790     Closed
13     0.2   0.7182647 8.074973e-02 0.5599982 0.8765313     Closed
14     0.4   0.8338575 4.592607e-02 0.7438441 0.9238709     Closed
15     0.6   0.9215916 1.929979e-02 0.8837647 0.9594185     Closed
16     0.8   0.9692426 6.264917e-03 0.9569636 0.9815216     Closed
17     1.0   0.9893733 1.730946e-03 0.9859807 0.9927659     Closed
18     1.2   0.9966174 4.813958e-04 0.9956739 0.9975609     Closed
19     1.4   0.9989804 1.650075e-04 0.9986570 0.9993038     Closed
20     1.6   0.9997042 6.534401e-05 0.9995762 0.9998323     Closed
21     1.8   0.9999166 2.567485e-05 0.9998663 0.9999669     Closed
22     2.0   0.9999770 9.604946e-06 0.9999582 0.9999958     Closed
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
1      0.0   0.7765330 5.585889e-02 0.6670516 0.8860144 Upper_cloud
2      0.2   0.8358241 3.923805e-02 0.7589189 0.9127293 Upper_cloud
3      0.4   0.8915125 2.360559e-02 0.8452464 0.9377786 Upper_cloud
4      0.6   0.9354738 1.208284e-02 0.9117918 0.9591557 Upper_cloud
5      0.8   0.9649731 5.329743e-03 0.9545270 0.9754192 Upper_cloud
6      1.0   0.9823020 2.076560e-03 0.9782321 0.9863720 Upper_cloud
7      1.2   0.9915212 7.615384e-04 0.9900286 0.9930138 Upper_cloud
8      1.4   0.9960932 3.214034e-04 0.9954633 0.9967232 Upper_cloud
9      1.6   0.9982514 1.802434e-04 0.9978982 0.9986047 Upper_cloud
10     1.8   0.9992348 1.100643e-04 0.9990190 0.9994505 Upper_cloud
11     2.0   0.9996711 6.447329e-05 0.9995447 0.9997975 Upper_cloud
12     0.0   0.7802244 4.147521e-02 0.6989345 0.8615144     Monsoon
13     0.2   0.8490462 2.680152e-02 0.7965162 0.9015762     Monsoon
14     0.4   0.9086970 1.435193e-02 0.8805678 0.9368263     Monsoon
15     0.6   0.9508262 6.415580e-03 0.9382519 0.9634006     Monsoon
16     0.8   0.9758389 2.455199e-03 0.9710268 0.9806510     Monsoon
17     1.0   0.9888942 8.351815e-04 0.9872573 0.9905311     Monsoon
18     1.2   0.9951304 2.858212e-04 0.9945702 0.9956906     Monsoon
19     1.4   0.9979362 1.296975e-04 0.9976820 0.9981904     Monsoon
20     1.6   0.9991474 7.363881e-05 0.9990031 0.9992917     Monsoon
21     1.8   0.9996547 4.140270e-05 0.9995736 0.9997359     Monsoon
22     2.0   0.9998625 2.185975e-05 0.9998196 0.9999053     Monsoon
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
2      0.2 0.7276194 0.019876494 0.6886622 0.7665766       Open Estimated 0.9937676
3      0.4 0.6008010 0.025913540 0.5500114 0.6515906       Open Estimated 0.9937676
4      0.6 0.5628529 0.025786372 0.5123126 0.6133933       Open Estimated 0.9937676
5      0.8 0.5801167 0.022667896 0.5356885 0.6245450       Open Estimated 0.9937676
6      1.0 0.6297250 0.018166007 0.5941202 0.6653297       Open Estimated 0.9937676
7      1.2 0.6941362 0.013468596 0.6677382 0.7205341       Open Estimated 0.9937676
8      1.4 0.7600773 0.009485700 0.7414856 0.7786689       Open Estimated 0.9937676
9      1.6 0.8190763 0.006654641 0.8060335 0.8321192       Open Estimated 0.9937676
10     1.8 0.8673381 0.004889330 0.8577551 0.8769210       Open Estimated 0.9937676
11     2.0 0.9044427 0.003818402 0.8969588 0.9119266       Open Estimated 0.9937676
12     0.0 1.0000000 0.000000000 1.0000000 1.0000000     Closed Estimated 0.9937676
13     0.2 0.6996114 0.035731801 0.6295783 0.7696444     Closed Estimated 0.9937676
14     0.4 0.5665175 0.045811273 0.4767291 0.6563059     Closed Estimated 0.9937676
15     0.6 0.5253117 0.045432919 0.4362649 0.6143586     Closed Estimated 0.9937676
16     0.8 0.5365988 0.039999708 0.4582008 0.6149967     Closed Estimated 0.9937676
17     1.0 0.5775431 0.032126724 0.5145758 0.6405103     Closed Estimated 0.9937676
18     1.2 0.6330777 0.023770021 0.5864893 0.6796661     Closed Estimated 0.9937676
19     1.4 0.6926130 0.016464519 0.6603431 0.7248828     Closed Estimated 0.9937676
20     1.6 0.7491390 0.011103175 0.7273772 0.7709008     Closed Estimated 0.9937676
21     1.8 0.7988183 0.007836255 0.7834595 0.8141770     Closed Estimated 0.9937676
22     2.0 0.8402207 0.006186862 0.8280947 0.8523467     Closed Estimated 0.9937676
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
2      0.2 0.8302829 0.0112153534 0.8083012 0.8522646 Upper_cloud Estimated 0.9925846
3      0.4 0.7565487 0.0149289533 0.7272885 0.7858089 Upper_cloud Estimated 0.9925846
4      0.6 0.7473532 0.0146001205 0.7187375 0.7759689 Upper_cloud Estimated 0.9925846
5      0.8 0.7795871 0.0120660777 0.7559380 0.8032361 Upper_cloud Estimated 0.9925846
6      1.0 0.8323270 0.0086369347 0.8153989 0.8492551 Upper_cloud Estimated 0.9925846
7      1.2 0.8869413 0.0054001960 0.8763571 0.8975255 Upper_cloud Estimated 0.9925846
8      1.4 0.9312253 0.0030042462 0.9253371 0.9371135 Upper_cloud Estimated 0.9925846
9      1.6 0.9614213 0.0015352568 0.9584122 0.9644303 Upper_cloud Estimated 0.9925846
10     1.8 0.9795937 0.0007486165 0.9781265 0.9810610 Upper_cloud Estimated 0.9925846
11     2.0 0.9895988 0.0003673000 0.9888789 0.9903187 Upper_cloud Estimated 0.9925846
12     0.0 1.0000000 0.0000000000 1.0000000 1.0000000     Monsoon Estimated 0.9925846
13     0.2 0.8270279 0.0076794775 0.8119764 0.8420794     Monsoon Estimated 0.9925846
14     0.4 0.7581724 0.0099141578 0.7387410 0.7776038     Monsoon Estimated 0.9925846
15     0.6 0.7554881 0.0094102379 0.7370443 0.7739318     Monsoon Estimated 0.9925846
16     0.8 0.7929638 0.0075350196 0.7781954 0.8077322     Monsoon Estimated 0.9925846
17     1.0 0.8482450 0.0052010173 0.8380512 0.8584388     Monsoon Estimated 0.9925846
18     1.2 0.9022994 0.0031161898 0.8961918 0.9084070     Monsoon Estimated 0.9925846
19     1.4 0.9437403 0.0016532387 0.9405000 0.9469806     Monsoon Estimated 0.9925846
20     1.6 0.9703214 0.0008043787 0.9687448 0.9718979     Monsoon Estimated 0.9925846
21     1.8 0.9853020 0.0003736373 0.9845697 0.9860343     Monsoon Estimated 0.9925846
22     2.0 0.9930196 0.0001728704 0.9926808 0.9933584     Monsoon Estimated 0.9925846
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
