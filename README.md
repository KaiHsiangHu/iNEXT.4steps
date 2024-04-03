<!-- README.md is generated from README.Rmd. Please edit that file -->

# iNEXT.4steps (R package)

<h5 align="right">
Latest version: 2024-04-03
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

<br> `iNEXT.4steps` (Four-Step Biodiversity Analysis based on iNEXT)
expands `iNEXT` (Chao et al. 2014) to include the estimation of sample
completeness and evenness under a unified framework of Hill numbers.
`iNEXT.4steps` links sample completeness, diversity estimation,
interpolation and extrapolation (`iNEXT`), and evenness in a fully
integrated approach. The pertinent background for the four-step
methodology is provided in Chao et al. (2020). The four-step procedures
are described in the following:

-   **Step 1: Assessment of sample completeness profile**

Before performing biodiversity analysis, it is important to first
quantify the sample completeness of a biological survey. Chao et
al. (2020) generalized the conventional sample completeness to a class
of measures parametrized by an order q ≥ 0. When q = 0, sample
completeness reduces to the conventional measure of completeness, i.e.,
the ratio of the observed species richness to the true richness
(observed plus undetected). When q = 1, the measure reduces to the
sample coverage (the proportion of the total number of individuals in
the entire assemblage that belong to detected species), a concept
original developed by Alan Turing in his cryptographic analysis during
WWII. When q = 2, it represents a generalized sample coverage with each
species being proportionally weighted by its squared species abundance
(i.e., each individual being proportionally weighted by its species
abundance); this measure thus is disproportionally sensitive to highly
abundant species. For a general order q ≥ 0 (not necessarily to be an
integer) , the sample completeness of order q quantifies the proportion
of the assemblage’s individuals belonging to detected species, with each
individual being proportionally weighted by the (q-1)th power of its
abundance. Sample completeness profile depicts its estimate with respect
to order q ≥ 0; this profile fully characterizes the sample completeness
of a biological survey.

`iNEXT.4steps` features the estimated profile for all orders of q ≥ 0
based on the methodology developed in Chao et al. (2020). All estimates
are theoretically between 0 and 1. If the estimated sample completeness
profile is a horizontal line at the level of unity for all orders of q ≥
0, then the survey is complete, implying there is no undetected
diversity. In most applications, the estimated profile increases with
order q, revealing the existence of undetected diversity. The sample
completeness estimate for q = 0 provides an upper bound for the
proportion of observed species; its complement represents a lower bound
for the proportion of undetected species. This interpretation is mainly
because data typically do not contain sufficient information to
accurately estimate species richness and only a lower bound of species
richness can be well estimated. By contrast, for q ≥ 1, when data are
not sparse, the sample completeness value for q ≥ 1 can be very
accurately estimated measures. The values for q ≥ 2 typically are very
close to unity, signifying that almost all highly abundant species (for
abundance data) or highly frequent species (for incidence data) had been
detected in the reference sample.

-   **STEP 2. Analysis of the size-based rarefaction and extrapolation
    sampling curves, and the asymptotic diversity profile for 0 ≤ q ≤
    2.**

**(STEP 2a)**. For each dataset, first examine the pattern of the
size-based rarefaction and extrapolation sampling curve up to double the
reference sample size for q = 0, 1 and 2. If the curve stays at a fixed
level (this often occurs for the measures of q = 1 and 2), then our
asymptotic estimate presented in Step 2b can be used to accurately infer
the true diversity of the entire assemblage. Otherwise, our asymptotic
diversity estimate represents only a lower bound (this often occurs for
the measures of q = 0).

**(STEP 2b)**. When the true diversity can be accurately inferred, the
extent of undetected diversity within each dataset is obtained by
comparing the estimated asymptotic diversity profile and empirical
profile; the difference in diversity between any two assemblages can be
evaluated and tested for significance.

-   **STEP 3. Analysis of non-asymptotic coverage-based rarefaction and
    extrapolation analysis for orders q = 0, 1 and 2.**

When sampling data do not contain sufficient information to accurately
infer true diversity, fair comparisons of diversity across multiple
assemblages should be made by standardizing the sample coverage (i.e.,
comparing diversity for a standardized fraction of an assemblage’s
individuals). This comparison can be done based on seamless integration
of coverage-based rarefaction and extrapolation sampling curves up to a
maximum coverage (Cmax = the minimum sample coverage among all samples
extrapolated to double reference sizes).

-   **STEP 4. Assessment of evenness profiles**

Chao and Ricotta (2019) developed five classes of evenness measures
parameterized by an order q \> 0. (For q = 0, species abundances are
disregarded, so it is not meaningful to evaluate evenness among
abundances specifically for q = 0. As q tends to 0, all evenness values
tend to 1 as a limiting value.) All classes of evenness measures are
functions of diversity and species richness, and all are standardized to
the range of \[0, 1\] to adjust for the effect of differing species
richness. Evenness profile depicts evenness estimate with respect to
order q ≥ 0. Because true species richness typically cannot be
accurately estimated, evenness profile typically can only be accurately
measured when both diversity and richness are computed at a fixed level
of sample coverage up to a maximum coverage Cmax defined in Step 3.
`iNEXT.4steps` shows, by default, the relevant statistics and plot for
only one class of evenness measure (based on the normalized slope of a
diversity profile), but all the five classes are featured.

<span style="color:blue;">NOTE 1</span>: Sufficient data are required to
perform the 4-step analysis. If there are only a few species in users’
data, it is likely that data are too sparse to use `iNEXT.4steps.`

<span style="color:blue;">NOTE 2</span>: The analyses in STEPs 2a, 2b
and 3 are mainly based on package `iNEXT` available from CRAN. Thus,
`iNEXT.4steps` expands `iNEXT` to include the estimation of sample
completeness and evenness.

<span style="color:blue;">NOTE 3</span>: As with `iNEXT`, `iNEXT.4steps`
only deals with taxonomic/species diversity. Researchers who are
interested in phylogenetic diversity and functional diversity should use
package `iNEXT.3D` available from CRAN and see the relevant paper (Chao
et al. 2021) for methodology.

<span style="color:blue;">NOTE 4</span>: `iNEXT.4steps` aims to compare
within-assemblage diversity. If the goal is to assess the extent of
differentiation among assemblages or to infer species compositional
shift and abundance changes, users should use `iNEXT.beta3D` available
from CRAN and see the relevant paper (Chao et al. 2023) for methodology.

## How to cite

If you publish your work based on results from iNEXT.4steps package, you
should make references to the following methodology paper and the
package:

-   Chao, A., Kubota, Y., Zelený, D., Chiu, C.-H., Li, C.-F., Kusumoto,
    B., Yasuhara, M., Thorn, S., Wei, C.-L., Costello, M. J. and
    Colwell, R. K. (2020). Quantifying sample completeness and comparing
    diversities among assemblages. Ecological Research, 35, 292-314.

-   Chao, A. and Hu, K.-H. (2024). The iNEXT.4steps package: Four-Step
    Biodiversity Analysis based on iNEXT. R package available from CRAN.

## SOFTWARE NEEDED TO RUN iNEXT.4STEPS IN R

-   Required: [R](https://cran.r-project.org/)
-   Suggested: [RStudio
    IDE](https://posit.co/products/open-source/rstudio/#Desktop)

## HOW TO RUN INEXT.4STEPS:

The `iNEXT.4steps` package can be downloaded from CRAN or Anne Chao’s
[iNEXT.4steps_github](https://github.com/AnneChao/iNEXT.4steps). For a
first-time installation, some additional packages must be installed and
loaded; see package manual.

``` r
## install iNEXT.4steps package from CRAN
install.packages("iNEXT.4steps")

## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('AnneChao/iNEXT.4steps')

## import packages
library(iNEXT.4steps)
```

An online version of iNEXT.4steps
(<https://chao.shinyapps.io/iNEXT_4steps/>) is also available for users
without an R background.

## <span style="color:red;">DATA INPUT FORMAT</span>

### Species abundance/incidence data format

For `iNEXT.4steps` package, pertinent information on species identity
(or any unique identification code) and assemblage affiliation is
required to be included in the input data for running iNEXT.4steps,
although species identity information is not taken into account in
inferring and comparing taxonomic/species diversity. Two types of
species taxonomic data are supported:

1.  Individual-based abundance data (`datatype = "abundance"`): When
    there are multiple assemblages, in addition to the assemblage/site
    names (as column names) and the species names (as row names),
    species abundance data (reference sample) can be input as a species
    (in rows) by assemblage (in columns) matrix/data.frame or a list of
    species abundance vectors. In the special case that there is only
    one assemblage, all data should be read in one column.

2.  Sampling-unit-based incidence data, i.e., Incidence or occurrence
    data (`datatype = "incidence_raw"`): for each assemblage, input data
    for a reference sample consist of a species-by-sampling-unit matrix,
    in addition to the sampling-unit names (as column names) and the
    species names (as row names). When there are N assemblages, input
    data consist of N lists of matrices, and each matrix is a
    species-by-sampling-unit matrix. Each element in the incidence raw
    matrix is 1 for a detection, and 0 for a non-detection. Input a
    matrix which combines data for all assemblages is allowed, but the
    argument `nT` must be specified to indicate the number of sampling
    units in each assemblage.

For example, the dataset `Data_spider` included in the `iNEXT.4steps`
package consists of species sample abundances of two assemblages/sites:
“Open” and “Closed”. Run the following code to view the data (only the
first 15 rows are shown below).

``` r
data("Data_spider")
Data_spider
```

                             Open Closed
    Pardosa_lugubris          350     10
    Alopecosa_taeniata        325     55
    Coelotes_terrestris       237    502
    Pardosa_riparia           102      1
    Haplodrassus_signifer      91      3
    Gnaphosa_badia             72      4
    Callobius_claustrarius     68    171
    Harpactea_lepida           61    140
    Tenuiphantes_tenebricola   53    180
    Cybaeus_angustiarum        50     24
    Inermocoelotes_inermis     50     76
    Pardosa_ferruginea         37     10
    Histopona_torpida          30     28
    Xerolycosa_nemoralis       29      0
    Amaurobius_fenestralis     17     35

We use incidence raw data (`Data_woody_plant`) collected from two sites,
`"Monsoon"` and `"Upper_cloud"`, as an example. Run the following code
to view the data (only the first 6 rows and first 3 columns for each
site are shown below).

``` r
data("Data_woody_plant")
Data_woody_plant
```

    $Monsoon
                                               Monsoon_unit_1 Monsoon_unit_2 Monsoon_unit_3
    Acacia_confusa_Merr.                                    0              0              0
    Acer_kawakamii_Koidzumi                                 0              0              0
    Acer_morrisonense_Hayata                                0              0              0
    Acer_serrulatum_Hayata                                  0              0              0
    Actinidia_chinensis_Planch._var._setosa_Li              0              0              0
    Adinandra_formosana_Hayata                              0              1              0

    $Upper_cloud
                                           Cloud_unit_1 Cloud_unit_2 Cloud_unit_3
    Abelia_chinensis_R._Br._var._ionandra             0            0            0
    Abies_kawakamii_(Hayata)_Ito                      0            0            0
    Acer_albopurpurascens_Hayata                      0            0            0
    Acer_kawakamii_Koidzumi                           0            0            0
    Acer_morrisonense_Hayata                          0            1            0
    Acer_palmatum_Thunb._var._pubescens_Li            0            0            0

## <span style="color:red;">MAIN FUNCTION iNEXT4steps()</span>

We first describe the main function `iNEXT4steps()` with default
arguments:

``` r
iNEXT4steps(data, q = seq(0, 2, 0.2), datatype = "abundance", 
            nboot = 30, conf = 0.95, nT = NULL, details = FALSE)
```

The arguments of this function are briefly described below, and will
explain details by illustrative examples in later text.

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
    (species by assemblages), or a list of species abundance vectors.
2.  For `datatype = 'incidence_raw'`, data can be input as a list of
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
    a numerical vector specifying the orders of q that will be used to
    compute sample completeness and evenness as well as plot the
    relevant profiles. Default is <code>seq(0, 2, by = 0.2)</code>.
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    datatype
    </td>
    <td style="text-align: left;">
    data type of input data: individual-based abundance data
    (`datatype = 'abundance'`) or species by sampling-units incidence
    matrix (`datatype = 'incidence_raw'`) with all entries being 0
    (non-detection) or 1 (detection)
    </td>
    </tr>
    <tr>
    <td style="text-align: left;">
    nboot
    </td>
    <td style="text-align: left;">
    a positive integer specifying the number of bootstrap replications
    when assessing sampling uncertainty and constructing confidence
    intervals. Enter 0 to skip the bootstrap procedures. Default is 30.
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
    (required only when <code>datatype = ‘incidence_raw’</code> and
    input data in a single matrix/data.frame) a vector of positive
    integers specifying the number of sampling units in each assemblage.
    If assemblage names are not specified (i.e., <code>names(nT) =
    NULL</code>), then assemblages are automatically named as
    ‘Assemblage1’, ‘Assemblage2’,…, etc.
    </td>
    </tr>
    <tr>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    details
    </td>
    <td style="border-bottom: 2px solid grey; text-align: left;">
    a logical variable to indicate whether the detailed numerical values
    for each step are displayed. Default is `FALSE`.
    </td>
    </tr>
    </tbody>
    </table>

The output of `iNEXT4steps` will have three parts (if `details = TRUE`):
`$summary`, `$figure`, and `$details`. It may take some time to compute
when data size is large or `nboot` is large.

## <span style="color:blue;">iNEXT.4steps VIA EXAMPLES</span>

First, we use the data `Data_spider` to illustrate the complete 4-step
analysis.

### EXAMPLE 1: Complete 4 steps for abundance data

In the spider data, species abundances of epigeal spiders were recorded
in two forest stands (“closed” and “open”). In the open forest, there
were 1760 individuals representing 74 species, whereas in the closed
forest, there were 1411 individuals representing 44 species. In the
pooled habitat, a total of 3171 individuals belonging to 85 species are
recorded.

Run the following code to obtain the numerical output and six figures
including five individual figures (for STEPs 1, 2a, 2b, 3 and 4,
respectively) and a complete set of five plots. (Here only show the
complete set of five plot; all five individual plots are omitted.)

``` r
data(Data_spider)
Four_Steps_out1 <- iNEXT4steps(data = Data_spider, datatype = "abundance")
Four_Steps_out1$summary
$`STEP 1. Sample completeness profiles`
  Assemblage q = 0 q = 1 q = 2
1     Closed  0.61  0.99     1
2       Open  0.77  0.99     1

$`STEP 2b. Observed diversity values and asymptotic estimates`
  Assemblage               qTD TD_obs TD_asy  s.e. qTD.LCL qTD.UCL
1     Closed  Species richness  44.00  72.11 23.16   26.72  117.49
2     Closed Shannon diversity  10.04  10.30  0.31    9.69   10.91
3     Closed Simpson diversity   5.71   5.73  0.20    5.33    6.12
4       Open  Species richness  74.00  96.31 12.46   71.89  120.73
5       Open Shannon diversity  16.34  16.84  0.59   15.68   18.00
6       Open Simpson diversity   9.41   9.46  0.32    8.83   10.09

$`STEP 3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.994 q = 0 q = 1 q = 2
1       Closed 55.62 10.18  5.72
2         Open 86.51 16.59  9.43

$`STEP 4. Evenness among species abundances of orders q = 1 and 2 at Cmax based on the normalized slope of a diversity profile`
  Cmax = 0.994 Pielou J' q = 1 q = 2
1       Closed      0.58  0.17  0.09
2         Open      0.63  0.18  0.10
Four_Steps_out1$figure[[6]]
```

<img src="README/README-unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" />

`$summary`: numerical tables for STEPs 1, 2b, 3 and 4.

-   `Assemblage` = the assemblage names.
-   `qTD` = ‘Species richness’ represents the taxonomic diversity of
    order q=0; ‘Shannon diversity’ represents the taxonomic diversity of
    order q=1, ‘Simpson diversity’ represents the taxonomic diversity of
    order q=2.
-   `TD_obs` = the observed taxonomic diversity value of order q.
-   `TD_asy` = the estimated asymptotic diversity value of order q.
-   `s.e.` = the bootstrap standard error of the estimated asymptotic
    diversity of order q.
-   `qTD.LCL`, `qTD.UCL` = the bootstrap lower and upper confidence
    limits for the estimated asymptotic diversity of order q at the
    specified level in the setting (with a default value of 0.95).
-   `Pielou J'` = a widely used evenness measure based on Shannon
    entropy.

`$figure`: six figures including five individual figures (for STEPS 1,
2a, 2b, 3 and 4 respectively) and a complete set of five plots.

`$details`: (only when `details = TRUE`). The numerical output for
plotting all figures.

### EXAMPLE 2: Complete 4 steps for incidence data

In the “Woody_plant” data, species incidence-raw data were recorded in
two forest vegetation types (“Monsoon” and “Upper_cloud” forest). In the
monsoon forest, 329 species and 6814 incidences were recorded in 191
plots. In the upper cloud forest, 239 species and 3371 incidences were
recorded in 153 plots (each 20×20-m plot is regarded as a sampling
unit). Because spatial clustering prevails in woody plants, individual
plants cannot be regarded as independent sampling units, violating the
basic sampling assumptions for the model based on abundance data. Thus,
it is statistically preferable to use incidence data to avoid this
violation.

Run the following code to obtain the numerical output and six figures
including five individual figures (for STEPs 1, 2a, 2b, 3 and 4,
respectively) and a complete set of five plots. (Here only show the
complete set of five plot; all five individual plots are omitted.)

``` r
data(Data_woody_plant)
Four_Steps_out2 <- iNEXT4steps(data = Data_woody_plant, datatype = "incidence_raw")
Four_Steps_out2$summary
$`STEP 1. Sample completeness profiles`
   Assemblage q = 0 q = 1 q = 2
1     Monsoon  0.78  0.99     1
2 Upper_cloud  0.78  0.98     1

$`STEP 2b. Observed diversity values and asymptotic estimates`
   Assemblage               qTD TD_obs TD_asy  s.e. qTD.LCL qTD.UCL
1     Monsoon  Species richness 329.00 421.67 20.20  382.07  461.27
2     Monsoon Shannon diversity 145.65 150.15  1.81  146.61  153.69
3     Monsoon Simpson diversity 102.33 103.35  1.54  100.33  106.36
4 Upper_cloud  Species richness 239.00 307.78 17.30  273.86  341.69
5 Upper_cloud Shannon diversity 105.53 110.50  1.88  106.81  114.19
6 Upper_cloud Simpson diversity  71.17  72.23  1.10   70.06   74.39

$`STEP 3. Non-asymptotic coverage-based rarefaction and extrapolation analysis`
  Cmax = 0.993  q = 0  q = 1  q = 2
1      Monsoon 359.80 147.29 102.67
2  Upper_cloud 278.96 108.52  71.69

$`STEP 4. Evenness among species abundances of orders q = 1 and 2 at Cmax based on the normalized slope of a diversity profile`
  Cmax = 0.993 Pielou J' q = 1 q = 2
1      Monsoon      0.85  0.41  0.28
2  Upper_cloud      0.83  0.39  0.25
Four_Steps_out2$figure[[6]]
```

<img src="README/README-unnamed-chunk-9-1.png" width="100%" style="display: block; margin: auto;" />

## <span style="color:red;">Completeness and ggCompleteness: MAIN FUNCTIONS FOR STEP 1</span>

Function `Completeness()` computes sample completeness estimates of
orders q = 0 to q = 2 in increments of 0.2 (by default), and function
`ggCompleteness` is used to plot the corresponding sample completeness
profiles. These two functions are specifically for users who only
require sample completeness estimates and profiles. The two functions
with arguments are described below:

``` r
Completeness(data, q = seq(0, 2, 0.2), datatype = "abundance", nboot = 30, 
             conf = 0.95, nT = NULL)
```

``` r
ggCompleteness(output)
```

All the arguments in these two functions are the same as those in the
main fnction `iNEXT4steps` for details.

### Sample completeness estimates and profiles for abundance data

Run the following code to obtain sample completeness estimates based on
the abundance data `Data_spider`:

``` r
data(Data_spider)
SC_out1 <- Completeness(data = Data_spider, datatype = "abundance")
SC_out1
   Order.q Estimate.SC         s.e.    SC.LCL    SC.UCL Assemblage
1      0.0   0.7683622 9.973432e-02 0.5728865 0.9638379       Open
2      0.2   0.8181065 7.486201e-02 0.6713796 0.9648333       Open
3      0.4   0.8768761 4.564035e-02 0.7874227 0.9663296       Open
4      0.6   0.9302044 2.156722e-02 0.8879335 0.9724754       Open
5      0.8   0.9664828 8.007891e-03 0.9507876 0.9821779       Open
6      1.0   0.9858045 2.445838e-03 0.9810107 0.9905982       Open
7      1.2   0.9944729 6.459167e-04 0.9932069 0.9957388       Open
8      1.4   0.9979624 1.730570e-04 0.9976232 0.9983016       Open
9      1.6   0.9992757 6.790607e-05 0.9991426 0.9994088       Open
10     1.8   0.9997489 3.316620e-05 0.9996839 0.9998140       Open
11     2.0   0.9999146 1.548913e-05 0.9998842 0.9999449       Open
12     0.0   0.6102206 1.617713e-01 0.2931547 0.9272866     Closed
13     0.2   0.7182647 1.283739e-01 0.4666565 0.9698729     Closed
14     0.4   0.8338575 7.549952e-02 0.6858812 0.9818339     Closed
15     0.6   0.9215916 3.005453e-02 0.8626858 0.9804974     Closed
16     0.8   0.9692426 8.931708e-03 0.9517368 0.9867484     Closed
17     1.0   0.9893733 2.393479e-03 0.9846821 0.9940644     Closed
18     1.2   0.9966174 6.939717e-04 0.9952573 0.9979776     Closed
19     1.4   0.9989804 2.335822e-04 0.9985226 0.9994382     Closed
20     1.6   0.9997042 8.447411e-05 0.9995387 0.9998698     Closed
21     1.8   0.9999166 3.060301e-05 0.9998566 0.9999766     Closed
22     2.0   0.9999770 1.087748e-05 0.9999557 0.9999983     Closed
```

-   `Order.q`: the order of sample completeness.
-   `Estimate.SC`: the estimated sample completeness of order q.
-   `s.e.`: standard error of sample completeness estimate.
-   `SC.LCL`, `SC.UCL`: the bootstrap lower and upper confidence limits
    for the sample completeness of order q at the specified level (with
    a default value of 0.95).
-   `Assemblage`: the assemblage name.

Run the following code to plot sample completeness profiles for q
between 0 to 2, along with confidence intervals.

``` r
ggCompleteness(SC_out1)
```

<img src="README/README-unnamed-chunk-13-1.png" width="70%" style="display: block; margin: auto;" />

### Sample completeness estimates and profiles for incidence data

Similar procedures can be applied to incidence data `Data_woody_plant`
to infer sample completeness.

``` r
data(Data_woody_plant)
SC_out2 <- Completeness(data = Data_woody_plant, datatype = "incidence_raw")
SC_out2
   Order.q Estimate.SC         s.e.    SC.LCL    SC.UCL  Assemblage
1      0.0   0.7802244 5.443597e-02 0.6735319 0.8869170     Monsoon
2      0.2   0.8490462 3.596453e-02 0.7785570 0.9195354     Monsoon
3      0.4   0.9086970 1.922901e-02 0.8710089 0.9463852     Monsoon
4      0.6   0.9508262 8.389296e-03 0.9343835 0.9672690     Monsoon
5      0.8   0.9758389 3.089211e-03 0.9697842 0.9818936     Monsoon
6      1.0   0.9888942 9.979319e-04 0.9869383 0.9908501     Monsoon
7      1.2   0.9951304 3.172420e-04 0.9945086 0.9957522     Monsoon
8      1.4   0.9979362 1.370679e-04 0.9976676 0.9982049     Monsoon
9      1.6   0.9991474 7.772191e-05 0.9989951 0.9992997     Monsoon
10     1.8   0.9996547 4.360856e-05 0.9995693 0.9997402     Monsoon
11     2.0   0.9998625 2.285966e-05 0.9998176 0.9999073     Monsoon
12     0.0   0.7765330 6.598825e-02 0.6471984 0.9058676 Upper_cloud
13     0.2   0.8358241 4.572396e-02 0.7462068 0.9254414 Upper_cloud
14     0.4   0.8915125 2.644153e-02 0.8396880 0.9433369 Upper_cloud
15     0.6   0.9354738 1.265255e-02 0.9106752 0.9602723 Upper_cloud
16     0.8   0.9649731 5.068840e-03 0.9550384 0.9749079 Upper_cloud
17     1.0   0.9823020 1.729246e-03 0.9789128 0.9856913 Upper_cloud
18     1.2   0.9915212 5.621355e-04 0.9904194 0.9926230 Upper_cloud
19     1.4   0.9960932 2.764244e-04 0.9955514 0.9966350 Upper_cloud
20     1.6   0.9982514 1.843399e-04 0.9978901 0.9986127 Upper_cloud
21     1.8   0.9992348 1.156412e-04 0.9990081 0.9994614 Upper_cloud
22     2.0   0.9996711 6.632610e-05 0.9995411 0.9998011 Upper_cloud
```

Run the following code to plot sample completeness profiles for q
between 0 to 2, along with confidence intervals.

``` r
ggCompleteness(SC_out2)
```

<img src="README/README-unnamed-chunk-15-1.png" width="70%" style="display: block; margin: auto;" />

## <span style="color:red;">Evenness and ggevenness: MAIN FUNCTIONS FOR STEP 4</span>

`Evenness()` computes standardized (or observed) evenness of order q = 0
to q = 2 in increments of 0.2 (by default) based on five classes of
evenness measures, and function `ggevenness` is used to plot the
corresponding evenness profiles. These two functions are specifically
for users who only require evenness estimates and profiles. The two
functions with arguments are described below:

``` r
Evenness(data, q = seq(0, 2, 0.2), datatype = "abundance", method = "Estimated",
         nboot = 30, conf = 0.95, nT = NULL, E.class = 1:5, SC = NULL)
```

``` r
ggEvenness(output)
```

There are only three arguments that are not used in the main function
`iNEXT4steps`; these three arguments are described below (see
`iNEXT4steps` for other arguments)

<table class="gmisc_table" style="border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;">
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
method
</td>
<td style="text-align: left;">
a binary selection of method with <code>‘Estimated’</code> (evenness is
computed under a standardized coverage value) or <code>‘Observed’</code>
(evenness is computed for the observed data).
</td>
</tr>
<tr>
<td style="text-align: left;">
E.class
</td>
<td style="text-align: left;">
an integer vector between 1 to 5 specifying which class(es) of evenness
measures are selected; default is 1:5 (select all five classes).
</td>
</tr>
<tr>
<td style="border-bottom: 2px solid grey; text-align: left;">
SC
</td>
<td style="border-bottom: 2px solid grey; text-align: left;">
(required only when <code>method = ‘Estimated’</code>) a standardized
coverage value for calculating estimated evenness. If
<code>SC=NULL</code>, then this function computes the diversity
estimates for the minimum sample coverage among all samples extrapolated
to double reference sizes (Cmax).
</td>
</tr>
</tbody>
</table>

Two simple examples for demonstrating functions `Evenness` and
`ggEvenness` are given below.

### Evenness estimates and profiles for abundance data with default standardized coverage value (Cmax)

The dataset `Data_spider` is used to estimate evenness at the default
standardized sample coverage (`SC = NULL`). Evenness estimates of order
q = 0 to q = 2 in increments of 0.2 (by default) will be computed based
on five classes of evenness measures. (For q = 0, species abundances are
disregarded, so it is not meaningful to evaluate evenness among
abundances specifically for q = 0. As q tends to 0, all evenness values
tend to 1 as a limiting value.) Function `ggevenness` is used to plot
the corresponding evenness profiles. Here only the evenness estimates
are shown for the first class of evenness measures. The corresponding
numerical tables for the other four classes of evenness measures are
omitted. Users also can use argument `E.class` to specify which class
(e.g., `E.class= 3`) or classes are required.

``` r
data(Data_spider)
Even_out1_est <- Evenness(data = Data_spider, datatype = "abundance", 
                          method = "Estimated", SC = NULL, E.class = 1:5)
Even_out1_est
```

    $E1
       Order.q  Evenness        s.e.  Even.LCL  Even.UCL Assemblage    Method        SC
    1      0.0 1.0000000 0.000000000 1.0000000 1.0000000       Open Estimated 0.9937676
    2      0.2 0.7276194 0.023738926 0.6810919 0.7741468       Open Estimated 0.9937676
    3      0.4 0.6008010 0.030661740 0.5407051 0.6608969       Open Estimated 0.9937676
    4      0.6 0.5628529 0.030325802 0.5034154 0.6222904       Open Estimated 0.9937676
    5      0.8 0.5801167 0.026422376 0.5283298 0.6319036       Open Estimated 0.9937676
    6      1.0 0.6297250 0.020846901 0.5888658 0.6705841       Open Estimated 0.9937676
    7      1.2 0.6941362 0.015075231 0.6645892 0.7236831       Open Estimated 0.9937676
    8      1.4 0.7600773 0.010243046 0.7400013 0.7801533       Open Estimated 0.9937676
    9      1.6 0.8190763 0.006881747 0.8055883 0.8325643       Open Estimated 0.9937676
    10     1.8 0.8673381 0.004872872 0.8577874 0.8768887       Open Estimated 0.9937676
    11     2.0 0.9044427 0.003733457 0.8971253 0.9117602       Open Estimated 0.9937676
    12     0.0 1.0000000 0.000000000 1.0000000 1.0000000     Closed Estimated 0.9937676
    13     0.2 0.6996114 0.033769758 0.6334239 0.7657989     Closed Estimated 0.9937676
    14     0.4 0.5665175 0.045015232 0.4782893 0.6547457     Closed Estimated 0.9937676
    15     0.6 0.5253117 0.046017907 0.4351183 0.6155052     Closed Estimated 0.9937676
    16     0.8 0.5365988 0.041426598 0.4554041 0.6177934     Closed Estimated 0.9937676
    17     1.0 0.5775431 0.033821777 0.5112536 0.6438325     Closed Estimated 0.9937676
    18     1.2 0.6330777 0.025399022 0.5832966 0.6828589     Closed Estimated 0.9937676
    19     1.4 0.6926130 0.017983239 0.6573665 0.7278595     Closed Estimated 0.9937676
    20     1.6 0.7491390 0.012664600 0.7243169 0.7739612     Closed Estimated 0.9937676
    21     1.8 0.7988183 0.009586218 0.7800296 0.8176069     Closed Estimated 0.9937676
    22     2.0 0.8402207 0.008084184 0.8243760 0.8560654     Closed Estimated 0.9937676

-   `Order.q`: the order of evenness.
-   `Evenness`: the computed evenness value of order q.
-   `s.e.`: standard error of evenness value.
-   `Even.LCL`, `Even.UCL`: the bootstrap lower and upper confidence
    limits for the evenness of order q at the specified level (with a
    default value of `0.95`).
-   `Assemblage`: the assemblage name.
-   `Method`: `"Estimated"` or `"Observed"`.
-   `SC`: the standardized coverage value under which evenness values
    are computed (only for `method = "Estimated"`)

The following commands plot the evenness profiles for all five classes
of even measures, along with their confidence intervals.

``` r
ggEvenness(Even_out1_est)
```

<img src="README/README-unnamed-chunk-21-1.png" width="100%" style="display: block; margin: auto;" />

### Evenness estimates and profiles for incidence data with user’s specified coverage value of 0.98.

In the function `Evenness`, users can specify a particular sample
coverage value under which all five classes of evenness measures will be
computed. For example, instead of using the default standardized
coverage value, if users want to compute evenness estimates for 0.98
based on the incidence dataset `Data_woody_plant`, then the argument
`SC=0.98` is used instead of `SC=NULL`, as shown below. Here only the
evenness estimates are displayed for the first class of evenness
measures. The corresponding numerical tables for the other four classes
of evenness measures are omitted.

``` r
data(Data_woody_plant)
Even_out2_est <- Evenness(data = Data_woody_plant, datatype = "incidence_raw", 
                          method = "Estimated", SC = 0.98, E.class = 1:5)
Even_out2_est
```

    $E1
       Order.q  Evenness         s.e.  Even.LCL  Even.UCL  Assemblage    Method   SC
    1      0.0 1.0000000 0.0000000000 1.0000000 1.0000000     Monsoon Estimated 0.98
    2      0.2 0.8756341 0.0040909678 0.8676160 0.8836523     Monsoon Estimated 0.98
    3      0.4 0.8206746 0.0056212690 0.8096572 0.8316921     Monsoon Estimated 0.98
    4      0.6 0.8142938 0.0055707897 0.8033753 0.8252124     Monsoon Estimated 0.98
    5      0.8 0.8394339 0.0045936444 0.8304305 0.8484373     Monsoon Estimated 0.98
    6      1.0 0.8797243 0.0032418153 0.8733705 0.8860782     Monsoon Estimated 0.98
    7      1.2 0.9206306 0.0019800198 0.9167499 0.9245114     Monsoon Estimated 0.98
    8      1.4 0.9530254 0.0010677797 0.9509326 0.9551183     Monsoon Estimated 0.98
    9      1.6 0.9745076 0.0005256135 0.9734775 0.9755378     Monsoon Estimated 0.98
    10     1.8 0.9870261 0.0002460964 0.9865437 0.9875084     Monsoon Estimated 0.98
    11     2.0 0.9936878 0.0001141111 0.9934641 0.9939114     Monsoon Estimated 0.98
    12     0.0 1.0000000 0.0000000000 1.0000000 1.0000000 Upper_cloud Estimated 0.98
    13     0.2 0.8606066 0.0065549808 0.8477590 0.8734541 Upper_cloud Estimated 0.98
    14     0.4 0.7970245 0.0090550828 0.7792769 0.8147721 Upper_cloud Estimated 0.98
    15     0.6 0.7868162 0.0091104141 0.7689601 0.8046723 Upper_cloud Estimated 0.98
    16     0.8 0.8119024 0.0077189656 0.7967735 0.8270313 Upper_cloud Estimated 0.98
    17     1.0 0.8550606 0.0056813182 0.8439255 0.8661958 Upper_cloud Estimated 0.98
    18     1.2 0.9007350 0.0036838034 0.8935148 0.9079551 Upper_cloud Estimated 0.98
    19     1.4 0.9385213 0.0021499901 0.9343074 0.9427352 Upper_cloud Estimated 0.98
    20     1.6 0.9648542 0.0011662102 0.9625685 0.9671399 Upper_cloud Estimated 0.98
    21     1.8 0.9810635 0.0006089086 0.9798701 0.9822570 Upper_cloud Estimated 0.98
    22     2.0 0.9902087 0.0003148737 0.9895916 0.9908259 Upper_cloud Estimated 0.98

The following commands plot the evenness profiles for five classes,
along with their confidence intervals.

``` r
ggEvenness(Even_out2_est)
```

<img src="README/README-unnamed-chunk-24-1.png" width="100%" style="display: block; margin: auto;" />

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

-   Chao, A., Kubota, Y., Zelený, D., Chiu, C.-H., Li, C.-F., Kusumoto,
    B., Yasuhara, M., Thorn, S., Wei, C.-L., Costello, M. J. and
    Colwell, R. K. (2020). Quantifying sample completeness and comparing
    diversities among assemblages. Ecological Research, 35, 292-314.

-   Chao, A. and Ricotta, C. (2019). Quantifying evenness and linking it
    to diversity, beta diversity, and similarity. Ecology, 100(12),
    e02852.

-   Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R.
    L., Wu, J., Magnago, L. F. S., Dornelas, M., Zelený, D., Colwell, R.
    K., and Magurran, A. E. (2023). Rarefaction and extrapolation with
    beta diversity under a framework of Hill numbers: the iNEXT.beta3D
    standardization. Ecological Monographs e1588.
