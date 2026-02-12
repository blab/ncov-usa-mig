## Introduction 

* Review of prior work at pandemic transmission trends over large areas  
* Phylo limitations/feasibility  
* Discussion of Cécile’s prior work with WA state  
* Outline goals to spread to larger dataset across all of NA and to look at variation over time

## Methods
* GISAID set
* idseq flows for curating sequence data
* Data Cleaning - Age and sex variables needed to be modified prior 
* RR calculation and other descriptive stats
* Time normalization
* Clustering methods
* Outside datasets
* Regression calculation

## Result Structure

#### Figure 1 - Exploratory analysis of sequencing metadat
* Maps of sequence counts and sequencing effort per capita
* Age/sex pyramid of sequence data
* Sequencing effort over time (by country and by state)

#### Figure 2 - Global geographic structure of RR
* Small set of example RR maps (NY,CA,ON,MX as origins?)
* Region Maps paired with Hclust Ordered Relative Risk Heatmap
* PCoA plots of RR

#### Figure 3 - Variation of geographic RRs over time
* Normalized RR scores by region membership over time
* Fold-change of nRR over time for domestic vs intl pairs (need to split int'l between border and non-border, possibly do new statistical analysis for sig testing)
* Map of strong connections (nRR > certain threshold) over time

#### Figure 4 - Relationship to Distance/Movement
* Scatterplots of RR versus distance/movement
* Plot of the correlations over time for these variable
* Possibly marginal effects plot for spline interactions

#### Table 1 - Regression of Predictors for RR
* Top part - Single variable predictors, spearman corr + spline R-sq
* Middle part - Focus on spline model interactions with distance (e.g. NHTS + Dist and Air Travel + Dist)
* Bottom - Final multivariable model performance

#### Figure 5 - Age relationship with RR
* Heatmap across all age groups and across all of time
* 3 small age heatmaps of within state pairs, intraregion pairs, and interregion pairs
* Focused age line nRR plots over time for particularly interesting pairings (likely childhood with other age groups and parents age group)
* Trying to look at states which had a big shift from virtual to in person schooling and seeing if the RR or nRR heatmaps changed significantly at the beginning of the school year versus the end

## Discussion