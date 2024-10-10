# **NCOV-USA-MIG Interval Report**

*Author: Amin Bemanian* <br />
*Date: October 7, 2024*

# Objectives

The goal of this document is to provide a summary of the analyses completed so far for this project and to help with writing an eventual manuscript. For this report, I have completed the basic age group and state-to-state relative risk calculations. Additionally, I looked at identifying state-level clusters and the relationship of travel data with state-level RRs.

# Findings

## Overview of Dataset

Data was pulled on July 16, 2024 from GISAID dataset. Criteria for inclusion were sequences with at least 90% coverage, from the United States, within the 50 states or the District of Columbia. Sequences with missing state or invalid state data or from other US territories (e.g. US Virgin Islands, Guam) were excluded from the analysis. This resulted in 4,902,211 sequences and 994.579,915 pairs of identical sequences being included in the final analysis. 

## Age Analysis

36% of the pairs of identical sequences had appropriate age data for both individual sequences. RR analysis was done with ages binned in 3 year groups (0-2 year olds, 3-5 year olds, etc.) up to 90 years of age. The age distributions of individuals and pairs with at least one sequence in an age group are shown below. In both of these we find that the largest age group of individuals and pairs is during the young adult years (20s to 30s). This age group is slightly more represented in pairs compared to the expedcted based on the individual age distribution.

<img src= '../figs/USA/age_seq_hist.jpg' alt = 'Age distribution of individuals' width = '900' />

<img src= '../figs/USA/age_pairs_hist.jpg' alt = 'Age distribution of pairs' width = '900' />

We then calculate the relative risks of identical sequences across these age groups. This is displayed below as a heatmap and then as a series of line plots stratified by each age group. Looking at these figures we find that the there is generally a trend of higher RR within the same age groups and neighboring age groups. There seems to be 3 large categories of age groups that stand out in this analysis:

- Pediatric (0-20 year olds)
- Early Adulthood (20-49 year olds)
- Later Adulthood (50+ year olds)

<img src= '../figs/USA/age_heatmap.jpg' alt = 'Age RR heatmap' width = '1500' />

<img src= '../figs/USA/age_lineplot_all.jpg' alt = 'Age RR heatmap' width = '1500' />

The pediatric age groups all have RRs > 1 with each other and then tend to have their peaks at with other groups in the same school structure (example 6-8 and 9-11 overlap in elementary school and are highly related). These groups also all have RRs < 1 with the later adulthood age groups, with the notable exception of the youngest group (0-2 yo). The youngest group has RRs > 1 with the later adulthood group with peaks in the 70s. This may be explained by the role of grandparents in helping with the care of infants and young babies.

The early adulthood age groups notably do not have strong peaks or valleys with any of the other age groups. They do have RRs > 1 on average with the pediatric group, but the absolute value of these peaks is often barely at 1.01. Unfortunately I was not able to run a bootstrap analysis for this dataset to check for significance. (see [Performance](#performance) under [Next Steps](#next-steps)) The later adulthood groups have a very strong relationship with eachother, with the highest RRs of any age being in the 80 and older age groups.

## State Analysis

For this section I calculated the relative risks across all 50 states and the District of Columbia. Earlier versions of this analysis included the other US territories, but due to how few sequences were available in the GISAID dataset for those territories, RRs were unable to be calculated with all pairs of states. State level relative risks are shown below in a heatmap. The states have been manually ordered by geographic division/neighbor adjacency:

<img src= '../figs/USA/state_heatmap.jpg' alt = 'State RR heatmap' width = '1500' />

Universally across all states, there is clear evidence that identical sequences are more likely to occur within a state compared to not. Additionally, there is evidence for several clusters of states which are more likely to conserve identical sequences. There is a particularly prominent cluster around the center of the figure spanning from Pennsylvania to Rhode Island. The states in this cluster make up 4 contigious major metropolitan areas: New York City (NY-PA-NJ-CT), Philadelphia (PA-NJ-DE-MD), Baltimore (MD), and Washington (DC-MD-VA). There is evidence for other clusters around several other multi-state metro areas including Boston (MA-NH-RI), Chicago (IL-WI-IN), Kansas City (MO-KS), and Portland (OR-WA); as well as other general regional clusters (e.g. Northern Mountain States (MT-WY-ID) and the Gulf Coast (TX-LA-AL-MS)). Maps of RR for each state was calculated as well and I have included a sample of 5 states here: Washington State, New York, Illinois, California, and Texas.

<img src= '../figs/USA/maps/rr_map_Washington.jpg' alt = 'Washington RR map' width = '900' />
<img src= '../figs/USA/maps/rr_map_New York.jpg' alt = 'New York RR map' width = '900' />
<img src= '../figs/USA/maps/rr_map_Illinois.jpg' alt = 'Illinois RR map' width = '900' />
<img src= '../figs/USA/maps/rr_map_California.jpg' alt = 'California RR map' width = '900' />
<img src= '../figs/USA/maps/rr_map_Texas.jpg' alt = 'Texas RR map' width = '900' />

Given that RR was high for neighboring states, I investigated the relationship of distance and state adjacency with the relative risk between states. The graphs for each are shown below as well as a boxplot comparing the RR with state adjacency. Again within state RRs were generally found to be higher than out-of-state RRs, and adjacent state RRs tended to be higher than non-adjacent states. The distance plots suggest that higher distance (either Euclidean spatial or neighbor adjacency) was associated with a reversion to an RR of 1.

<img src = '../figs/USA/state_adj_plot.jpg' alt = 'State RR Box Plot by Adjacency' width = '900' />
<img src = '../figs/USA/state_euclid_dist_plot.jpg' alt = 'State RR vs Euclidean Distance' width = '900' />
<img src = '../figs/USA/state_nb_dist_plot.jpg' alt = 'State RR vs Neighbor Adjacency' width = '900' />

## Cluster and Travel Analyses

Since there was evidence of state-level clustering, the next analysis I completed was trying several unsupervised cluster analyses. I used a negative-exponential function of the RR between two states as a distance function (i.e. $D(i,j) = e^{-RR_{i,j}}$). I tested hierarchical cluster analysis (HCA) and K-means clustering of a principle coordinate analysis (PCoA) and a non-metric multidimensional scaling analysis (nMDS). The results of the HCA are shown below with the dendogram, map. and the clusters plotted on the PCoA and nMDS dimensions.

<img src = '../figs/USA/state_dendrogram.jpg' alt = 'State HCA Dendrogram' width = 1500 />
<img src = '../figs/USA/maps/state_hclust.jpg' alt = 'State HCA Map' width = 1500 />
<img src = '../figs/USA/state_PCoA_hclust.jpg' alt = 'State PCoA w/ HCA plot' width = 1500 />
<img src = '../figs/USA/state_NMDS_hclust.jpg' alt = 'State PCoA w/ HCA plot' width = 1500 />

Looking at the HCA clusters we do find several clusters that are geographically concentrated: Cluster 5 which spans from North Carolina to Rhode Island centered around the DC-NY metro corridor; Cluster 10 which includes the remaidner of the New England states; Cluster 7 which includes all the states around Lake Michigan, Cluster 8 which includes most of the Upper Great Plains and Mountain States, and Cluster 3 which includes the West Coast and adjacent states (except for Washington). Cluster 5 is particularly prominent on the dendrogram as it is the first branch point. Additionally when the states are plotted in the PCoA-space, Cluster 5 is completely separated from the other clusters on Component #1. There are several other clusters that do not seem to have any clear geographical relationship such as Cluster #2 (AK-UT-MN-MS) or Cluster #4 (some Ohio River Valley states but also OK-AR and Florida). The HCA clusters do tend to gnerally stick together in the PCoA and nMDS plots as well.

Next, I wanted to test how predictive travel patterns were for RR between states. I used publically available [travel trip estimates](https://www.bts.gov/browse-statistical-products-and-data/state-transportation-statistics/interstate-passenger-trips) from the Bureau of Transportation Statistics (part of the USDOT). These were calculated from the [2022 National Household Travel Survey](https://nhts.ornl.gov/). The trips are partially broken down by mode. For my analysis I calculated travel RR between states and plotted these against identical sequence RR. Three graphs are presented below: all trip modalities between all pairs of states (including intra-state), all trip modalities for only interstate pairs, and air trips only for interstate pairs. For the air trips, DC was excluded since it does not have an airport physically located in the District.

<img src = '../figs/USA/travel_seq_all.jpg' alt = 'Travel RR vs Identical Sequence RR for All Trips and States' width = 1500 />
<img src = '../figs/USA/travel_seq_out.jpg' alt = 'Travel RR vs Identical Sequence RR for All Trips and Different States' width = 1500 />
<img src = '../figs/USA/travel_seq_air.jpg' alt = 'Travel RR vs Identical Sequence RR for Air Trips and Different States' width = 1500 />

For the first figure of all trips with all pairs of states, I used color to distinguish between interstate and intrastate RRs. Unsurprisingly, intrastate points had much higher travel and identical sequence RR scores, which is why I repeated the analysis with them removed. In the other two graphs the colors distinguish between pairs that are in the same HCA cluster versus pairs that are in different HCA clusters. Pairs within the same HCA cluster had higher identical sequence RRs as expected. The pairs that are in within the same HCA cluster do have higher travel RRs compared to pairs not in the same cluster for both all modes of travel and air travel specifically. Preliminary boxplots of the travel RRs stratified by cluster are below.

<img src = '../figs/USA/travel_cluster_out.jpg' alt = 'Travel RR vs Clusters - All Modes' width = 1500 />
<img src = '../figs/USA/travel_cluster_air.jpg' alt = 'Travel RR vs Clusters - Air' width = 1500 />

# General Conclusions

This report shows that we are able to leverage the full set of almost 5 million sequences into a single analysis to look at some of the social and geographic forces that influence COVID-19 transmission. There remains some computational challenges, but even with those there are some interesting conclusions. We managed to redemonstrate with the age analysis that identical sequences tend to be most likely conserved within the same and close age groups. From a geographic standpoint, there is evidence of clustering in spread between states. Additioanlly, there is evidence that travel behaviors drive this at least partially, but there may still be more to investigate within the clusters that were identified. 

# Next Steps

## Performance 
One major next step should be performance optimization. Interestingly the mutation distance calculation between sequences is not the rate limiting step. On the HPC cluster, completion through the distance step takes only about 24-36 hours with 32 threads requested. The next step of calculating the clusters of identical sequences still has not been able to be completed for the full USA scenario. This is likely due to the inefficiencny of memory structures in R. 

Another issue is with the computation time of relative risk. Rhe age and state RRs each take about 8-12 hours to calculate. This should be a fairly quick operation of filtering and restructuring the table to count the number of pairs of sequence with zero distance. I wrote this code in R since computationally it was not too intensive. This was also done with the goal of implementing bootstrapping with parallelization fairly easily. However, with the length per RR calculation, it hasn't been feasible to actually run the bootstrapping for confidence intervals for the full US dataset. Additionally, the parallelization code no longer works correctlylma, with parallel threads inconsistently reporting back to the main thread and failing to re-aggregate the results of bootstrapping. I suspect this has something to do with the memory management on the cluster nodes, especially over longer periods. I anticipate I will have to rewrite the code to use tsv-utils for more of the heavy table filtering and restructring, and have the bootstrapped simulations saved to disk to prevent memory loss.

## Age Analyses
Once the performance issues are handled, adding back confidence intervals for the line plots will be a top priority. We could try adjusting the binning to be down to every 2 years or even every year although splitting up age groups will likely come with further performance costs. Additionally with age groups that are this small, it is no longer practical to stratify the analysis by each state. It still would be worthwhile to stratify the age analyis between within state and out of state pairs to see if these age relationships are mainly driven by local interactions.

## State Clustering, Geography, Travel
One of the next directions for this section is determining how best to visually represent this information. There are some intersting stories here but it would be helpful to hear what other people think in terms of what stands out and if there are other specific relationships I should further pursue. I also was wanting input on my use of HCA and other clustering techniques, specifically if the findings seem meaningful and if there are other clustering techniques I should look at it.

## Temporality and Directionality
Big picture wise, I think this is probably the ultimate direction of interest. We have talked in the past about several different approaches to move forward with this (looking at the arrival time of sequences or strains, CTMC models, etc.) Additionally in this category I do think it may be worthwhile to do some sub-analysis on specfic time intervals (e.g pre and post vaccine roll-out) or specific strains of interest and see how relationships change.