# Potential Duplicate Pairs Analysis

Analysis of potential duplicate sequence pairs in the CAM_1000 dataset with varying day interval thresholds.

## Summary Table

**Total pairs in database: 36,141,038**

| Day Interval | Number of Potential Duplicate Pairs | Percent of Total Pairs |
|------------------|---------------------------------|----------------------|
| 0            | 18,830                              | 0.052%                 |
| 3            | 42,617                              | 0.118%                 |
| 4            | 49,286                              | 0.136%                 |
| 5            | 55,236                              | 0.153%                 |
| 6            | 60,605                              | 0.168%                 |
| 7            | 66,857                              | 0.185%                 |
| 14           | 80,589                              | 0.223%                 |
| 21           | 86,049                              | 0.238%                 |
| 28           | 89,324                              | 0.247%                 |

## Methodology

Potential duplicate pairs are defined as sequence pairs meeting ALL of the following criteria:

-   Identical sequences (0 mutations, `n_mutations == 0`)
-   Same geographic division (state/province)
-   Same age (age_adj)
-   Collection dates within the specified day interval
-   Matching sex (or both NA)
-   Matching sub-state/province location (or both NA)

## Observations

-   The number of potential duplicates increases substantially with longer day intervals
-   From day 0 to day 3, there is a 126% increase (23,787 additional pairs)
-   The rate of increase slows considerably after day 14
-   Between days 21 and 28, only 3,275 additional pairs are identified (3.8% increase)
-   This suggests that most potential duplicates occur within a 2-3 week window

## Analysis Date

Generated: 2026-01-05