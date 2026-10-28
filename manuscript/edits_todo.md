# Edits To-Do (from `\ctkc{}` and `\ctk{}` comments)

Line numbers for open items refer to `main.tex` as of 2026-10-16 (latest sync); struck items keep the line numbers they had when struck.

## Global / manuscript-wide

- [x] ~~Thorough typo and tense-consistency (past/present) pass across the whole manuscript, especially Discussion, Methods, and the simulation section (L406, L505, L668)~~ — done: tense convention + typo/grammar pass, Fig./Table capitalization
- [x] ~~Use "North American subdivisions" rather than "states/provinces" everywhere, except where Mexico is genuinely excluded (L582)~~ — done: line-by-line pass in text + figure labels (age, geo_time)
- [x] ~~Rename the SafeGraph metric to "cellphone mobility" / "cellphone RR" everywhere (Fig. B labels, Results, Methods) because all the metrics measure mobility (L346, L525)~~ — done: renamed "smartphone mobility" throughout
- [x] ~~Rename "high-ranking across region pairs" to something easier to interpret (e.g. "persistent cross-region connections") and use the new name consistently in Methods, Results, and Discussion (L594)~~ — done: "persistent cross-region connections" + "persistent seasonal connections"; supp fig uses "percentile of cross-region connections"
- [ ] Make sure every supplementary figure and analysis is referenced in the main text (L720)

## Results

### Temporal and Seasonal Relationships for Geographic RR
- [x] ~~**L235** Clarify figure legends: in B, the color ramp should say nRR fold change "compared with 2023 Q1"; in C, "# of time points" should become something like "number of quarters with high connection", and "# of years" should say what is being plotted~~ — done: B legend "nRR Fold Change Compared with 2023 Q1"; C legends "Number of Quarters with High Connection" / "Number of Quarters"; caption explains seasonal stratification

### Geographical Predictors of Transmission
- [x] **L343** The section is long and hard to follow. Shorten it and make the paragraphs more direct
- [x] ~~**L346** (Cécile likes the new Fig. B.) Relabel SafeGraph as "cellphone mobility" (see Global)~~ — done: relabeled "Smartphone Mobility RR" in Fig. 2A/2B

### Age Relationships with Identical Sequence Transmission
- [x] ~~**L359** Figure caption: say what the shaded area in panel D represents~~ — done: removed GAM confidence bands from panel D instead
- [x] ~~**L357** *FYI:* Cécile rephrased to say transmission *rates* are highest between similar ages (absolute counts may be higher elsewhere). Check that the edit reads well~~ — approved Cécile's rewording
- [x] ~~**L381** *FYI:* praise for the contact-matrix correlation paragraph. No action needed~~
- [x] ~~**L371** Note that strong diagonal enrichment may come from repeated sampling of the same individual (Trevor's point too). Consider a supplementary heatmap of age RR restricted to pairs of different sex to support this~~ — done: Results pointer to Discussion + Supp Fig age_gender; filter scope/lineage wording fixed
- [x] ~~**L383** Deviance definition: drop "or any other deme definition" unless deviance was actually computed for other demes~~ — done: phrase and comment removed

## Discussion

- [x] ~~**L394** Open the Discussion with a short paragraph summarizing the main findings, then move on to interpretation and implications~~ — done: summary paragraph added
- [x] ~~**L406** Fix tense inconsistencies (see Global)~~ — done
- [x] ~~**L401** Clarify the "preprocessing" step: computing identical-sequence pairs / the distance matrix and storing them in DuckDB. Note that this is a one-time step and the database can then be reused for many analyses~~ — done: preprocessing step explained as one-time pair list + relational DB
- [x] ~~**L411** Add a reference (the claim that studies usually focus on a smaller area or subsample and aggregate)~~ — done: Matteson 2023, Deng 2020, Lemey 2021, Tsui 2023
- [x] ~~**L415** Add a figure reference (Southeast/Far West states with more cross-region than within-region transmission)~~ — done
- [x] ~~**L450** Direction of transmission: at least mention using temporal information, and possibly using pairs one mutation apart and identifying the more ancestral sequence~~ - added
- [x] ~~**L461** Add a reference to the supplementary figure (correlation with synthetic contact matrices)~~ - done
- [x] !!**L464** *FYI:* Cécile removed the comment on the slope being <1, because the slope mostly reflects transmission intensity and signal dilution. Confirm you agree~~ - agree
- [x] ~~**L479** The Florida/seasonal paragraph is too niche to end the paper on (agreeing with your `\abc` note). Replace it with a short concluding paragraph that restates the main contributions and ties back to the big-picture points from the Introduction~~ — done: Florida paragraph replaced with limitations + outlook close

## Methods

### Data Sources
- [x] ~~**L496** Add all the GISAID citations required (check the GISAID website)~~ — done: added Elbe 2017 and Shu 2017 alongside Khare 2021
- [x] ~~**L497** Say that aligned genomes came from the Nextstrain private ncov pipeline (see how Katie/Cécile worded it)~~ — done: Nextstrain nCoV ingest pipeline sentence + Hadfield 2018
- [x] ~~**L499** Say where coverage came from (Nextclade annotations?)~~ — done: coverage and Pango lineages assigned by Nextclade
- [x] ~~**L517** Explain why age results stop at 80 (binning? cutoff?)~~ - done: will include full age results in RR matrix in repo
- [ ] **L496** Describe how the NHTS ground-travel RR was computed, or cite the WA paper and say the same method was used
- [x] ~~**L525** Rename the SafeGraph metric (see Global)~~ — done: "smartphone mobility" throughout; SafeGraph only named as data source in Methods

### Identical Sequence Relative Risk Calculation
- [x] ~~**L559** Add citations: https://doi.org/10.1126/science.abe8372 and https://doi.org/10.1038/s41467-021-27163-1 (young adults disproportionately driving transmission)~~ — done: cited Monod et al. 2021 and Tran Kiem et al. 2021
- [ ] **L552** Reorder: describe geographic RR normalization first, then age RR (age is supplementary)

### Geographic Clustering and Regional Analysis
- [ ] **L558** Check the edit saying the dissimilarity matrix is also used for dimensionality reduction (Cécile suggested wording)
- [ ] **L565** Briefly explain why only pairs of subdivisions in *different* regions are ranked
- [x] ~~**L594** Rename "high-ranking across region pairs" (see Global)~~ — done
- [ ] **L569** Fill in the "xx" blanks in the seasonal-connection example sentence with a real example

### Distance/Mobility Regression Analysis
- [ ] **L576** Rephrase the exclusion criteria, e.g. "We restricted this comparison to pairs of states with reported travel or movement. As movement scores weren't available within the same state, we excluded all same-state pairs from this analysis."

### Acknowledgements
- [ ] **L583** Add the link to the GISAID acknowledgement table and the GISAID citation

## Supplemental Information

### Computational workflow
- [ ] **L606** Add a reference for `pairSNP`
- [ ] **L608** Mention the CLI option `-d 0`
- [ ] **L610** Say what the pre-processing steps were

### Simulation study
- [ ] **L631** Say how the negative binomial (contacts per day) was parameterized, and how the Beta-distributed infection probability was parameterized and used
- [ ] **L632** Add parameter values for the Gamma incubation period
- [ ] **L634** Say how the binomial draws for infection duration were parameterized
- [ ] **L639** Make sure the analysis can be fully reproduced: give parameters for every distribution used (check `tab:sim_params`)
- [x] ~~**L668** Polish typos and tenses in the simulation section~~ — done
- [ ] **L650** Add a note on which scenario the movement network applies to
- [ ] **L677** Clarify "scaled down by 0.10": does it mean multiplied by 0.10 (90% reduction) or reduced to 0?
- [ ] **L680** Confirm time units are days, use the unit consistently, and say the simulation uses a daily time step
- [ ] **L684** Add a reference for phastSim and say what it is used for
- [ ] **L686** Clarify "baseline sequence" (Wuhan-Hu-1): is it the sequence used to seed the outbreak?
- [ ] **L689** Clarify how 8000 sequences per simulation were subsampled (fixed sequencing probability, or exactly 8000?)
- [ ] **L697** Define the phases and which sequences were used for each one. Consider a figure showing how transmission intensity changes across scenarios
- [ ] **L698** Say which nRR was calculated and what the reference is
- [ ] **L705** Figure caption: say what the box and whiskers represent

### Additional Geographic Analyses
- [x] ~~**L719** Check whether the Oklahoma–Mexico link (BA.1 oversampling) is mentioned in the main text~~ — done: Oklahoma/BA.1 oversampling now discussed as a limitation in the main text
- [ ] **L720** Make sure all supplementary figures and analyses are referenced in the main text (see Global)

### Additional Age Analyses
- [x] ~~**L763** Reconsider whether the between-region heatmap finding (00–04y RR decreasing with other pediatric groups) is worth reporting, since little age structure is expected for non-local sequences~~ — moot: Florida/Mideast seasonal supplemental figure removed
