# Auxological predictors of intrauterine fetal growth restriction (FGR)
The repository with working scripts, figures, and raw data used for testing auxological FGR prediction <br>
<br/><br/>
Team Lead: Dmitrii Zhakota

Team Members:
Ruslan Alagov
Yury Malovichko
Anton Shikov
<br/><br/>
<br/><br/>
Bioinformatics Hackathon'2023


## Project objectives
1. Estimate the level of comparability between obstetric, neonatal, and pathological evidence for fetal growth dynamics
2. Compare an actual dataset of newborn autopsy records to the external evidence
3. Assess the adequacy of the existing approaches to FGR diagnostics
4. Test whether auxological FGR prediction can be enhanced with novel statistical approaches


## Data
Analyzed raw data are included in the `data` directory.

* `data/raw/` folder contains raw data from existing datasets and experimental data used for statistical analysis.
* `data/output/` folder contains the results of gamma-regression applied on the Intergrowth21 dataset regarding body length and weight for boys, girls, or combined weight data.

## Scripts
The `src/` directory includes all code used for comparing methods, imputing, and calculating feature importance.
* `src/HackFGR.Rproj` folder contains raw data from existing datasets and experimental data used for statistical analysis.
* `src/output/` folder contains the results of gamma-regression applied on the Intergrowth21 dataset regarding body length and weight for boys, girls, or combined weight data.

## Plots
The `pics/` directory contains graphical representation of the obtained results. 

