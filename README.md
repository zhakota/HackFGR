# Auxological predictors of intrauterine fetal growth restriction (FGR)
The repository with working scripts, figures, and raw data used for testing auxological FGR prediction <br>
<br/><br/>
Team Lead: Dmitrii Zhakota

Team Members:Ruslan Alagov, Yury Malovichko, Anton Shikov
<br/><br/>
<br/><br/>
Bioinformatics Hackathon'2023


## Project objectives
<ol>
<li>Estimate the level of comparability between obstetric, neonatal, and pathological evidence for fetal growth dynamics</li>
<li>Compare an actual dataset of newborn autopsy records to the external evidence</li>
<li>Assess the adequacy of the existing approaches to FGR diagnostics and suggest alternative approaches</li>
</ol>

## Data
Analyzed raw data are included in the `data` directory.

* `data/raw/` folder contains raw data from existing datasets and experimental data used for statistical analysis.
* `data/output/` folder contains the results of gamma-regression applied on the Intergrowth21 dataset regarding body length and weight for boys, girls, or combined weight data.

## Scripts
The `src/` directory includes all code used for comparing methods, imputing, and calculating feature importance.
* `src/HackFGR.Rproj` file determines basic options within R scripts;
* `src/to_15.ipynb` script builds regressions based on Intergrow data. Data are available for different characteristics (weight, height, weight-for-height ratio) and different genders (boys, girls, combined). In each dataset, there are averages for fetal ages 24-42 weeks. For each characteristic, a regression is built for each gender, allowing us to recreate the picture for up to 15 weeks. Then, graphs are built for each regression, and a comparison is made with the data on Kiserud
*  `src/Data_preparation.R` script is used for comparing differences between genders for each time point of gestation weeks as well as between methods for assessing the body weight of a fetus using ultrasound, pathology analysis, and neonatological observations. In the script, comparisons between genders are calculated using pair-wise t-tests for a certain gestation week. Both raw and FDR-adjusted (false discovery rate) p-values are presented. For comparing diagnostic methods, three approaches were applied, namely, pair-wise t-test (with or without correction for multiple comparisons), meta-analysis (umbrella analysis with uniting datasets belonging to the same diagnostic method and one-sample mean meta-analysis per study), and time-series Chow test.
*  `src/FGR_filter.Rmd` script includes approaches to both cleaning the initial source dataset from outliers and calculating metrics used for diagnosing FGR based on body length and weight of certain organs.

## Plots
The `pics/` directory contains graphical representations of the obtained results. 

## Conclusions
<ol>
<li>Auxological observations from pathologists, neonatologists, and obstetricians correspond to each other and can be plausibly compared</li>
<li>Several assumed patterns and approaches do not seem to hold upon multiple sources’ comparison</li>
  <ol>
<li>The assumption of lower weight-by-gestation the pathological evidence is only partially supported </li>
<li>Despite a long-running tradition, gender does not affect general auxological parameters</li>
    </ol>
<li>Advanced auxological indices may serve as plausible FGR predictors (needs further testing)</li>
</ol>
