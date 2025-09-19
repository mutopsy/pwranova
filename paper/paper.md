---
title: 'pwranova: An R package for power analysis for flexible ANOVA designs and related tests'
authors:
  - name: "Hiroyuki Muto"
    orcid:  0000-0002-0007-6019
    affiliation: "1"
affiliations:
  - name: Osaka Metropolitan University
    index: 1
date: 2025-09-20
bibliography: paper.bib
---

# Summary

Power analysis is a critical step in the design of psychological and behavioral experiments, 
yet existing tools often provide limited flexibility in handling complex ANOVA designs. 
*pwranova* is an R package that implements power analysis for between-, within-, and mixed-factor ANOVA designs, 
with full support for main effects, interactions, and planned contrasts (custom contrasts with user-defined weights).  
The package allows researchers to calculate statistical power, 
required total sample size, significance level, or minimal detectable effect sizes 
expressed as partial eta squared or CohenÅfs *f*. 

In addition to ANOVA terms, *pwranova* provides complementary functions for 
common related tests, including *t*-tests (one-sample, paired, and two-sample) 
and tests of PearsonÅfs correlation. 
This makes the package a convenient toolkit for planning experimental studies 
in psychology and related fields.

# Statement of need

Researchers in psychology and the behavioral sciences frequently rely on ANOVA 
to analyze factorial designs with multiple between- and within-subject factors. 
However, existing tools such as G\*Power have limited flexibility 
for complex factorial structures and do not always allow direct specification of contrasts 
with user-defined weights, or minimal detectable effect sizes in terms of partial eta squared or CohenÅfs *f*. 

*pwranova* addresses these limitations by providing:
- Support for arbitrary factorial ANOVA terms, including interactions.  
- Planned contrast power analysis with flexible weight specification.  
- Exact and approximate methods for power analysis using the noncentral *F* distribution.  
- Integrated functions for related *t*-tests and Pearson correlation.  
- A unified and extensible R implementation suitable for reproducible workflows.  

This combination of flexibility and reproducibility makes *pwranova* 
especially useful for experimental psychologists, 
cognitive scientists, and other researchers designing studies 
with complex factorial structures.

# Usage

Detailed examples, tutorials, and vignettes are available on the package website:  
<https://mutopsy.github.io/pwranova/>

# Acknowledgements

This work was supported by a JSPS Grant-in-Aid for Early-Career Scientists (number 21K13750). The author used ChatGPT to obtain suggestions for English phrasing and to improve clarity and readability during the writing process. All content was critically reviewed and finalized by the author.

# References
