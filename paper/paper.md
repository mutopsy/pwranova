---
title: 'pwranova: An R package for power analysis of flexible ANOVA designs and related tests'
authors:
  - name: "Hiroyuki Muto"
    orcid:  0000-0002-0007-6019
    affiliation: "1"
affiliations:
  - name: Osaka Metropolitan University
    index: 1
date: 2025-09-28
bibliography: paper.bib
---

# Summary

Power analysis is a critical step in the design of psychological and behavioral experiments, 
yet existing tools often lack the flexibility to accommodate complex ANOVA designs. 
`pwranova` is an R package that performs power analysis for between-, within-, and mixed-factor ANOVA designs, 
with full support for main effects, interactions, and planned contrasts (custom contrasts with user-defined weights).  
The package allows researchers to calculate statistical power, 
required total sample size, significance level, or minimal detectable effect sizes 
expressed as partial eta squared or CohenÅfs *f*. 

In addition to ANOVA, `pwranova` provides complementary functions for 
common related tests, including *t*-tests (one-sample, paired, and two-sample) 
and tests of PearsonÅfs correlation (using either the *t*-distribution or Fisher's *z*-transformation approach). 
This makes the package a convenient toolkit for planning experimental studies 
in psychology and related fields.

# Statement of need

Researchers in psychology and the behavioral sciences frequently rely on 
analysis of variance (ANOVA) to analyze factorial designs with multiple 
between- and within-participant factors. However, existing tools such as G*Power [@faul2007gpower] and the `pwr` R package [@champely2020pwr] 
offer only limited flexibility when it comes to handling such complex designs. 
For example, specifying interactions in multi-factor mixed designs is difficult or not directly supported in these tools. 
They also generally do not allow direct specification of user-defined contrasts, and while effect sizes can be 
specified via CohenÅfs *f* [@cohen1988statistical], they do not directly support partial eta squared, 
which is more commonly reported in psychological research as a standard effect size index.

`pwranova` addresses these limitations by providing:  
- Support for between-, within-, and mixed-factor ANOVA designs, including both main effects and interactions.  
- Power analysis for planned contrasts with flexible, user-defined weight specification.  
- Methods based on the noncentral *F*-distribution.  
- Integrated functions for related *t*-tests and Pearson correlations.  
- A unified and extensible R implementation designed for reproducible research workflows.  

In addition, `pwranova` not only extends power analysis to complex factorial ANOVA designs but also incorporates related t-tests and correlation tests within the same framework. This integration allows researchers to conduct power analysis for a wide range of commonly used statistical tests in a consistent and reproducible way.

This combination of flexibility and reproducibility makes `pwranova` 
especially useful for experimental psychologists and cognitive scientists, as well as researchers in the behavioral, social, and biological sciences designing studies with complex factorial structures.
Detailed examples and tutorials are available on the package website:  
<https://mutopsy.github.io/pwranova/>

# Acknowledgements

This work was supported by a JSPS Grant-in-Aid for Early-Career Scientists (number 21K13750). The author used ChatGPT to obtain suggestions for English phrasing and to improve clarity and readability during the writing process. All content was critically reviewed and finalized by the author.

# References
