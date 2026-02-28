---
title: 'pwranova: An R package for power analysis of flexible ANOVA designs and related tests'
authors:
  - name: "Hiroyuki Muto"
    orcid:  0000-0002-0007-6019
    affiliation: "1"
affiliations:
  - name: Osaka Metropolitan University
    index: 1
date: 2025-10-14
bibliography: paper.bib
---

# Summary

Power analysis is a critical step in the design of psychological and behavioral experiments, 
yet existing tools often lack the flexibility to accommodate complex ANOVA designs. 
`pwranova` is an R package that performs power analysis for between-, within-, and mixed-factor ANOVA designs, 
with full support for main effects, interactions, and planned contrasts (custom contrasts with user-defined weights).  
The package allows researchers to calculate statistical power, 
required total sample size, significance level, or minimal detectable effect sizes 
expressed as partial eta squared ($\eta^2_p$) or Cohen's *f* [@cohen1988statistical]. 

In addition to ANOVA, `pwranova` provides complementary functions for 
common related tests, including *t*-tests (one-sample, paired, and two-sample) 
and tests of Pearson correlation (using either the *t*-distribution or Fisher's *z*-transformation approach). 
This makes the package a convenient toolkit for planning experimental studies 
in psychology and related fields.

# Statement of need

Researchers in psychology and the behavioral sciences frequently rely on 
analysis of variance (ANOVA) to analyze factorial designs with multiple 
between- and within-participant factors. However, existing tools such as G\*Power [@faul2007gpower] and the `pwr` R package [@champely2020pwr] 
offer only limited flexibility when it comes to handling such complex designs. 
For example, specifying interactions in multi-factor mixed designs is difficult or not directly supported in these tools. 
They also generally do not allow direct specification of user-defined contrasts. In addition, although effect sizes 
can be specified via Cohen's f [@cohen1988statistical], they do not directly support partial eta squared, 
which is more commonly reported in psychological research.

`pwranova` addresses these limitations by providing:  
- Support for between-, within-, and mixed-factor ANOVA designs, including both main effects and interactions.  
- Power analysis for planned contrasts with flexible, user-defined weight specification.  
- Methods based on standard noncentral *F*-distribution power calculations.  
- Integrated functions for related *t*-tests and Pearson correlations.  
- A unified and extensible R implementation designed for reproducible research workflows.  

In addition, `pwranova` not only extends power analysis to complex factorial ANOVA designs but also incorporates related *t*-tests and correlation tests within the same framework. This integration allows researchers to conduct power analysis for a wide range of commonly used statistical tests in a consistent and reproducible way.

This combination of flexibility and reproducibility makes `pwranova` 
especially useful for experimental psychologists and cognitive scientists, as well as researchers in the behavioral, social, and biological sciences designing studies with complex factorial structures.

# Example use case

To illustrate a realistic application, consider a visual search experiment investigating age-related differences in search efficiency. In this type of task, participants search for a target item among distractors on a screen. They indicate whether the target is present or absent by pressing a key, and response time is recorded. The study includes one between-participant factor, age group (young vs. older adults), and two within-participant factors: target presence (present vs. absent) and set size (8, 16, or 24 items displayed), forming a 2 $\times$ 2 $\times$ 3 mixed design. Suppose the researcher wishes to detect a three-way interaction effect size of *f* = 0.25 with 80% power at $\alpha = .05$. The required total sample size can be estimated as follows:

```r
pwranova(
  nlevels_b = 2,
  nlevels_w = c(2, 3),
  cohensf   = 0.25,
  alpha     = 0.05,
  power     = 0.80,
  target    = "B1:W1:W2"
)
```

The output indicates that a total sample size of 156 is required to achieve the desired statistical power (assuming equal group sizes). 
If the `target` argument is omitted, the function returns results for all main effects and interactions under the specified design. 
Although Cohen's *f* is used here, $\eta^2_p$ (`peta2`) can also be specified.a 
More detailed examples and tutorials are available on the package website: 
<https://mutopsy.github.io/pwranova/>

# Acknowledgements

This work was supported by a JSPS Grant-in-Aid for Early-Career Scientists (number 21K13750). The author used ChatGPT to obtain suggestions for English phrasing and to improve clarity and readability during the writing process. All content was critically reviewed and finalized by the author.

# References
