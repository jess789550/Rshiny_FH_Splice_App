# Splice apps
These apps allows the user to view splice site prediction results from SpliceAI, MaxEntScan, MMSplice, GeneSplicer, and SQUIRLS. The data must originate from the OUH FH or Exome pipeline. The FH splice app demonstrates the FH validation data and the Exome splice app demonstrates the Exome validation data. The generic splice app allows the user to upload their own CSV file with columns: 

- file_id
- CHROM
- POS
- REF
- ALT
- SYMBOL
- HGVSc
- gnomAD_AF
- SpliceAI_DS_AG
- SpliceAI_DS_AL
- SpliceAI_DS_DG
- SpliceAI_DS_DL
- mmsplice_delta_logit_psi
- MaxEntScan_alt
- MaxEntScan_diff
- GeneSplicer_score
- SQUIRLS
- Type
- QUAL


The user can adjust cutoff thresholds for each tool and the app will show how this affects the false discovery rate.

Default settings are the optimum values for OUH FH and Exome data.

## FH app worklists

Results for the following worklists can be viewed:
- 2004442 (32 samples containing controls and GIAB samples, 3 positive)
- 2005265 (similar to above)
- 2005267 (similar to above)
- 2005745 (similar to above)
- 2324533 (64 samples, 2 positive)
- 2400720 (64 samples, 2 positive)
- 2322015 (64 samples, 2 positive)
- 2127291 (64 samples, 2 positive)
- 2226732 (64 samples, 2 positive)
- 2327211 (64 samples, 1 positive)
- 2330804 (64 samples, 1 positive)
- 2331473 (64 samples, 1 positive)
- mix (64 samples, 19 positive)

## Exome app worklists

Results for the following validation worklists:
- 1234567
- 1234568
- 1234569
- 1234570
- 1234571
- 1234572
- 1234573
- 1234574

Each worklist contains 4 positive samples, 59 negative samples, and 1 GIAB sample.

## Cutoffs and thresholds

The SpliceAI delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), ranges from 0 to 1 and can be interpreted as the probability of the variant being splice-altering. In the scores table, 0.2-0.5 is colored red, 0.5-0.8 is colored amber, and 0.8-1 is colored green.

MaxEntScan (MES) outputs a diff and alt score which are used to assess the potential of disrupting a native splice site. In the scores table, Low potential is colored red, Medium is colored amber, and High potential is colored green. The MES cutoff filters are:
- None = Unfiltered
- Low = potential of disrupting native splice sites is diff > 0 and alt < 8.5; or diff < 0 and alt > 6.2
- High = potential of disrupting native splice sites is diff > 0 and alt < 6.2; or diff < 0 and alt > 8.5

The SQUIRLS score denotes a probability of the variant being splice-altering ranging from 0 to 1. In the scores table, 0-0.5 is colored red, 0.5-0.9 is colored amber, and 0.9-1 is colored green.

Delta_logit_psi is the main score predicted by MMSplice, which shows the effect of the variant on the inclusion level (PSI percent spliced in) of the exon. 0 is the most lenient threshold and 2 is the strictest. In the scores table, 0.5 < |delta_logit_psi| < 1 is colored red, 1 < |delta_logit_psi| < 2 is colored amber, and |delta_logit_psi| > 2 is colored green.

Although GeneSplicer has low recall, it has been added at a clinical scientist's request. A context of 100 base pairs (default) was used, which means it is the amount of sequence included either side of the variant. In the scores table, the absence of a score is colored red, and the presence of a score is colored green.

Other tools considered were either too slow (Spliceator, DSSP, Pangolin, DNABERT) or were deprecated (KipoiSplice, PresPSI-SVR, Splice2Deep, S-CAP).

The gnomAD allele frequency can be filtered to less than or equal to 0, 0.0001, 0.0002 (PM2), 0.002 (BS1), 0.005 (BA1), or None (no filter).


## How to use
1. Choose/upload CSV file.

2. Use sliders to filter data then press "Submit". 

3. To filter the splice variants and prediction scores table Please press "Filter" after using splice variant table column filters.

4. Navigate through the different tabs to see the results.



## Run app locally
```
cd <path_to_app>
R
library(shiny)
runApp("splice")
```


## Deployment command example
```
R
library(rsconnect)
rsconnect::deployApp('/data/jess_tmp/fh/Rshiny/fh_splice') 
```


## View apps online
https://jessicakan.shinyapps.io/fh_splice/ 
https://jessicakan.shinyapps.io/exome_splice/ 
https://jessicakan.shinyapps.io/splice/ 




## Author
Jessica Kan



## References
Jaganathan, K., Kyriazopoulou Panagiotopoulou, S., McRae, J. F., Darbandi, S. F., Knowles, D., Li, Y. I., Kosmicki, J. A., Arbelaez, J., Cui, W., Schwartz, G. B., Chow, E. D., Kanterakis, E., Gao, H., Kia, A., Batzoglou, S., Sanders, S. J. and Farh, K. K. (2019) 'Predicting Splicing from Primary Sequence with Deep Learning', Cell, 176(3), pp. 535-548.e24.

Yeo, G. and Burge, C. B. (2004) 'Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals', J Comput Biol, 11(2-3), pp. 377-94.

Cheng, J., Nguyen, T. Y. D., Cygan, K. J., Çelik, M. H., Fairbrother, W. G., Avsec, Ž. and Gagneur, J. (2019) 'MMSplice: modular modeling improves the predictions of genetic variant effects on splicing', Genome Biol, 20(1), pp. 48.

Danis, D., Jacobsen, J. O. B., Carmody, L. C., Gargano, M. A., McMurry, J. A., Hegde, A., Haendel, M. A., Valentini, G., Smedley, D. and Robinson, P. N. (2021) 'Interpretable prioritization of splice variants in diagnostic next-generation sequencing', Am J Hum Genet, 108(11), pp. 2205.

Shamsani, J.,, Kazakoff, S.H.,, Armean, I.M.,, McLaren, W.,, Parsons, M.T.,, Thompson, B.A.,, O’Mara, T.A.,, Hunt, S.E.,, Waddell, N., and Spurdle, A.B., 2018. A plugin for the Ensembl variant effect predictor that uses MaxEntScan to predict variant spliceogenicity. Bioinformatics, 35(13), pp.2315–2317. 

Pertea, M., Lin, X. and Salzberg, S. L. (2001) 'GeneSplicer: a new computational method for splice site prediction', Nucleic Acids Res, 29(5), pp. 1185-90.
