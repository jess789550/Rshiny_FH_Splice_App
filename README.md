# FH splice app Description
This app allows the user to view splice site prediction results from SpliceAI, MaxEntScan, MMSplice, and SQUIRLS.

The user can adjust cutoff thresholds for each tool and the app will show how this affects the false discovery rate.

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

The SpliceAI delta score of a variant, defined as the maximum of (DS_AG, DS_AL, DS_DG, DS_DL), ranges from 0 to 1 and can be interpreted as the probability of the variant being splice-altering. 

MES cutoff:
- Low potential of disrupting native splice sites is diff > 0 and alt < 8.5; or diff < 0 and alt > 6.2
- High potential of disrupting native splice sites is diff > 0 and alt < 6.2; or diff < 0 and alt > 8.5

The SQUIRLS score denotes a probability of the variant being splice-altering ranging from 0 to 1.

Delta_logit_psi is the main score is predicted by MMSplice, which shows the effect of the variant on the inclusion level (PSI percent spliced in) of the exon. 0 is the most lenient threshold and 2 is the strictest.

Other tools considered were either of low recall (GeneSplicer), too slow (Spliceator, DSSP, Pangolin, DNABERT) or were deprecated (KipoiSplice, PresPSI-SVR, Splice2Deep, S-CAP).



# How to use
Use sliders to filter data then press "submit".

Navigate through the different tabs to see the results.



# Run app locally
```
cd <path_to_app>
R
library(shiny)
runApp("splice")
```



# View app online
https://jessicakan.shinyapps.io/splice/ 




# Author
Jessica Kan



# References
Jaganathan, K., Kyriazopoulou Panagiotopoulou, S., McRae, J. F., Darbandi, S. F., Knowles, D., Li, Y. I., Kosmicki, J. A., Arbelaez, J., Cui, W., Schwartz, G. B., Chow, E. D., Kanterakis, E., Gao, H., Kia, A., Batzoglou, S., Sanders, S. J. and Farh, K. K. (2019) 'Predicting Splicing from Primary Sequence with Deep Learning', Cell, 176(3), pp. 535-548.e24.

Yeo, G. and Burge, C. B. (2004) 'Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals', J Comput Biol, 11(2-3), pp. 377-94.

Cheng, J., Nguyen, T. Y. D., Cygan, K. J., Çelik, M. H., Fairbrother, W. G., Avsec, Ž. and Gagneur, J. (2019) 'MMSplice: modular modeling improves the predictions of genetic variant effects on splicing', Genome Biol, 20(1), pp. 48.

Danis, D., Jacobsen, J. O. B., Carmody, L. C., Gargano, M. A., McMurry, J. A., Hegde, A., Haendel, M. A., Valentini, G., Smedley, D. and Robinson, P. N. (2021) 'Interpretable prioritization of splice variants in diagnostic next-generation sequencing', Am J Hum Genet, 108(11), pp. 2205.

Shamsani, J.,, Kazakoff, S.H.,, Armean, I.M.,, McLaren, W.,, Parsons, M.T.,, Thompson, B.A.,, O’Mara, T.A.,, Hunt, S.E.,, Waddell, N., and Spurdle, A.B., 2018. A plugin for the Ensembl variant effect predictor that uses MaxEntScan to predict variant spliceogenicity. Bioinformatics, 35(13), pp.2315–2317. 
