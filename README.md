# FH splice app Description
This app allows the user to view splice site prediction results from SpliceAI, MaxEntScan, GeneSplicer, MMSplice, Pangolin, and SQUIRLS.

The user can adjust cutoff thresholds for each tool (future work: show how this affects the precision and recall rate of each tool).

Results for the following worklists can be viewed:
2004442 (32 samples containing controls and GIAB samples, 3 positive)
2005265 (similar to above)
2005267 (similar to above)
2005745 (similar to above)
2324533 (64 samples, 2 positive)
2400720 (64 samples, 2 positive)
2322015 (64 samples, 2 positive)
2127291 (64 samples, 2 positive)
2226732 (64 samples, 2 positive)
2327211 (64 samples, 1 positive)
2330804 (64 samples, 1 positive)
2331473 (64 samples, 1 positive)
mix (64 samples, 19 positive)

# Run tool
```
cd /data/jess_tmp/fh/Rshiny
R
library(shiny)
runApp("splice")
```

# Author
Jessica Kan
