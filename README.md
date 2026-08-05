# RSV antibiotics use and vaccine impact in children under 5
This repository contains the code for the research: "Antibiotics use in children under five years with respiratory syncytial virus infection: a global systematic analysis and assessment of the impact of maternal vaccination"

Authors: Caini Wang; Bingbing Cong; Nana Wang; You Li

Corresponding author: You Li, You.Li@njmu.edu.cn

All the analysis and figure generation can be run from "main.R". This loads in all the dependencies and functions. All the code and data (which is publicly available) are in the repo, and can be run directly from the master script.

### IMPORTANT

- Please confirm that Rtools is installed on your computer, otherwise the code cannot run.
- Please install the necessary packages at first.
- No non-standard hardware is required.

### Input files

- Data used in this analysis are extracted from published studies and are provided in the **Supplementary Appendix** of the manuscript.
- Data/data.xlsx: All collected data extracted from included research papers.
- Location Data: Contains coordinate information for each study site. You can use the provided file directly, or generate coordinates from the source code (requires a Google Developer account and Google Maps API configuration).

### Output files
When the code runs, it will generate a file path named Output, and all running results will be saved to this path.
All analyses were completed on a Lenovo ThinkBook 14 G6+ IMH with an Intel Core Ultra 7 155H processor and 32GB of memory. The session info was as follows:

```
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  
[2] LC_CTYPE=Chinese (Simplified)_China.utf8    
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                                
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] sp_2.2-0

loaded via a namespace (and not attached):
 [1] generics_0.1.3      xml2_1.3.8          stringi_1.8.4       lattice_0.22-6     
 [5] lme4_1.1-36         hms_1.1.3           digest_0.6.37       magrittr_2.0.3     
 [9] grid_4.4.2          RColorBrewer_1.1-3  meta_8.0-2          CompQuadForm_1.4.3 
[13] Matrix_1.7-1        backports_1.5.0     forestplot_3.1.7    purrr_1.0.2        
[17] scales_1.4.0        numDeriv_2016.8-1.1 abind_1.4-8         reformulas_0.4.4   
[21] Rdpack_2.6.3        cli_3.6.3           rlang_1.1.4         rbibutils_2.3      
[25] splines_4.4.2       tools_4.4.2         tzdb_0.5.0          nloptr_2.2.1       
[29] checkmate_2.3.2     minqa_1.2.8         metafor_4.8-0       dplyr_1.1.4        
[33] ggplot2_3.5.2       mathjaxr_1.6-0      boot_1.3-31         vctrs_0.6.5        
[37] R6_2.5.1            lifecycle_1.0.4     stringr_1.5.1       MASS_7.3-61        
[41] pkgconfig_2.0.3     pillar_1.10.1       gtable_0.3.6        glue_1.8.0         
[45] Rcpp_1.0.13-1       tibble_3.2.1        tidyselect_1.2.1    rstudioapi_0.17.1  
[49] dichromat_2.0-0.1   farver_2.1.2        nlme_3.1-166        readr_2.1.5        
[53] compiler_4.4.2      metadat_1.4-0      
```   

See LICENSE file for licensing details. Code will be uploaded by August 7.
