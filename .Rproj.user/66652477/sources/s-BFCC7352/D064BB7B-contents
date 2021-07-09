# mibcCPIclass - Classify MIBC samples according to CPI-response related classes

Install the package directly from Github:
```r
# install.packages("remotes")
remotes::install_github("csgroen/mibcCPIclass")
```
For classification, a RNA-seq gene expression matrix normalized with log2(FPKM)
is recommended, with samples on columns and genes on rows. Supported IDs are HGNC
symbol, ENTREZ id and Ensembl gene IDs for the classification genes.

## Example

We include some example data from the TCGA-BLCA cohort:
```r
library(mibcCPIclass)
data(tcga_blca_ex)
classify_mibcCPIclass(tcga_blca_ex)
# A tibble: 81 x 8
   id         prediction class_name probability_s1 probability_s2 probability_s3 probability_s4 probability_s5
   <chr>      <fct>      <fct>               <dbl>          <dbl>          <dbl>          <dbl>          <dbl>
 1 TCGA-5N-A… S4         MycU           0.0959       0.262              0.195            0.447      0.000400  
 2 TCGA-GD-A… S4         MycU           0.243        0.0364             0.00851          0.712      0.0000648 
 3 TCGA-G2-A… S1         LumE           0.767        0.0658             0.0192           0.144      0.00377   
 4 TCGA-BT-A… S4         MycU           0.00371      0.000743           0.256            0.740      0.0000706 
 5 TCGA-FD-A… S2         LumR           0.235        0.310              0.264            0.189      0.00162   
 6 TCGA-DK-A… S1         LumE           0.602        0.157              0.00122          0.239      0.0000204 
 7 TCGA-DK-A… S5         StroR          0.00205      0.00800            0.00922          0.280      0.701     
 8 TCGA-XF-A… S1         LumE           0.977        0.00677            0.000474         0.0156     0.00000343
 9 TCGA-C4-A… S3         ImmBas         0.00000267   0.0000000627       0.973            0.0269     0.00000400
10 TCGA-2F-A… S2         LumR           0.142        0.482              0.0614           0.310      0.00528   
# … with 71 more rows
```

If you're having problems, feel free to open an issue here on Github.
