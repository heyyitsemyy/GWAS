# GWAS
This is meant for inspiration. Copying is at your own discretion.

NOTE: PLINK binary files were too large for upload. Only code is uploaded

PROJECT DESCRIPTION
1) Perform linear regression in PLINK
2) Perform linear regression using R with given snpdata and phenotype data
3) Construct a manhattan plot with the PLINK linear regression data
4) Make a table of the minimum p-values in the PLINK linear regression data and compare with data from the UK bio-bank 
	* To download UK bio-bank data, go to this site, the Manifest 201807 tab, and look for HDL cholesterol (quantile) for both_sexes: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291 

------------------------------------------------------------------------------------
PLINK Commands:
To generate the linear regression data from PLINK:
  plink --bfile snpdata --linear --pheno phenotype.txt
	  output file name: plink.assoc.linear
To generate the transposed snpdata file:
  plink --bfile snpdata --recode 
	  output file name: plink.tped

Project_Script.R:
This R script contains all the code on how I implemented each part of the project.

Part 2: Linear Regression
  Required Libraries: data.table
  Required Files: phenotype.txt, plink.tped (transposed snpdata binary file)
  All p-values are stored in the object lin_outputs
	lin_outputs length, head, and tail are printed to show accuracy

Part 3: Manhattan Plot
  The given script qqman.r is copied and pasted
  Required Libraries: data.table
  plink.assoc.linear is assigned to plink_linear
  The manhattan function is called on plink_linear after data refining

Part 4: Minimum P-values
  Required Libraries: data.table, R.utils
  Required Files: plink.assoc.linear, 30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
        **30760_irnt.gwas.imputed_v3.both_sexes.tsv.bgz is unzipped into UK.txt**
  Table is stored in data frame snpdata
