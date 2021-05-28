# ImprovingCloningEfficacy
This code is for finding small molecules that can increase Cloning efficacy:


Please follow the following steps before running the code:
- Download gsea software from the following link:
http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/resources/software/gsea2-2.2.4.jar

- Place the downloaded jar file in the gsea folder.

- java 8 installation is a requirement.

- java should be added to the environmental variables in windows or added to the Path in linux.

- Run Cloning.r

- After running the R file, the rank ordered list of compounds would be generated in the results folder. One file takes into account both IVF embryo samples (SCNT2IVF) and one file takes into account one IVF embryo samples (SCNT2IVF1) . Some detailed gsea output would be generated in the output directory which can be ignored. Best compound candidate names would be printed.

- Refer to following link, under “Preparing Data Files for GSEA”, for description of file types in the data folder that are input to gsea software (gct, cls, chip).
https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Other_Ways_to

