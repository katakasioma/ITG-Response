# ITG-Response

## Description/Instructions

Please find below a few questions that mimic some common problems we encounter at ITG. They are grouped by broad theme. You will notice there are many questions, the goal is not to answer them all but to pick a few questions to focus on (10 is a good number, but pick as many as you want). You should pick from all three categories, but there are many more bioinformatics questions so you should mainly pick from those. We encourage you to choose your questions according to your areas of expertise but also to try and answer questions that are as varied as possible.

For programmatic questions, you can use the language and libraries of your choice, but we will assess whether your choice of language was optimal. Try and aim for a minimal solution in terms of code length. If you use a shell script, you can assume that common non-core packages will be installed (e.g. `awk`, `sed`, `perl`, `python`, `sponge`, `wget` or `jq`). You can use the shell of your choice, if not otherwise specified we will assume `bash`. Assume that all common bioinformatics tools `bcftools`, `bedtools`, `vcftools`, `plink` and others are all installed.

We are primarily interested in how you would solve these problems if you encountered them in real life. Whenever the command line or programming is used, please include your code along with your answer. Not all questions are programmatic, and some are open-ended. Feel free to include details and to discuss potential issues if you don't think there is a clear-cut answer.

To submit your results, please clone this repository and make your edits. Once you're done, send us a link to your repository, or compress it and send us a link to the archive.

## Questions

### Support/resource management/Shell
1. A user has several versions of R installed in their path. Each version of R has a number of locally installed libraries. The user is confused, and would like to know which library is installed for each version of R. Can you write a command to help them out?.
  ```
  ##Assuming uer has a linux OS, Start R from commandline, check R versions installed, check installed packages
    R
    version
    installed.packages()
  ```
    
 ###   
2. A common problem with shared filesystems is a disk quota overflow. This can be due to 1) a large combined size on disk or 2) too many files present on disk. We would like to help users who encounter this problem to locate the problematic files/directories. Write a command to sort all subdirectories of level `n` (`n` determined by the user) by their human-readable size. Write another command to sort all subdirectories of level `n` according to the number of files they contain.

  ```
  # assuming n is 100GB, we can list all files => 100GB
    find . -type f -size +100G
  ```
    
    
3. A user wants to install an `R` package and gets the following [error log](data/error.log). What is likely to cause the error and how can they solve it?

  ```bash
    I will first check the g++ version and then change the -std=c++11 to one that is compatible to the g++ version 
  ```
    
4. A user is running commands like this one `cat file1 <(cut -d " " -f 1-15,17,18 file2) > file3`. What does this command do? It runs fine on the command line, but then the user includes it into a file with other commands, saves it and runs `chmod +x` on it. However, that line of code throws the following error : `syntax error near unexpected token '('`. What has the user forgotten?

  ```
    This command takes contents of columns 1-15, 17, 18 of file 2 and merges them to file1 contents, finally it saves the merged data in file 3
    To make the file excecutable, the user must escape the brackets
  ```
    
5. A collaborator has sent you [this script](data/EasyQCWrapper.sh). It is a wrapper for a bioinformatics software called `EasyQC`.  Running it, you get the following error: 

    ```bash
    ./test.EasyQC-START.R: line 6: syntax error near unexpected token 'EasyQC'
    ./test.EasyQC-START.R: line 6: 'library(EasyQC)'
    ```

     You need to run this script now, but your collaborator is unavailable for a few days. What is causing the error? (Hint: Nothing is wrong with the `.ecf` EasyQC script.)
     
    ```
    The user forgot to uncomment the variables $1-$4 thus they cannot be accessed
    ```
    
6. Programmatic download
    - You have to download all autosomal files from this location: [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/) onto **your server**. You connect to the server via SSH. Using only the command line, how do you perform this download?
    
    ```
    Step 1. SSH to server
    Step 2 . Data download using the command :
    wget -m -e robots=off --no-parent http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
    ```
- You are at a conference abroad and you quickly realise that your connection is unstable. You get disconnected constantly, which interrupts the download. How do you ensure the download survives these disconnections?
    
```
Step1. Create a job ID in screen (in the background)
Step2. Run download
wget -m -e robots=off --no-parent http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/
```
    
7. Bioinformaticians often work on a computing cluster. The cluster runs a software called a job scheduler that attributes resources to users depending on the requirements of their jobs. In this case, let's imagine the cluster is running IBM LSF. You do not need to know it to answer this question. The `bjobs` command lists all jobs submitted by the user (manpage [here](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bjobs.1.html)). It gives this kind of output:
    ```
    JOBID   USER             STAT  QUEUE      FROM_HOST EXEC_HOST JOB_NAME SUBMIT_TIME
    9670681 current_user     RUN   basement   head_node node1     job1     Oct 24 10:24
    9740051 current_user     RUN   basement   head_node node1     job2     Oct 24 17:41
    9670681 current_user     RUN   normal     head_node node2     job3     Oct 24 10:24
    9740981 current_user     PEND  basement   head_node           job4     Oct 24 17:44

    ```
     - Given the [following output](data/farm-snapshot.txt) of `bjobs -all`, which users are the top 5 users of the cluster?
     - How many jobs does the user `pathpip` have running in all queues?
     - A user wants to know how many jobs they have pending (`PEND`) and running (`RUN`) in each queue. Write a command line to do that (You can use the log above to check your command line). How would they display this on their screen permanently, in real time?
     
    ```bash
    Top 5 users:
    awk -F" " '{print $2}' farm-snapshot.txt | uniq | sort | uniq -c  | sort | tail -5
    153 so11
    201 km18
    409 ro4
    418 igs
    489 pathpip

    Running jobs fro pathpip:
    cat farm-snapshot.txt | grep "pathpip" | grep "RUN" | wc -l 
    1138

   Assuming user kk8:
   awk '!seen[$2]++ {print "Feature",$2} {print $2,$1,$3,$4}' farm-snapshot.txt | sort | grep "kk8" | wc -l
   124
   ```
   
   ```
   To display On a screen user should just type 'bjobs' to see his running and pending jobs. If user is kkb
   bjobs | grep "kkb"
   ```
    
8. An analysis you need to run on the cluster requires a particular python library, but you do not have administrator rights. IT is on holiday. What do you do?

  ```
  I install the python library locally and give the path of the library when running job on cluster

  ```
    
9. All major computational tasks in your lab are done via SSH connection to mainframe servers or HPC clusters. A user comes from a Linux (mostly command-line) background but IT only support Windows 10 for laptops. How would you advise them to configure their laptop to make their transition easier?

  ```
  I install a virtual machine and then Windows OS on their linux laptop to show them how to do the jobs in Windows.
  ```

### Bioinformatics
1. The [VCF format](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/) is a popular format to describe genetic variations in a study group. It is often used in sequencing projects. Due to size concerns, it is often compressed using `gzip` and indexed using `tabix`. A binary version, BCF, also exists.
    - Write a command or script to remove duplicate positions in a VCF such as [this one](data/duplicates.vcf.gz), independently of their alleles. The positions can be duplicated an arbitrary number of times. Write code to keep the first, last and a random record among each set of duplicated records.
    
  ```
  bcftools norm -D ./duplicates.vcf.gz -o Non_duplicates.vcf.gz
  ```
 
 - Same question, but make duplicate detection allele-specific. When it finds such an exact duplicate, your code should remove all of the corresponding records.
    
 ```
 bcftools view -v snps ALL.chr21_GRCh38_sites.20170504.vcf.gz | grep -v "^#" | cut -f2 | sort -u > Non_duplicated.vcf
 ```
  
2. From an existing VCF with an arbitrary number of samples, how do you produce a VCF file without any samples using `bcftools`?

  ```
  ```
 
  
3. You are the curator of a genotype dataset with a very strict privacy policy in place. In particular, it should be impossible to tell, given access to a person's genetic data, whether they were part of your study by looking at a dataset you provided. A collaborator is asking you for some data to run tests on their code. What information can you safely contribute from your study?

  ```
  I will anonymise all the patient names and then provide the data to the collaborator
  ```
  
4. How do you convert a gzipped VCF to the `bimbam` format? (you may choose to script a solution yourself, or not)

  ```
  plink --vcf duplicates.vcf --recode bimbam --out duplicates.bimbam
  ```
  
5. A user sends you a small number of chromosome and positions in build 38 that they want to know the rsID of. 
    - What is missing from their request? What kind of unexpected output can they expect?
    - Given [this file](data/rand.chrpos.txt), honour their request using the Ensembl REST API.
    - Do the same, but offline, using the dbSNP r.150 VCF file.
    - What would change if these positions were in build 37?
    - If the user sends you 7,000 such chromosome positions, how would the above methods perform? Do you know of any alternatives?
    
  ```
  The user has neither provided the starnd information nor the chromosaomal start and end positions.
  
 Using Biomart in R we can query the Ensembl REST API
 coords=read.csv("rand.chrpos.txt", stringsAsFactors = F, strip.white = T, header = F, sep="\t")
 colnames(coords) <- c("chr","start")
 coords$stop <- coords$start
 library(biomaRt)
 ensembl=useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
 coords$chr=gsub("chr","",coords$chr)
 
 results=data.frame(matrix(ncol=5,nrow=0))

 for (i in seq (1,nrow(coords))){
  ens=getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), 
            filters = c('chr_name','start','end'), values = as.list(coords[i,]), mart = mart)
  results=rbind(results,ens)
  }
  
  
  To find the rsIDs offline: download vcf file and search iteratively for matching coordinates
  
  $ wget -qO- https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/common_all_20170710.vcf.gz | gunzip -c | convert2bed --input=vcf --output=bed --sort-tmpdir=${PWD} - > Gr38.snp150.bed
  awk 'FNR==NR{a[$1,$2,$3];next} (($1, $2,$3) in a)' coords Gr38.snp150.bed > output.txt
  
 In the case of 7000 cases, the above approach will still work, albeit a bit longer.
  ```
  
6. How would you change the chromosome numbers in the file above to chromosome names (e.g. "chr1" instead of "1")?
  ```
  I would apend the "chr" character before the chromosomal number on the particular column.
  ```
- How would you change the names back to the original? Would your solution work if an additional column containing text of arbitrary length and content is appended at the left of the file?
    
  ```
  I would delete the "chr" character before the chromosomal number. This will work since i know which column contains the chromosome numbers.
  ```
- These positions are extracted from a VCF. Convert this file to the BED format.
  
```
vcf2bed --keep-header < foo.vcf
```
    
7.	Download the 1000 Genomes sites VCF file for chromosome 21 [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38_sites.20170504.vcf.gz). We want to compare it to [a locally stored file](data/ALL.chr21_GRCh38_sites.20170504.vcf.gz).
    - What is the fastest way to check the integrity of, or compare, any such downloaded file?
```
bcftools stats ALL.chr21_GRCh38_sites.20170504.vcf.gz  ~/ITG_Responses/challenge/data/ALL.chr21_GRCh38_sites.20170504.vcf.gz > comparison    ```
```
- If you find that the files are indeed different, how do you find their differences? Keep in mind that this kind of file can be very large (>100Gb), your solution should be as fast and memory-efficient as possible.
    
    ```
    vcftools --gzvcf ALL.chr21_GRCh38_sites.20170504.vcf.gz --gzdiff /home/data/ALL.chr21_GRCh38_sites.20170504.vcf.gz --diff-site --out in1_v_in2
    ```
- If you found no differences in the end, what could cause a false alarm?
    
    ```
    ```
      
8.	What is the p-value corresponding to standard normal z-scores of 10.35, 29.7, 45.688 and 78.1479?

    ```
    ```
    
9.	We want to round a column of numbers to `n` decimal places, with values with 5 as their rightmost significant digit rounded up. Use the language of your choice.

``` 
Assuming a dataframe called roundoff with column zip as interest, we can roundoff any digit ending with 5 to 1 decimal place
roundoff$zip <- as.numeric(roundoff$zip)
roundoff$zip <- ifelse(grepl("5$", roundoff$zip), round(roundoff$zip, 1), 
                             roundoff$zip)
```

10.  Is [this HRC-imputed file](https://drive.google.com/open?id=1dOYlwIlAyz9-i4gVy2sgpQv_4YX9Y5nU) missing any chromosomes? Try to find out in seconds if you can.

``` 
cut -d$'\t' -f2  hrc.positions.txt.bgz | sort -u 

```
11.  Find out the coordinates of the _ADIPOQ_ gene. Your method should be generalisable to a list of genes and take a few seconds to run (take inspiration from question 5). Then, please report the following:
  
 -Find out the coordinates of the _ADIPOQ_ gene.
    ```
    wget -q --header='Content-type:application/json' 'http://rest.ensembl.org/overlap/id/ENSG00000181092?feature=gene'  -O -
    This gene maps to Chromosome 3: 186,842,704-186,858,463 forward strand.
    ```
 
 - the coordinates of the exons of its canonical transcript.
  
  ```
[{"db_type":"core","end":186560569,"start":186560463,"assembly_name":"GRCh37","species":"homo_sapiens","version":1,"id":"ENSE00002239164","strand":1,"object_type":"Exon","seq_region_name":"3"},

{"assembly_name":"GRCh37","start":186570814,"end":186571061,"db_type":"core","seq_region_name":"3","id":"ENSE00001799162","strand":1,"object_type":"Exon","version":2,"species":"homo_sapiens"},

{"assembly_name":"GRCh37","start":186571973,"end":186576252,"db_type":"core","seq_region_name":"3","object_type":"Exon","strand":1,"id":"ENSE00001765971","version":1,"species":"homo_sapiens"}],"biotype":"protein_coding","end":186576252,"is_canonical":1,"source":"ensembl"}]

 ```

- all documented variants in this gene.

```
wget -q --header='Content-type:application/json' 'http://rest.ensembl.org/overlap/id/ENSG00000181092?feature=variation'  -O - >> variants.txt

```  
  - all phenotype-associated variants. 
  
  - all documented loss-of-function (LoF) variants in this gene. How did you define LoF?
  - If asked to find all regulatory variants that potentially affect this gene, which questions would you ask and how would you proceed?
    

12. How would you convert a VCF file to the Plink binary format? How would you do the reverse, and what kind of problems do you anticipate?

``` 
plink --vcf myvcf.vcf.gz --make-bed --out myfileout

```

13. Write a snippet to reformat a PED file so as to output a file with the following header `sample_name genotype_SNP1 genotype_SNP2 ...` where genotypes are coded VCF-style (e.g `A/C`, the order of the alleles in the output is not important).
``` 
```

14. A genetic association pipeline starts with a VCF and produces summary statistics for every position. The VCF contains multiallelics and indels. Unfortunately, a program in the pipeline trims all alleles to their first character. Why might allele frequencies not always be equal for a given variant? Find a way to correct the alleles in the association file by using the information from the VCF. Select columns are provided for [the association file](https://github.com/hmgu-itg/challenge/raw/master/data/association.txt.gz). We also provide [a file](https://github.com/hmgu-itg/challenge/raw/master/data/fromVCF.txt.gz) that was created from the VCF using `bcftools query -f '%CHROM %POS %REF %ALT %AN %AC\n'`.
15. [This file](https://github.com/hmgu-itg/challenge/raw/master/data/mock.eQTL.forChallenge.txt) contains eQTL overlap data for SNPs that arise as signals in GWAS for several phenotypes. Reformat this file to have one line per SNP/phenotype pair, and two additional columns formatted as such : `GENE1(tissue1, tissue2),GENE2(tissue1, tissue3)`, and `GENE1(2),GENE2(2)`. Each line should contain the SNP/phenotype pair, all genes found overlapping and their respective tissues, and all genes found overlapping with the number of tissues.
16. A researcher wants to conduct a disease association study. However, colleagues warn him that the dataset contains related individuals. He would like to remove relatedness in his dataset, but given his disease is rare, he would also like to maximise the number of cases kept in. Using [a list of samples with disease status](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.pheno.tsv) and [a file containing pairs of individuals above a relatedness threshold](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.tsv), create an exclusion list of samples to remove to help the researcher achieve their goal.


### Statistical genetics
1. You sample at random 10,000 variants from a deep (50x) whole-genome sequencing variant call file describing 1,000 individuals. What do you expect the distribution of minor allele frequency to look like? In particular, which minor allele counts are likely to be most frequent?

```
Since the sequencing depth is relatively high, I would expect to see an exponential distribution, with the rarest variants (MAF=0.1%, i.e. unique to each individual) being the most common and the majority of variants falling within a MAF value below 5%.
```

2. You are running a single-point association study with a quantitative phenotype on the dataset above. Which filters, if any, would you apply to the dataset? 
```
Since the QTLs need to be robustly associated (with high confidence i.e. low p-value), their effect must be present in a considerable number of individuals, therefore I would disregard SNPs with MAF<5%. Furthermore, I would run all the standard QC filters i.e. remove any SNPs with Hardy Weinberg Equilibrium Chi^2 p-value <10e-6; remove any SNPs with a call rate < 0.98; remove SNPs that deviate strongly (> 3 standard deviations) from the sample's heterozygosity rate mean.
```

3. A common practice when performing genetic association studies is to perform an ethnicity check as a quality control step. Can you explain how this is done?
    - You are dealing with whole-genome sequencing data in 2,326 Bulgarian samples. How would you perform such an ethnicity check, and which projection dataset would you use? 
```
The ethnicity check is typically done by running a Principle Component Analysis of the study samples together with individuals from a validated cohort (for example the 1000 Genomes Project Cohort) to obtain a projection of each individual onto the plane described by the first two principal components, which are orthogonal vectors along which all samples exhibit the highest proportion variance. All Bulgarian individuals are expected to cluster together, within the European clusters; any individuals falling outside of the main cluster should be removed, particularly if they are detected to be of non-European ethnicity.
```

4. You perform a single-point association with blood lipids and find a variant with MAF=0.7% associated at p=1e-14. Do you expect the effect size to be large or small? What would be your next steps investigating this signal?
```
Since MAF for this variant is very small and so is the p-value I would expect the effect size to be very large. I would try to make sure that the result makes biological sense so I would use Ensembl's VEP tool to identify possible molecular consequences of the variant (does it affect a protein coding gene, or a regulatory region of a gene? what is the function of the nearby gene, does it make sense that it could affect blood lipid levels?)
```

5. You are running an inverse-variance based single-point meta-analysis of the above dataset together with a UK study of 12,400 individuals. The variant above is found in the UK dataset, but its association is weak (1e-3). However, another variant located 1kb downstream is strongly associated with the same trait (p=1e-15) but weakly associated in your study (p=5e-4). Which of these 2 variants do you expect will have the strongest meta-analysis p-value? What do you think is happening in this region, how can you test it, and which model could you apply if it is the case?

```
I think this may be a situation where the two variants are in linkage disequilibrium so their alleles are non-randomly associated (one could say they are correlated). I would expect that the variant from the British cohort will have a stronger p-value, primarily due to the cohort's size.
```
6. An analyst studies a population of remote villages in Eastern Europe. They are interested in a particular variant, and compare the frequency in their villages (3.5%) to the EUR population frequency in the 1000 Genomes (0.03%). They conclude that the variant has increased in frequency in their villages. Do you agree, and if not, what would your advice be?

```
The catch here is that the 1000 Genomes Projects includes individuals from all across the globe. The 'control' MAF frequency should obtained only from eastern European samples, and then the variant's enrichment in the villages can be assessed.
```

7.  The same analyst sends you association summary statistics for random glucose.
    - Which checks would you perform on such a dataset?
    - You wish to include this dataset in a meta-analysis. Which additional information should you ask for in your next email to your colleague?
    - In this dataset, you observe  &#955;=1.25. The analyst has adjusted for age, age-squared, sex, and fasting status. What would you suggest they do?
    
```
Compare against results obtained for another cohort which had already been published. Maybe correct for additional confounding factors such as diet/lifestyle?
```

8. You are a co-author on a manuscript. You receive a draft, in which the main author has used the traditional &#945;=5e-8 as the significance threshold. The paper describes an analysis of 10 related blood phenotypes (platelet count, platelet volume, immature platelet fraction ...) using the fixed markers of the Infinium ImmunoArray on 897 individuals. What do you think about the chosen threshold, and what would you suggest to the first author? What would be your comments on the design of the study? 

```
Since multiple related blood phenotypes are tested for, multiple hypothesis testing occurs so I would suggest using a more stringent significant threshold reduced by a factor of 10 (Bonferroni correction). 
```
