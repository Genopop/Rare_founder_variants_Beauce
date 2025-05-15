module load StdEnv/2020 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16
module load java/11.0.2
module load StdEnv/2023 java/21.0.1
export JAVA_TOOL_OPTIONS=-Xmx20g

for chr in {1..22}; do
  ## Phase with Beagle
  java -Xmx29127m -jar ${soft_dir}beagle.18May20.d20.jar gt=${vcf_in}${chr}.vcf out=${vcf_in}${chr}
  ## Split into clusters
  vcftools --gzvcf ${vcf_in}${chr}.vcf.gz --keep ${cluster1} --recode --recode-INFO-all --out ${outfile}${chr}_beauce_cluster
  bgzip -f ${outfile}${chr}_beauce_cluster.recode.vcf
  vcftools --gzvcf ${vcf_in}${chr}.vcf.gz --keep ${cluster2} --recode --recode-INFO-all --out ${outfile}${chr}_urbanqc_cluster
  bgzip ${outfile}${chr}_urbanqc_cluster.recode.vcf
  vcftools --gzvcf ${vcf_in}${chr}.vcf.gz --keep ${cluster3} --recode --recode-INFO-all --out ${outfile}${chr}_slsj_cluster
  bgzip ${outfile}${chr}_slsj_cluster.recode.vcf
  ## run Hap-IBD #hap-ibd.jar  [ version 1.0, 15Jun23.92f ]
  java -jar hap-ibd.jar gt=${outfile}${chr}_beauce_cluster.recode.vcf.gz map=${genmap}${chr}.GRCh38.map out=${outfile}${chr}_beauce_cluster_hap-ibd.out
  java -jar hap-ibd.jar gt=${outfile}${chr}_urbanqc_cluster.recode.vcf.gz map=${genmap}${chr}.GRCh38.map out=${outfile}${chr}_urbanqc_cluster_hap-ibd.out
  java -jar hap-ibd.jar gt=${outfile}${chr}_slsj_cluster.recode.vcf.gz map=${genmap}${chr}.GRCh38.map out=${outfile}${chr}_slsj_cluster_hap-ibd.out
done

## Combine genetic maps
cat ${genmap}{1..22}.GRCh38.map > ${combined_genmap}
## Compute Ne
gunzip ${outfile}{1..22}_beauce_cluster_hap-ibd.out.ibd.gz
gunzip ${outfile}{1..22}_urbanqc_cluster_hap-ibd.out.ibd.gz
gunzip ${outfile}{1..22}_slsj_cluster_hap-ibd.out.ibd.gz
cat ${outfile}{1..22}_beauce_cluster_hap-ibd.out.ibd | java -jar ${soft_dir}ibdne.23Apr20.ae9.jar map=${combined_genmap} out=${outfile}_beauce_cluster_hap-ibd_ne nthreads=8
cat ${outfile}{1..22}_urbanqc_cluster_hap-ibd.out.ibd | java -jar ${soft_dir}ibdne.23Apr20.ae9.jar map=${combined_genmap} out=${outfile}_urbanqc_cluster_hap-ibd_ne nthreads=8
cat ${outfile}{1..22}_slsj_cluster_hap-ibd.out.ibd | java -jar ${soft_dir}ibdne.23Apr20.ae9.jar map=${combined_genmap} out=${outfile}_slsj_cluster_hap-ibd_ne nthreads=8


  
