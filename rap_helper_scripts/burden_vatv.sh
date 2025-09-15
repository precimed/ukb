
# Specify phenotype file and BGEN genotype file prefix
phenotype_file1="/zillur/bodymri/pheno/VATV.csv"
covariate_file1="/zillur/bodymri/pheno/covariates_age_sex_bmi_pc10_caucasian.csv"
genotype_prefix="/ukb_c1-22_genome.filtered3_merged"
bgen_prefix="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release"
#ldl_custom="project-G64PFY0JyZk3Gq1xBQYfk45Q:/norment/metabolimics/ldl_custom_masks.txt"

data_field_bgen="ukb23159"
data_out1="/zillur/bodymri/gwas/vatv/"

run_plink="plink2 --bfile ukb_c1-22_genome.filtered3_merged \
 --no-pheno --keep VATV.csv --autosome \
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
 --write-snplist --write-samples --no-id-header \
 --out qc_vatv_pass --threads 32 --memory 16000"

#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_file1}" \
# -icmd="${run_plink}" --tag="qc_vatv" --instance-type "mem2_ssd1_v2_x32" \
# --destination="${project}:${data_out1}" --brief --yes

all_sample_snplist="/zillur/bodymri/gwas/vatv/qc_vatv_pass.snplist"
all_sample_id="/zillur/bodymri/gwas/vatv/qc_vatv_pass.id"

regenie_step1="regenie --step 1 \
 --lowmem --out step1_vatv --bed ukb_c1-22_genome.filtered3_merged \
 --phenoFile VATV.csv \
 --covarFile covariates_age_sex_bmi_pc10_caucasian.csv \
 --extract qc_vatv_pass.snplist \
 --bsize 1000 --loocv --gz --threads 96"

#dx run swiss-army-knife -iin="${genotype_prefix}.bed" \
# -iin="${genotype_prefix}.bim" \
# -iin="${genotype_prefix}.fam" \
# -iin="${phenotype_file1}"  \
# -iin="${covariate_file1}"  \
# -iin="${all_sample_snplist}" \
# -icmd="${regenie_step1}" \
# --tag="burden_step1_vatv" --instance-type "mem1_hdd1_v2_x96" \
# --destination="${project}:${data_out1}" --brief --yes

pred_folder="/zillur/bodymri/gwas/vatv"
path_to_500kwes_helper_files="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/helper_files"
ldl_custom="/zillur/norment/metabolimics/ldl_custom_masks.txt"
ready_snp_list="/zillur/norment/metabolimics/gene_sets_all.txt"

pred1="/zillur/bodymri/gwas/vatv/step1_vatv_pred.list"

#dx run swiss-army-knife -icmd="cp /mnt/project/zillur/bodymri/gwas/vatv/*loco.gz . ; zip loco.zip *loco.gz ; rm *loco.gz" --destination /zillur/bodymri/gwas/vatv/

for i in {1..22}; do
    regenie_burden="regenie \
    --step 2 --pred step1_vatv_pred.list \
    --bgen ukb23159_c${i}_b0_v1.bgen \
    --ref-first \
    --sample ukb23159_c${i}_b0_v1.sample \
    --phenoFile VATV.csv \
    --covarFile covariates_age_sex_bmi_pc10_caucasian.csv \
    --set-list ukb23158_500k_OQFE.sets.txt.gz \
    --anno-file ukb23158_500k_OQFE.annotations.txt.gz \
    --mask-def ldl_custom_masks.txt \
    --aaf-bins 0.005 --nauto 23 \
    --bsize 200 --extract-sets gene_sets_all.txt \
    --write-mask-snplist --write-mask \
    --vc-tests skato,acato-full \
    --out vatv_burden_chr${i}"
  
 dx run swiss-army-knife -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen" \
     -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.sample" \
     -iin="${bgen_prefix}/${data_field_bgen}_c${i}_b0_v1.bgen.bgi" \
     -iin="${ldl_custom}" \
     -iin="${ready_snp_list}" \
     -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.sets.txt.gz" \
     -iin="${path_to_500kwes_helper_files}/ukb23158_500k_OQFE.annotations.txt.gz" \
     -iin="${phenotype_file1}" \
     -iin="${covariate_file1}" \
     -iin="${pred1}" \
     -iin="${pred_folder}/loco.zip" \
     -icmd="unzip loco.zip ;  ${regenie_burden}" \
     --tag="burden_step2_vatv" --instance-type "mem2_ssd1_v2_x96" \
     --destination="${project}:${data_out1}" --brief --yes
done
