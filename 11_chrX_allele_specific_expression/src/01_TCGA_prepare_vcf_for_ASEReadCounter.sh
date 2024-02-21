cd /xchip/beroukhimlab/kei/project/TangentXY_script/11_chrX_allele_specific_expression

vcf_file='data/PanCan-Germline/PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz'
out_dir='output/01_TCGA_prepare_vcf_for_ASEReadCounter'

## Extract sample list from the vcf file
bcftools query -l $vcf_file > $out_dir/'sample.list.txt'

## Subset biallelic SNPs on the master VCF file
gatk \
  SelectVariants \
  -R /xchip/beroukhimlab/kei/resource/Homo_sapiens_assembly19.fasta \
  -V $vcf_file \
  -O $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.vcf.gz" \
  -L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y \
  --restrict-alleles-to BIALLELIC \
  --select-type-to-include SNP

## Convert VCF to BCF
bcftools \
  view \
  -Ou $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.vcf.gz" \
  -o $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.bcf"
bcftools index $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.bcf"

## Subset BCF by chromosomes
for chr in {1..22} X Y; do
  bcftools view -Ou -r ${chr} -o $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs_chr"${chr}".bcf" $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.bcf" 
  bcftools index $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs_chr"${chr}".bcf"
done

## Subset chr1-22 to construct bcf for autosomes
bcftools view -Ou -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -o $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs_autosomes.bcf" $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.bcf"
bcftools index $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs_autosomes.bcf"

## Extract one sample from the master VCF file
sample_name="CESC.FU.A3WB"
sample_id="TCGA-FU-A3WB-10A-01D-A22X-09"
vcf=$out_dir"/TCGA-FU-A3WB-10A-01D-A22X-09.vcf.gz"
bam="/xchip/beroukhimlab/kei/project/Tangent/20230130_ASEReadCounter/data/bam/9f90820c-7a48-42b0-9ad0-1db04eb273a0/UNCID_1724291.A1DBD227-73B7-466F-B766-063A3006A206.sorted_genome_alignments.bam"

bcftools view -Oz -s $sample_id -o $out_dir/$sample_id"_Biallelic_SNPs.vcf.gz" $out_dir/"PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389_Biallelic_SNPs.bcf"
bcftools index -t $out_dir/$sample_id"_Biallelic_SNPs.vcf.gz"
