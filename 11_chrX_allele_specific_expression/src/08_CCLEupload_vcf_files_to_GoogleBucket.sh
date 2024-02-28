use Google-Cloud-SDK

File='./11_chrX_allele_specific_expression/output/06_CCLE_prepare_for_ASEReadCounter/gSNPs_PASS_VCFs_toBeUploadedToGoogleBucket.txt'
Lines=$(cat $File)

for Line in $Lines
do
  echo $Line
  gsutil -u gdan-copy-number-gdac cp $Line gs://fc-secure-516f8144-b7c9-458a-934f-2cc78f09cfb3
done
