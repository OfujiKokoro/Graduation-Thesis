merge vcfはgzで圧縮し、indexをつける

bcftools merge --threads 10 -l Sample_List_row_131_vcfgz.txt > bed_vcf_131_GT105.vcf
