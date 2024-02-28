workflow ASEReadCounterWorkflow {
    String sample_id
	File input_vcf
	File input_vcf_index
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_fasta_dict
    
    Int min_mapping_quality = 20
    Int min_base_quality = 30
    Int min_depth_of_filtered_base = 20

    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size
    
    call FilterBiallelicSNPs {
    	input:
        	sample_id = sample_id,
        	input_vcf = input_vcf,
        	input_vcf_index = input_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_index = ref_fasta_index,
            preemptible = preemptible,
            mem_size_gb = mem_size_gb,
            cpu_num = cpu_num,
            disk_size = disk_size
    }
    call FilterHeterozygousSNPs {
    	input:
        	sample_id = sample_id,
        	input_vcf = FilterBiallelicSNPs.biallelic_vcf,
            preemptible = preemptible,
            mem_size_gb = mem_size_gb,
            cpu_num = cpu_num,
            disk_size = disk_size
    }
    call ASEReadCounter {
    	input:
        	sample_id = sample_id,
            input_bam = input_bam,
            input_bam_index = input_bam_index,
        	input_vcf = FilterHeterozygousSNPs.hetero_vcf,
        	input_vcf_index = FilterHeterozygousSNPs.hetero_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_index = ref_fasta_index,
            min_mapping_quality = min_mapping_quality,
            min_base_quality = min_base_quality,
            min_depth_of_filtered_base = min_depth_of_filtered_base,
            preemptible = preemptible,
            mem_size_gb = mem_size_gb,
            cpu_num = cpu_num,
            disk_size = disk_size
    }
}
task FilterBiallelicSNPs {
	String sample_id
    File input_vcf
    File input_vcf_index
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size
    
    command {
    	gatk SelectVariants \
        	-R ${ref_fasta} \
            -V ${input_vcf} \
            -O ${sample_id}_Biallelic_SNPs.vcf.gz \
            --restrict-alleles-to BIALLELIC \
            --select-type-to-include SNP
    }
    
    output {
    	File biallelic_vcf = "${sample_id}_Biallelic_SNPs.vcf.gz"
    }
	runtime {
        docker: "broadinstitute/gatk:4.3.0.0"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}

task FilterHeterozygousSNPs {
    String sample_id
	File input_vcf
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size
    
	command {
    	bcftools view -Oz -o ${sample_id}_Biallelic_SNPs_Hetero.vcf.gz -i 'GT="0/1"' ${input_vcf}
        bcftools index -t ${sample_id}_Biallelic_SNPs_Hetero.vcf.gz
    }
    
    output {
    	File hetero_vcf = "${sample_id}_Biallelic_SNPs_Hetero.vcf.gz"
    	File hetero_vcf_index = "${sample_id}_Biallelic_SNPs_Hetero.vcf.gz.tbi"
    }
    
	runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}

task ASEReadCounter {
    String sample_id
	File input_bam
	File input_bam_index
	File input_vcf
	File input_vcf_index
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index
    Int min_mapping_quality
    Int min_base_quality
    Int min_depth_of_filtered_base
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size
    
	command {
    	gatk ASEReadCounter \
        	-R ${ref_fasta} \
            -I ${input_bam} \
            -V ${input_vcf} \
            -O ${sample_id}_ASE.txt \
            --read-index ${input_bam_index} \
            --min-mapping-quality ${min_mapping_quality} \
            --min-base-quality ${min_base_quality} \
            --min-depth-of-non-filtered-base ${min_depth_of_filtered_base} \
            -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY
    }
    
    output {
    	File ase_output = "${sample_id}_ASE.txt"
    }
    
	runtime {
        docker: "broadinstitute/gatk:4.3.0.0"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}
