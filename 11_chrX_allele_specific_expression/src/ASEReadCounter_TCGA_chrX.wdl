workflow ASEReadCounterWorkflow {
    String vcf_sample_id
    String sample_id
	File input_bcf
    File chain_file
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
    
    call SubsetSample {
    	input:
        	vcf_sample_id = vcf_sample_id,
            input_bcf = input_bcf,
            preemptible = preemptible,
            mem_size_gb = mem_size_gb,
            cpu_num = cpu_num,
            disk_size = disk_size
    }
    call FilterHeterozygousSNPs {
    	input:
        	vcf_sample_id = vcf_sample_id,
        	input_vcf = SubsetSample.subset_vcf,
            preemptible = preemptible,
            mem_size_gb = mem_size_gb,
            cpu_num = cpu_num,
            disk_size = disk_size
    }
    call LiftOverVCF {
    	input:
        	vcf_sample_id = vcf_sample_id,
        	input_vcf = FilterHeterozygousSNPs.hetero_vcf,
            ref_fasta = ref_fasta,
            ref_fasta_dict = ref_fasta_dict,
            ref_fasta_index = ref_fasta_index,
            chain_file = chain_file,
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
        	input_vcf = LiftOverVCF.lifted_vcf,
        	input_vcf_index = LiftOverVCF.lifted_vcf_index,
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

task SubsetSample {
    String vcf_sample_id
	File input_bcf
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size

	command {
    	bcftools view -Oz -s ${vcf_sample_id} -o ${vcf_sample_id}_Biallelic_SNPs.vcf.gz ${input_bcf}
        #bcftools index -t ${vcf_sample_id}_Biallelic_SNPs.vcf.gz
    }
    
    output {
    	File subset_vcf = "${vcf_sample_id}_Biallelic_SNPs.vcf.gz"
    	#File subset_vcf_index = "${vcf_sample_id}_Biallelic_SNPs.vcf.gz.tbi"
    }

	runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}

task FilterHeterozygousSNPs {
    String vcf_sample_id
	File input_vcf
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size
    
	command {
    	bcftools view -Oz -o ${vcf_sample_id}_Biallelic_SNPs_Hetero.vcf.gz -i 'GT="0/1"' ${input_vcf}
    }
    
    output {
    	File hetero_vcf = "${vcf_sample_id}_Biallelic_SNPs_Hetero.vcf.gz"
    }
    
	runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}

task LiftOverVCF {    
    String vcf_sample_id
	File input_vcf
    File chain_file
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_index
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size

	command {
    	gatk LiftoverVcf \
        	-I ${input_vcf} \
            -O ${vcf_sample_id}_Biallelic_SNPs_Hetero_LO.vcf.gz \
            -CHAIN ${chain_file} \
            -REJECT ${vcf_sample_id}_rejected_SNPs.vcf.gz \
            -R ${ref_fasta}
    }
    
    output {
    	File lifted_vcf = "${vcf_sample_id}_Biallelic_SNPs_Hetero_LO.vcf.gz"
    	File lifted_vcf_index = "${vcf_sample_id}_Biallelic_SNPs_Hetero_LO.vcf.gz.tbi"
    }
    
	runtime {
        docker: "broadinstitute/gatk:4.3.0.0"
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
            -O ${sample_id}_ASE_chrX.txt \
            --read-index ${input_bam_index} \
            --min-mapping-quality ${min_mapping_quality} \
            --min-base-quality ${min_base_quality} \
            --min-depth-of-non-filtered-base ${min_depth_of_filtered_base} \
            -L chrX
    }
    
    output {
    	File ase_output = "${sample_id}_ASE_chrX.txt"
    }
    
	runtime {
        docker: "broadinstitute/gatk:4.3.0.0"
    	preemptible: "${preemptible}"
        memory: "${mem_size_gb} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
	}
}
