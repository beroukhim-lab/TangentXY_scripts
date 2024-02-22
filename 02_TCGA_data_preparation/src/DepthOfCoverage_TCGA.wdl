## Modified from DepthOfCoverage.5.wdl that was created by Coyin
## https://portal.firecloud.org/?return=terra#methods/coyinoh/DepthOfCoverage/5/wdl
workflow DepthOfCoverageWorkflow {
	File input_bam
	File input_bam_index
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File interval_list
	String sample_id

  	String gatk_docker = "broadinstitute/gatk:4.2.6.1"
  	String gatk_path = "/gatk/gatk"
    Int preempitble = 1
  	Int mem_size_gb = 20
    Int cpu_num = 2
    Int disk_size = 50

	call DepthOfCoverage {
		input: 
			input_bam = input_bam,
			input_bam_index = input_bam_index,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_fasta_dict = ref_fasta_dict,
			interval_list = interval_list,
			sample_id = sample_id,
			docker_image = gatk_docker
  	}
}

task DepthOfCoverage {
	File input_bam
	File input_bam_index
	File ref_fasta
	File ref_fasta_index
	File ref_fasta_dict
	File interval_list
	String sample_id

	String gatk_docker
	String gatk_path
    String docker_image
    Int preemptible
	Int mem_size_gb
    Int cpu_num
  	Int disk_size

	command {
		${gatk_path} \
			DepthOfCoverage \
			-R ${ref_fasta} \
			-I ${input_bam} \
			-O ${sample_id} \
			-L ${interval_list} \
            --read-index ${input_bam_index} \
            --omit-depth-output-at-each-base true \
            --omit-locus-table true \
            --omit-per-sample-statistics true
	}

	output {
		#File output_sample_cumulative_coverage_counts = "${sample_id}.sample_cumulative_coverage_counts"
		#File output_sample_cumulative_coverage_proportions = "${sample_id}.sample_cumulative_coverage_proportions"
		File output_sample_interval_summary = "${sample_id}.sample_interval_summary"
		File output_sample_interval_statistics = "${sample_id}.sample_interval_statistics"
		#File output_sample_statistics = "${sample_id}.sample_statistics"
		#File output_sample_summary = "${sample_id}.sample_summary"
	}

    runtime {
      preemptible: preemptible
      docker: docker_image
      memory: "${mem_size_gb} GB"
      cpu: "${cpu_num}"
      disks: "local-disk " + disk_size + " HDD"
      }
}
