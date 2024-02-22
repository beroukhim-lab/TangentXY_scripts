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
  	Float mem_size_gb = 20
  	Int preemptible_tries = 3
  	Int compression_level = 5
  	Int flowcell_medium_disk = 200

	call DepthOfCoverage {
		input: 
			input_bam = input_bam,
			input_bam_index = input_bam_index,
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_fasta_dict = ref_fasta_dict,
			interval_list = interval_list,
			sample_id = sample_id,
			docker_image = gatk_docker,
			preemptible_tries = preemptible_tries,
			compression_level = compression_level,
			disk_size = flowcell_medium_disk
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
	Float mem_size_gb
	Int command_mem_gb
  	Int compression_level
  	Int preemptible_tries
  	Int disk_size

	command {
		${gatk_path} --java-options "-Dsamjdk.compression_level=${compression_level} -Xms${command_mem_gb}G" \
			DepthOfCoverage \
			-R ${ref_fasta} \
			-I ${input_bam} \
			-O ${sample_id} \
			-L ${interval_list} \
	}

	output {
		File output_sample_cumulative_coverage_counts = "${sample_id}.sample_cumulative_coverage_counts"
		File output_sample_cumulative_coverage_proportions = "${sample_id}.sample_cumulative_coverage_proportions"
		File output_sample_interval_summary = "${sample_id}.sample_interval_summary"
		File output_sample_interval_statistics = "${sample_id}.sample_interval_statistics"
		File output_sample_statistics = "${sample_id}.sample_statistics"
		File output_sample_summary = "${sample_id}.sample_summary"
	}

  runtime {
  	preemptible: preemptible_tries
    docker: docker_image
    memory: "${mem_size_gb} GB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
	}
}
