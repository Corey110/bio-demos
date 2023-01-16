version 1.0
# 动态interval版
import './gatk_germline_variant_calling_tasks.wdl' as tasks
workflow wgs {

    input {
        String sample_id
    
        File fastq1
        File fastq2

        Boolean if_down_fastq
        Int? target_depth
        Float? origin_depth

        Boolean phred64 
        Boolean fix_mgi_id
        String? adapter_sequence
        String? adapter_sequence_r2
        Int? reads_to_process 

        Array[File] genome_indexes

        String docker_img_germline_tools = 'cr-cn-beijing.volces.com/popgenomics/germline_tools:v3'
        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
        String docker_img_seqtk = 'cr-cn-beijing.volces.com/popgenomics/seqtk:1.3'
        String docker_img_mosdepth = 'cr-cn-beijing.volces.com/popgenomics/mosdepth:latest'

        Int scatter_count
    }

    if (if_down_fastq){
        call tasks.down_sample_fastq {
            input:
                sample_id = sample_id,
                fastq1 = fastq1,
                fastq2 = fastq2,
                origin_depth = select_first([origin_depth]),
                target_depth = select_first([target_depth]),
                docker_img = docker_img_seqtk
        }
    }

    File run_fastq1 = select_first([down_sample_fastq.output_fastq1, fastq1])
    File run_fastq2 = select_first([down_sample_fastq.output_fastq2, fastq2])

    call tasks.fastp_clean {
        input:
            sample_id = sample_id,
            fastq1 = run_fastq1,
            fastq2 = run_fastq2,
            phred64 = phred64,
            fix_mgi_id = fix_mgi_id,
            adapter_sequence = adapter_sequence,
            adapter_sequence_r2 = adapter_sequence_r2,
            reads_to_process = reads_to_process,
            docker_img = docker_img_germline_tools
    }

    call tasks.bwa_align {
        input: 
            sample_id = sample_id,
            fastq1 = fastp_clean.cleaned_fq.left,
            fastq2 = fastp_clean.cleaned_fq.right,
            genome_indexes = genome_indexes,
            docker_img = docker_img_germline_tools
    }
    
    call tasks.dedup {
        input:
            sample_id = sample_id,
            sorted_bam = bwa_align.sorted_bam_output,
            docker_img = docker_img_gatk4
    }

    call tasks.mosdepth {
        input:
            sample_id = sample_id,
            aligned_bam = dedup.dedup_bam_output,
            docker_img = docker_img_mosdepth
    }

    call tasks.split_intervals {
        input: 
            genome_indexes = genome_indexes,
            scatter_count = scatter_count,
            docker_img = docker_img_gatk4
    }

    scatter (intervals in split_intervals.scatter_interval_list_file) {
        call tasks.haplotype_caller {
            input:
                sample_id = sample_id,
                input_bam = dedup.dedup_bam_output,
                genome_indexes = genome_indexes,
                calling_region = intervals,
                docker_img = docker_img_gatk4
        }
    }

    call tasks.merge_vcfs {
        input:
            sample_id = sample_id,
            vcfs = haplotype_caller.hc_block_gvcf,
            docker_img = docker_img_gatk4
            
    }

    call tasks.genotype_gvcfs {
        input: 
            sample_id = sample_id,
            genome_indexes = genome_indexes,
            gvcf = merge_vcfs.gvcf_output,
            docker_img = docker_img_gatk4
    }
    
    output {
        Array[File] fastp_report = [fastp_clean.fastp_report.left, fastp_clean.fastp_report.right]
        Array[File] mosdepth_res = mosdepth.mosdepth_res
        Array[File] output_bam = [dedup.dedup_bam_output.left, dedup.dedup_bam_output.right]
        File dedup_metrics = dedup.dedup_metrics_output
        Array[File] output_gvcf = [merge_vcfs.gvcf_output.left, merge_vcfs.gvcf_output.right]
        Array[File] output_vcf = [genotype_gvcfs.vcf_output.left, genotype_gvcfs.vcf_output.right]
        File bam_flagstat = mosdepth.flagstate_res
        Array[String] stat_summary = mosdepth.stat_summary
        String depth = mosdepth.depth
        String mapped_rate = mosdepth.mapped_rate
    }

}
