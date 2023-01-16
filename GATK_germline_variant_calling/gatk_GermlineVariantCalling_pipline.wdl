version 1.0
# 动态interval版
import './gatk_germline_variant_calling_tasks.wdl' as tasks
workflow wgs {

    input {
        String sample_id
    
        File fastq1
        File fastq2

        Boolean phred64 
        Boolean fix_mgi_id
        String? adapter_sequence
        String? adapter_sequence_r2
        Int? reads_to_process 

        Array[File] genome_indexes

        String docker_img_germline_tools # 修改自己的镜像仓库地址
        String docker_img_gatk4 # 修改自己的镜像仓库地址

        Int scatter_count
    }

    call tasks.fastp_clean {
        input:
            sample_id = sample_id,
            fastq1 = fastq1,
            fastq2 = fastq2,
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
    
    call tasks.mosdepth {
        input:
            sample_id = sample_id,
            aligned_bam = bwa_align.sorted_bam_output,
            docker_img = docker_img_germline_tools
    }

    call tasks.dedup {
        input:
            sample_id = sample_id,
            sorted_bam = bwa_align.sorted_bam_output,
            docker_img = docker_img_gatk4
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
        Pair[File,File] fastp_report = fastp_clean.fastp_report
        Array[File] mosdepth_res = mosdepth.mosdepth_res
        Pair[File,File] output_bam = dedup.dedup_bam_output
        File dedup_metrics = dedup.dedup_metrics_output
        Pair[File,File] output_gvcf = merge_vcfs.gvcf_output
        Pair[File,File] output_vcf = genotype_gvcfs.vcf_output
        File bam_flagstat = mosdepth.flagstate_res
    }

}
