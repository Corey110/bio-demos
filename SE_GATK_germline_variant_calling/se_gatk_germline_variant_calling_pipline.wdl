version 1.0
# 动态interval版
import './se_gatk_germline_variant_calling_tasks.wdl' as tasks
workflow wgs {

    input {
        String sample_id
    
        File fastq1

        Boolean phred64 
        Boolean fix_mgi_id
        String? adapter_sequence
        String? adapter_sequence_r2
        Int? reads_to_process 

        Array[File] genome_indexes
        
        Int bwa_cpu
        Int bwa_memory

        String docker_img_germline_tools = 'cr-cn-beijing.volces.com/popgenomics/germline_tools:v3'
        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'


        Int scatter_count
    }

    call tasks.fastp_clean {
        input:
            sample_id = sample_id,
            fastq1 = fastq1,
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
            fastq1 = fastp_clean.cleaned_fq,
            genome_indexes = genome_indexes,
            cpu = bwa_cpu,
            memory = bwa_memory,
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