version 1.0
# 动态interval版
import './gatk_germline_variant_calling_tasks.wdl' as tasks
workflow wgs {

    input {
        String sample_id

        Array[File] input_bam

        Boolean if_down_bam
        Int? target_depth
        Float? origin_depth

        Array[File] genome_indexes

        String docker_img_germline_tools = 'cr-cn-beijing.volces.com/popgenomics/germline_tools:v3'
        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
        String docker_img_mosdepth = 'cr-cn-beijing.volces.com/popgenomics/mosdepth:latest'

        Int scatter_count
    }

    Pair[File, File] input_bams = (input_bam[0], input_bam[1])

    if (if_down_bam) {
        call tasks.down_sample_bam {
            input:
                bam_file = input_bams,
                sample_id = sample_id,
                origin_depth = select_first([origin_depth]),
                target_depth = select_first([target_depth]),
                docker_img = docker_img_germline_tools
        }

        call tasks.mosdepth {
            input:
                sample_id = sample_id,
                aligned_bam = down_sample_bam.down_bam,
                docker_img = docker_img_mosdepth
        }
    }

    call tasks.split_intervals {
        input: 
            genome_indexes = genome_indexes,
            scatter_count = scatter_count,
            docker_img = docker_img_gatk4
    }

    Pair[File, File] call_snp_bam = select_first([down_sample_bam.down_bam, input_bams])

    scatter (intervals in split_intervals.scatter_interval_list_file) {
        call tasks.haplotype_caller {
            input:
                sample_id = sample_id,
                input_bam = call_snp_bam,
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
        Array[File]? mosdepth_res = mosdepth.mosdepth_res
        File? bam_flagstat = mosdepth.flagstate_res
        Pair[File, File]? output_bam = down_sample_bam.down_bam
        Array[File] output_gvcf = [merge_vcfs.gvcf_output.left, merge_vcfs.gvcf_output.right]
        Array[File] output_vcf = [genotype_gvcfs.vcf_output.left, genotype_gvcfs.vcf_output.right]
        Array[String]? stat_summary = mosdepth.stat_summary
        String? depth = mosdepth.depth
        String? mapped_rate = mosdepth.mapped_rate
    }

}
