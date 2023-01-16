version 1.0
 
import "./filter_variant_tasks.wdl" as tasks

workflow filter_variant {
    input {
        # 直接输入多样本VCF需要将其index文件一起作为Array输入
        Array[File] vcf_to_filter
        
        Array[File] genome_indexes

        String output_vcf_basename

        Int disk_size_gb_snp = 50
        Int disk_size_gb_indel = 50
        Int disk_size_gb_filter_snp = 50
        Int disk_size_gb_filter_indel = 50

        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
        String docker_img_bcftools = 'cr-cn-beijing.volces.com/popgenomics/bcftools:1.9'
        String docker_img_plink = 'cr-cn-beijing.volces.com/popgenomics/plink:1.9'

        Boolean if_select_snp
        Boolean if_select_indel

        Boolean if_filter_snp
        Boolean if_filter_indel

        Int scatter_count

    }

    call tasks.bcftools_stats {
        input:
            input_vcfs = vcf_to_filter,
            output_basename = output_vcf_basename,
            docker_img = docker_img_bcftools,
            sufix = 'total'
    }

    call tasks.split_intervals {
            input: 
                genome_indexes = genome_indexes,
                scatter_count = scatter_count,
                docker_img = docker_img_gatk4
    }

    if (if_select_snp) {
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.select_snp {
                input:
                    input_vcfs = vcf_to_filter,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals,
                    disk_size_gb = disk_size_gb_snp
            }        
        }

        call tasks.merge_vcfs as merge_snp_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = select_snp.vcf_snp,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'snp'

        }

        call tasks.bcftools_stats as select_snp_stats{
            input:
                input_vcfs = merge_snp_vcfs.vcf_output,
                output_basename = output_vcf_basename,
                docker_img = docker_img_bcftools,
                sufix = 'select_snp'
        }
    }

    Array[File] snp_to_filter = select_first([merge_snp_vcfs.vcf_output, vcf_to_filter])

    if (if_select_indel) {
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.select_indel {
                input:
                    input_vcfs = vcf_to_filter,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals,
                    disk_size_gb = disk_size_gb_indel
            }        
        }

        call tasks.merge_vcfs as merge_indel_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = select_indel.vcf_indel,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'indel'
        }

        call tasks.bcftools_stats as select_indel_stats{
            input:
                input_vcfs = merge_indel_vcfs.vcf_output,
                output_basename = output_vcf_basename,
                docker_img = docker_img_bcftools,
                sufix = 'select_indel'
        }
    }

    Array[File] indel_to_filter = select_first([merge_indel_vcfs.vcf_output, vcf_to_filter])

    if (if_filter_snp){
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.filter_snp {
                input:
                    snp_to_filter = snp_to_filter,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals,
                    disk_size_gb = disk_size_gb_filter_snp
            }
        }

        call tasks.merge_vcfs as merge_filtered_snp_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = filter_snp.filtered_snp_vcf,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'snp.filtered'
        }

        call tasks.bcftools_stats as snp_after_filter_stat {
            input:
                input_vcfs = merge_filtered_snp_vcfs.vcf_output,
                output_basename = output_vcf_basename,
                docker_img = docker_img_bcftools,
                sufix = 'snp_filtered'
        }

        call tasks.plink_filter {
            input:
                input_vcf = merge_filtered_snp_vcfs.vcf_output,
                output_vcf_basename = output_vcf_basename,
                docker_img = docker_img_plink 
        }
    }

    if (if_filter_indel){
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.filter_indel {
                input:
                    indel_to_filter = indel_to_filter,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals,
                    disk_size_gb = disk_size_gb_filter_indel
            }
        }

        call tasks.merge_vcfs as merge_filtered_indel_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = filter_indel.filtered_indel_vcf,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'indel.filtered'
        }

        call tasks.bcftools_stats as indel_after_filter_stat {
            input:
                input_vcfs = merge_filtered_indel_vcfs.vcf_output,
                output_basename = output_vcf_basename,
                docker_img = docker_img_bcftools,
                sufix = 'indel_filtered'
        }

    }
    

    output {
        File vcf_stat = bcftools_stats.vcf_stat

        Array[File]? select_snp_to_filter = merge_snp_vcfs.vcf_output
        File? select_snp_stat = select_snp_stats.vcf_stat
        Array[File]? filtered_snp = merge_filtered_snp_vcfs.vcf_output
        File? filtered_snp_stat = snp_after_filter_stat.vcf_stat

        Array[File]? plink_file = plink_filter.plink_file
        Array[File]? plink_filter_file = plink_filter.plink_filter_file
        Array[File]? plink_filter_vcf = plink_filter.plink_filter_vcf

        Array[File]? select_indel_to_filter = merge_indel_vcfs.vcf_output
        File? select_indel_stat = select_indel_stats.vcf_stat
        Array[File]? filtered_indel = merge_filtered_indel_vcfs.vcf_output
        File? filtered_indel_stat = indel_after_filter_stat.vcf_stat
    }

    
}
