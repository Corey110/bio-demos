version 1.0
 
import "./imputation_tasks.wdl" as tasks

workflow imputation {
    input {
        # 直接输入多样本VCF需要将其index文件一起作为Array输入
        # 单样本bcftools merge 需要将VCF和index的Array分开输入
        Array[File]? multi_sample_vcf
        Array[File]? single_sample_vcfs
        Array[File]? single_sample_vcf_indices
        
        Array[File] genome_indexes

        String output_vcf_basename

        # Int? disk_size_gb_snp
        # Int? disk_size_gb_indel
        # Int? disk_size_gb_filter

        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
        String docker_img_bcftools = 'cr-cn-beijing.volces.com/popgenomics/bcftools:1.9'
        String docker_img_plink = 'cr-cn-beijing.volces.com/popgenomics/plink:1.9'
        String docker_img_shapeit4 = 'cr-cn-beijing.volces.com/popgenomics/shapeit:4.2.2'
        String docker_img_beagle = 'cr-cn-beijing.volces.com/popgenomics/beagle:5.2'
        String docker_img_germline_tools = 'cr-cn-beijing.volces.com/popgenomics/germline_tools:v3'

        Boolean if_select_snp
        Boolean if_select_indel
        Boolean if_filter_snp

        File ref_panel
        Boolean choose_beagle

        Int scatter_count

        File contig_to_phasing
    }

    if (defined(single_sample_vcfs)) {
        call tasks.bcftools_merge_vcf {
            input:
                input_vcfs = select_first([single_sample_vcfs]),
                input_vcf_indices = select_first([single_sample_vcf_indices]),
                output_vcf_basename = output_vcf_basename,
                docker_img = docker_img_bcftools
        }
    }

    Array[File] total_vcf = select_first([multi_sample_vcf, bcftools_merge_vcf.output_vcf])

    call tasks.bcftools_stats {
        input:
            input_vcfs = total_vcf,
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
                    input_vcfs = total_vcf,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals
                    # disk_size_gb = select_first([disk_size_gb_snp, 50])
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

    Array[File] vcf_to_impute = select_first([merge_snp_vcfs.vcf_output, total_vcf])

    if (if_select_indel) {
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.select_indel {
                input:
                    input_vcfs = total_vcf,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    docker_img = docker_img_gatk4,
                    intervals = intervals
                    # disk_size_gb = select_first([disk_size_gb_indel, 50])
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
    if (if_filter_snp){
        scatter (intervals in split_intervals.scatter_interval_list_file) {
            call tasks.filter_variation {
                input:
                    vcf_to_impute = vcf_to_impute,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    intervals = intervals,
                    docker_img = docker_img_gatk4
                    # disk_size_gb = select_first([disk_size_gb_filter, 50])
            }
        }

        call tasks.merge_vcfs as merge_filtered_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = filter_variation.filtered_snp_vcf,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'snp.filtered'
        }

        call tasks.bcftools_stats as vcf_stat_after_filter {
            input:
                input_vcfs = merge_filtered_vcfs.vcf_output,
                output_basename = output_vcf_basename,
                docker_img = docker_img_bcftools,
                sufix = 'snp_filtered'
        }

        call tasks.plink_filter {
            input:
                input_vcf = merge_filtered_vcfs.vcf_output,
                output_vcf_basename = output_vcf_basename,
                docker_img = docker_img_plink  
        }
    }

    Array[File] vcf_to_phasing = select_first([plink_filter.plink_filter_vcf, vcf_to_impute])

    Array[String] scatter_contig_to_phasing = read_lines(contig_to_phasing)
    scatter (contig in scatter_contig_to_phasing){
        call tasks.shapit4_phasing {
            input:
                contig = contig,
                input_vcf = vcf_to_phasing,
                output_vcf_basename = output_vcf_basename,
                docker_img = docker_img_shapeit4,
                genome_indexes = genome_indexes
        }

        call tasks.reheader {
            input:
                contig = contig,
                vcf_to_reheader = shapit4_phasing.phased_vcf,
                output_vcf_basename = output_vcf_basename,
                genome_indexes = genome_indexes,
                docker_img = docker_img_germline_tools
        }
    }

    call tasks.merge_vcfs as merge_phased_vcfs {
        input:
            output_vcf_basename = output_vcf_basename,
            vcfs = reheader.phased_vcf,
            docker_img = docker_img_gatk4,
            output_vcf_sufix = 'snp.filtered.plink.phased'
    }

    if (choose_beagle) {
        scatter (contig in scatter_contig_to_phasing){
            call tasks.beagle_impute {
                input:
                    contig = contig,
                    input_vcf = merge_phased_vcfs.vcf_output,
                    output_vcf_basename = output_vcf_basename,
                    genome_indexes = genome_indexes,
                    ref_panel = ref_panel,
                    docker_img = docker_img_beagle
            }
        }

        call tasks.merge_vcfs as merge_imputed_vcfs {
            input:
                output_vcf_basename = output_vcf_basename,
                vcfs = beagle_impute.imputed_vcf,
                docker_img = docker_img_gatk4,
                output_vcf_sufix = 'snp.filtered.plink.phased.beagle'
        }
    }

    output {
        Array[File]? vcf_snp = select_snp.vcf_snp
        Array[File]? vcf_indels = merge_indel_vcfs.vcf_output

        File total_vcf_stat = bcftools_stats.vcf_stat
        File? snp_selected_stat = select_snp_stats.vcf_stat
        File? indel_selected_stat = select_indel_stats.vcf_stat

        Array[File]? hard_filtered_vcf = merge_filtered_vcfs.vcf_output
        File? hard_filtered_vcf_stat = vcf_stat_after_filter.vcf_stat

        Array[File]? plink_file = plink_filter.plink_file
        Array[File]? plink_filter_file = plink_filter.plink_filter_file
        Array[File]? plink_filter_vcf = plink_filter.plink_filter_vcf

        Array[File] phased_vcf = merge_phased_vcfs.vcf_output
        
        Array[File]? beagle_impute_vcf = merge_imputed_vcfs.vcf_output
    }

    
}
