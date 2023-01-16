version 1.0

import './gatk_joint_calling_DBI_tasks.wdl' as tasks

workflow JointGenotyping {
    input {
        Array[File] single_sample_gvcfs
        Array[File] single_sample_gvcf_indices

        Array[File] genome_indexes
        File? interval_list
        Int scatter_count

        String output_vcf_basename

        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
    }
    
    call tasks.split_intervals {
        input:
            scatter_count = scatter_count,
            genome_indexes = genome_indexes,
            intervals = interval_list,
            docker_img = docker_img_gatk4
    }

    scatter (intervals in split_intervals.scatter_interval_list_file) {
        call tasks.ImportGVCFs {
            input:
                single_sample_gvcfs = single_sample_gvcfs,
                single_sample_gvcf_indices = single_sample_gvcf_indices,
                workspace_dir_name = output_vcf_basename,
                intervals = intervals,
                docker_img = docker_img_gatk4
        }

        call tasks.GenotypeGVCFs {
            input:
                workspace_tar = ImportGVCFs.output_genomicsdb,
                intervals = intervals,
                output_vcf_filename = output_vcf_basename,
                genome_indexes = genome_indexes,
                docker_img = docker_img_gatk4
        }
    }

    call tasks.merge_vcfs {
        input:
            output_vcf_basename = output_vcf_basename,
            output_vcf_sufix = 'importdb',
            vcfs = GenotypeGVCFs.output_vcf,
            docker_img = docker_img_gatk4
    }
    
    output {
        Array[File] joint_vcf = merge_vcfs.vcf_output
    }
}