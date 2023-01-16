version 1.0

task split_intervals {
    input {
        Int scatter_count
        Array[File] genome_indexes
        File? intervals

        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = 40
        String docker_img
    }

    command <<<
        gatk SplitIntervals \
            -R ~{genome_indexes[0]} \
            --scatter-count ~{scatter_count} \
            -O intervals \
            ~{if defined(intervals) then "-L ~{intervals}" else ""}

    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Array[File] scatter_interval_list_file = glob("intervals/*.interval_list")
    }
}

task ImportGVCFs {
    input {
        Array[File] single_sample_gvcfs
        Array[File] single_sample_gvcf_indices

        String workspace_dir_name
        File intervals

        Int cpu = 4
        Int memory = 32
        Int disk_size_gb = ceil((size(single_sample_gvcfs, "GB") + size(single_sample_gvcf_indices, "GB")) * 2) + 20
        String docker_img
    }

    Int max_mem = (memory - 2) * 1000 

    command <<<
        set -euo pipefail
        
        gatk --java-options "-Xms8000m -Xmx~{max_mem}m" \
                GenomicsDBImport \
                --genomicsdb-workspace-path ~{workspace_dir_name} \
                --batch-size 50 \
                -L ~{intervals} \
                --reader-threads 5 \
                --merge-input-intervals \
                --consolidate \
                -V ~{sep = " -V " single_sample_gvcfs}

        tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
        
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: "~{memory} GB"
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        File output_genomicsdb = "~{workspace_dir_name}.tar"
    }

}

task GenotypeGVCFs {
    input {
        File workspace_tar
        File intervals

        String output_vcf_filename
        Array[File] genome_indexes

        Int cpu = 2
        Int memory = 16
        Int disk_size_gb = ceil(size(workspace_tar, "GB") * 5) + 20
        String docker_img
    }

    Int max_mem = (memory - 2) * 1000 

    command <<<
        set -euo pipefail

        tar -xf ~{workspace_tar}
        WORKSPACE=$(basename ~{workspace_tar} .tar)

        gatk --java-options "-Xms8000m -Xmx~{max_mem}m" \
            GenotypeGVCFs \
            -R ~{genome_indexes[0]} \
            -O ~{output_vcf_filename}.vcf.gz \
            -G StandardAnnotation -G AS_StandardAnnotation \
            --only-output-calls-starting-in-intervals \
            -V gendb://$WORKSPACE \
            -L ~{intervals} \
            --merge-input-intervals
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: "~{memory} GB"
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        File output_vcf = "~{output_vcf_filename}.vcf.gz"
        File output_vcf_index = "~{output_vcf_filename}.vcf.gz.tbi"
    }
}

task merge_vcfs {

    input {
        String output_vcf_basename
        String output_vcf_sufix
        Array[File] vcfs
        
        # Resource
        Int cpu = 2
        Int memory = 16
        Int disk_size_gb = 2 * ceil(size(vcfs, "GB")) + 50
        
        # docker 
        String docker_img    
    }

    String vcf = output_vcf_basename + ".merged." + output_vcf_sufix + ".vcf.gz"
    String vcf_idx = vcf + ".tbi"

    command <<<
        gatk MergeVcfs \
            -I ~{sep=" -I " vcfs} \
            -O ~{vcf}
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Array[File] vcf_output = [vcf, vcf_idx]
    }
}
