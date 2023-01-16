version 1.0

workflow joint_calling{
    input{
        Array[File] single_sample_vcfs
        Array[File] single_sample_vcf_indices

        Array[File] genome_indexes

        String output_vcf_basename

        String docker_img_gatk4 = 'cr-cn-beijing.volces.com/popgenomics/gatk:4.3.0.0'
    }

    call Combine_GVCFs {
        input:
            single_sample_vcfs = single_sample_vcfs,
            single_sample_vcf_indices = single_sample_vcf_indices,
            genome_indexes = genome_indexes,
            output_vcf_basename = output_vcf_basename,
            docker_img = docker_img_gatk4
    }

    call Genotype_GVCFs {
        input:
            sample_id = output_vcf_basename,
            gvcf = Combine_GVCFs.combined_gvcf,
            genome_indexes = genome_indexes,
            docker_img = docker_img_gatk4
    }

    output {
        Pair[File, File] combined_gvcf = Combine_GVCFs.combined_gvcf
        Pair[File, File] combined_vcf = Genotype_GVCFs.vcf_output
    } 
}

task Combine_GVCFs {
    input {
        Array[File] single_sample_vcfs
        Array[File] single_sample_vcf_indices
        Array[File] genome_indexes
        String output_vcf_basename

        Int cpu = 4
        Int memory = 32
        Int disk_size_gb = ceil(3 * ceil(size(single_sample_vcfs, 'GB')) /10) * 10
        String docker_img
    }

    command <<<
        gatk CombineGVCFs \
            --java-options "-Xmx~{memory}g" \
            -R ~{genome_indexes[0]} \
            --variant ~{sep=" --variant " single_sample_vcfs} \
            -O ~{output_vcf_basename}.g.vcf.gz
    >>>

    runtime{
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disks: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File, File] combined_gvcf = ('~{output_vcf_basename}.g.vcf.gz', '~{output_vcf_basename}.g.vcf.gz.tbi')
    }
}

task Genotype_GVCFs {
    input {
        String sample_id
        Pair[File, File] gvcf
        
        Array[File] genome_indexes
       
        # Resource
        Int cpu = 4
        Int memory = 16
        Int disk_size_gb = ceil(2 * ceil(size(gvcf.left, 'GB')) /10) * 10 + 20
        String docker_img
    }

    String vcf = sample_id + ".vcf.gz"
    String vcf_idx = sample_id + ".vcf.gz.tbi"


    command <<<
        gatk GenotypeGVCFs \
            --java-options "-Xmx~{memory}g" \
            -R ~{genome_indexes[0]} \
            -V ~{gvcf.left} \
            -O ~{vcf} 
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disks: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File, File] vcf_output = (vcf, vcf_idx)
    }
}