version 1.0
#WORKFLOW DEFINITION

task HaplotypeCaller_GVCF {
    input {
        String java_opt = " -Xmx2048m "
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_bam
        File input_bai
        String docker_image = "broadinstitute/gatk:4.3.0.0"
        Int mem_size = 20
        Int cpu_size = 8
        Int disk_size = 50
        String sample_name

    }

    command {
        gatk --java-options ${java_opt} HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${sample_name}.gvcf \
            -ERC GVCF
    }
    
    runtime {
        docker: docker_image
        cpu: cpu_size
        memory: mem_size + " GB"
        disk: disk_size + " GB" 
    }

    output {
        File output_gvcf = "${sample_name}.gvcf"
    } 
}

task CramToBamTask{
    input{
        #Command parameters
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_cram
        String sample_name

        # Runtime parameters
        # Int addtional_disk_size = 20 
        Int machine_mem_size = 4 
        Int disk_size = 50
        String docker_image = "biocontainers/samtools:v1.7.0_cv4"
    }
       # Float output_bam_size = size(input_cram, "GB") / 0.60
       # Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
       # Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + addtional_disk_size
#Calls samtools view to do the conversion
    command {
        set -eo pipefail

        samtools view -h -T ~{ref_fasta} ~{input_cram} |
        samtools view -b -o ~{sample_name}.bam -
        samtools index -b ~{sample_name}.bam
        mv ~{sample_name}.bam.bai ~{sample_name}.bai
    }

    runtime {
        docker: docker_image
        memory: machine_mem_size + " GB"
        disk: disk_size + " GB"
    }
    #Outputs a BAM and BAI with the same sample name
    output {
     File outputBam = "~{sample_name}.bam"
     File outputBai = "~{sample_name}.bai"
    }
}


workflow Cram2GVCF {
    input {
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File input_cram
        String sample_name
    }
    #converts CRAM to SAM to BAM and makes BAI
    call CramToBamTask {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            input_cram = input_cram,
            sample_name = sample_name
    }

    #call gatk function
    call HaplotypeCaller_GVCF{
        input:
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            sample_name = sample_name,
            input_bam = CramToBamTask.outputBam,
            input_bai = CramToBamTask.outputBai
    }



    #Outputs Bam, Bai, and validation report to the FireCloud data model
    output {
        File outputBam = CramToBamTask.outputBam
        File outputBai = CramToBamTask.outputBai
        File outputVcf = HaplotypeCaller_GVCF.output_gvcf
    }
}
