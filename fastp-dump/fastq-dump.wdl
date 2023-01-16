version 1.0

workflow sra2fastq {
    input {
        String sample_id
        File sra_file
        String docker_img = 'cr-cn-beijing.volces.com/popgenomics/sra-tools:latest'
    }

    call fastq_dump {
        input:
            sample_id=sample_id,
            sra_file=sra_file,
            docker_img=docker_img
    }

    output {
        File fastq1 = fastq_dump.fastq1
        File fastq2 = fastq_dump.fastq2
    }
}

task fastq_dump {
    input {
        String sample_id
        File sra_file

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil((3 * ceil(size(sra_file, "GB"))) / 10) * 10 + 20
        String docker_img
    }
    
    String out_fastq1 = sample_id + "_1.fastq.gz"
    String out_fastq2 = sample_id + "_2.fastq.gz"

    command <<<
        fastq-dump --split-files --gzip ~{sra_file}
    >>>

    runtime {    
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File fastq1 = out_fastq1
        File fastq2 = out_fastq2
    }
}