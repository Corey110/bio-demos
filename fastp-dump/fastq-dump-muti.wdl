version 1.0

workflow sra2fastq {
    input {
        String sample_name
        File sra_file_path
        String docker_img = 'cr-cn-beijing.volces.com/popgenomics/sra-tools:latest'
    }
    Array[File] sra_files = read_lines(sra_file_path)
    if (length(sra_files) == 1){
        File sra_file = sra_files[0]
        call fastq_dump as single_fastq_dump{
            input:
                sample_name = sample_name,
                sra_file = sra_file,
                docker_img = docker_img
        }
    }
    if (length(sra_files) > 1){
        scatter (idx in range(length(sra_files))) {
            call fastq_dump as multi_fastq_dump {
                input:
                    sample_name=sample_name,
                    sra_file=sra_files[idx],
                    docker_img=docker_img
            }
        }

        call merge_fastq {
            input:
                sample_name=sample_name,
                fastq1_files=multi_fastq_dump.fastq1,
                fastq2_files=multi_fastq_dump.fastq2,
                docker_img=docker_img
        }

    }

    output {
        File fastq1 = select_first([single_fastq_dump.fastq1, merge_fastq.fastq1])
        File fastq2 = select_first([single_fastq_dump.fastq2, merge_fastq.fastq2])
    }
}

task fastq_dump {
    input {
        String sample_name
        File sra_file

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(4 * size(sra_file, "GB")) + 20
        String docker_img
    }
    
    String out_fastq1 = sample_name + "_1.fastq.gz"
    String out_fastq2 = sample_name + "_2.fastq.gz"

    command <<<
        set -euo pipefail
        
        RUN=$(basename ~{sra_file} | awk -F. '{print $1}')

        fastq-dump --split-files --gzip ~{sra_file}

        mv ${RUN}_1.fastq.gz ~{out_fastq1}
        mv ${RUN}_2.fastq.gz ~{out_fastq2}

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

task merge_fastq {
    input {
        String sample_name
        Array[File] fastq1_files
        Array[File] fastq2_files

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(3 * (size(fastq1_files, "GB") + size(fastq2_files, "GB"))) + 20
        String docker_img
    }

    command <<<
        cat ~{sep=" " fastq1_files} > ~{sample_name}_1.fastq.gz &
        cat ~{sep=" " fastq2_files} > ~{sample_name}_2.fastq.gz &
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File fastq1 = "~{sample_name}_1.fastq.gz"
        File fastq2 = "~{sample_name}_2.fastq.gz"
    }
}