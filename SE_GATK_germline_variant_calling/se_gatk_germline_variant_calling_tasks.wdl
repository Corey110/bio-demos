version 1.0
# 动态interval版
task fastp_clean {
    input {
        File fastq1
        String sample_id

        Boolean phred64 
        Boolean fix_mgi_id

        String? adapter_sequence
        String? adapter_sequence_r2

        # specify how many reads/pairs to be processed. Default 0 means process all reads.
        Int? reads_to_process 

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(2 * ceil(size(fastq1, 'GB')) /10) * 10 + 40
        String docker_img
    }

    # reporting options
    String json = sample_id + "_fastp.json"
    String html = sample_id + "_fastp.html"
    String report_title = "\'~{sample_id} fastp report\'"

    String out1 = sample_id + '_1.clean.fq.gz'

    command <<<
        # basic command
        fastp \
            --in1 ~{fastq1} \
            --out1 ~{out1} \
            --json ~{json} \
            --html ~{html} \
            --report_title ~{report_title} \
        
        # options 
            ~{ true="--phred64 " false="" phred64 } \
            ~{ "--reads_to_process " + reads_to_process } \
            ~{ true="--fix_mgi_id " false="" fix_mgi_id } \
            ~{ "--adapter_sequence " + adapter_sequence } \
            ~{ "--adapter_sequence_r2 " + adapter_sequence_r2 }

    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }
    
    output {
        File cleaned_fq = out1
        Pair[File,File] fastp_report = (json, html)
    }
}

task bwa_align {
    input {
        File fastq1

        Array[File] genome_indexes

        String sample_id

        Int cpu
        Int memory
        Int disk_size_gb = 50
        String docker_img
    }
    String reads_group = "'@RG\\tID:~{sample_id}\\tSM:~{sample_id}\\tPL:ILLUMINA'"
    String sorted_bam = sample_id + '.sorted.bam'
    String sorted_bam_index = sorted_bam + '.bai'

    command <<<
       cpu_cores=$(nproc)

       bwa mem \
            -Y -K 100000000 \
            -t $cpu_cores \
            -R ~{reads_group} \
            ~{genome_indexes[0]} \
            ~{fastq1} \
       | samtools sort -@ $cpu_cores -o ~{sorted_bam}

       samtools index ~{sorted_bam}

    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File,File] sorted_bam_output = (sorted_bam, sorted_bam_index)
    }
}

task mosdepth {
    input {
        String sample_id
        Pair[File, File] aligned_bam

        Int cpu = 2
        Int memory = 4
        Int disk_size_gb = 40
        String docker_img
    }
    String mosdepth_global_dist = sample_id + '.mosdepth.global.dist.txt'
    String mosdepth_summary = sample_id + '.mosdepth.summary.txt'
    String mosdepth_per_base = sample_id + '.per-base.bed.gz'
    String mosdepth_per_base_csi = sample_id + '.per-base.bed.gz.csi'
    String sammtools_flagstat = sample_id + '.flagstat'


    command <<< 
        cpu_cores=$(nproc)
        mosdepth -t $cpu_cores ~{sample_id} ~{aligned_bam.left}
        samtools flagstat ~{aligned_bam.left} > ~{sammtools_flagstat}
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Array[File] mosdepth_res = [mosdepth_global_dist, mosdepth_summary, mosdepth_per_base, mosdepth_per_base_csi]
        File flagstate_res = sammtools_flagstat
    }
}

task dedup {
    input {
        String sample_id
        Pair[File, File] sorted_bam

        # Resource
        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(2 * ceil(size(sorted_bam.left, 'GB')) /10) * 10 + 40
        String docker_img
    }

    String dedup_bam = sample_id + ".deduplicated.bam"
    String dedup_bam_index = sample_id + ".deduplicated.bai"
    String dedup_metrics = sample_id + ".deduplicated.metrics"

    command <<<
        gatk MarkDuplicates \
            -I ~{sorted_bam.left} \
            -O ~{dedup_bam} \
            -M ~{dedup_metrics} \
            --CREATE_INDEX true
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File,File] dedup_bam_output = (dedup_bam, dedup_bam_index)
        File dedup_metrics_output = dedup_metrics
    }
}

task split_intervals {
    input {
        Int scatter_count
        Array[File] genome_indexes

        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = 40
        String docker_img
    }

    command <<<
        gatk SplitIntervals \
            -R ~{genome_indexes[0]} \
            --scatter-count ~{scatter_count} \
            -O intervals 

    >>>

    runtime{
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Array[File] scatter_interval_list_file = glob("intervals/*.interval_list")
    }
}

task haplotype_caller {
    input {
        String sample_id
        Pair[File,File] input_bam

        Array[File] genome_indexes
        File calling_region 
        
        # Resource
        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = 50
        String docker_img


    }

    String gvcf_name = sample_id + "_hc.g.vcf.gz"    

    command <<<
        gatk HaplotypeCaller \
            -R ~{genome_indexes[0]} \
            -I ~{input_bam.left} \
            -L ~{calling_region} \
            --ERC GVCF \
            -O ~{gvcf_name}
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        File hc_block_gvcf = gvcf_name
    }
}

task merge_vcfs {

    input {
        String sample_id
        Array[File] vcfs
        
        # Resource
        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(2 * ceil(size(vcfs, 'GB')) /10) * 10 + 40
        String docker_img
    }

    String gvcf = sample_id + ".g.vcf.gz"
    String gvcf_idx = sample_id + ".g.vcf.gz.tbi"

    command <<<
        gatk MergeVcfs \
            -I ~{sep=" -I " vcfs} \
            -O ~{gvcf}
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File, File] gvcf_output = (gvcf, gvcf_idx)
    }
}

task genotype_gvcfs {
    input {
        String sample_id
        Pair[File, File] gvcf
        
        Array[File] genome_indexes
       
        # Resource
        Int cpu = 2
        Int memory = 4
        Int disk_size_gb = ceil(2 * ceil(size(gvcf.left, 'GB')) /10) * 10 + 40
        String docker_img
    }

    String vcf = sample_id + ".vcf.gz"
    String vcf_idx = sample_id + ".vcf.gz.tbi"


    command <<<
        gatk GenotypeGVCFs \
            -R ~{genome_indexes[0]} \
            -V ~{gvcf.left} \
            -O ~{vcf} 
    >>>

    runtime {
        cpu: "~{cpu}"
        memory: '~{memory} GB'
        disk: "~{disk_size_gb} GB"
        docker: docker_img
    }

    output {
        Pair[File, File] vcf_output = (vcf, vcf_idx)
    }

}