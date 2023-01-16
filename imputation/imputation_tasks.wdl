version 1.0

task bcftools_merge_vcf {
    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices
        String output_vcf_basename
        
        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil((3 * ceil(size(input_vcfs, "GB") + ceil(size(input_vcf_indices, "GB"))) + 20) /10) * 10 + 50
        String docker_img
    }

    command <<<
        bcftools merge ~{sep=' ' input_vcfs} -O z -o ~{output_vcf_basename}.vcf.gz
        bcftools index -t ~{output_vcf_basename}.vcf.gz
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: "~{memory} GB"
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        Array[File] output_vcf = ["~{output_vcf_basename}.vcf.gz", "~{output_vcf_basename}.vcf.gz.tbi"]
    }
}

task bcftools_stats {
    input {
        Array[File] input_vcfs
        String output_basename
        String sufix

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil((3 * ceil(size(input_vcfs, "GB"))) / 10) * 10 + 20
        String docker_img
    }

    command <<<
        bcftools stats ~{input_vcfs[0]} > ~{output_basename}.~{sufix}.stat
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: "~{memory} GB"
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File vcf_stat = '~{output_basename}.~{sufix}.stat'
    }
}

task select_snp {
    input {
        Array[File] input_vcfs
        String output_vcf_basename
        Array[File] genome_indexes
        File intervals

        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = ceil((1.2 * ceil(size(input_vcfs, "GB"))) / 10) * 10 + 20
        String docker_img
    }

    String out_snp = output_vcf_basename + ".snp.vcf.gz"
    String out_snp_index = out_snp + ".tbi"

    command <<<
        gatk SelectVariants \
            -R ~{genome_indexes[0]} \
            -V ~{input_vcfs[0]} \
            -select-type SNP \
            -O ~{out_snp} \
            -L ~{intervals}

    >>>

    runtime {
        cpu: '~{cpu}'
        memory: "~{memory} GB"
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File vcf_snp = out_snp
        File vcf_snp_index = out_snp_index
    }
}

task select_indel {
    input {
        Array[File] input_vcfs
        String output_vcf_basename
        Array[File] genome_indexes
        File intervals

        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = ceil((1.2 * ceil(size(input_vcfs, "GB"))) / 10) * 10 + 20
        String docker_img
    }

    String out_indel = output_vcf_basename + ".indel.vcf.gz"
    String out_indel_index = out_indel + ".tbi"

    command <<<
        gatk SelectVariants \
            -R ~{genome_indexes[0]} \
            -V ~{input_vcfs[0]} \
            -select-type INDEL \
            -O ~{out_indel} \
            -L ~{intervals}

    >>>

    runtime {
        cpu: '~{cpu}'
        memory: "~{memory} GB"
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File vcf_indel = out_indel
        File vcf_indel_index = out_indel_index
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
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        Array[File] scatter_interval_list_file = glob("intervals/*.interval_list")
    }
}


task filter_variation{
    input {
        Array[File] vcf_to_impute
        Array[File] genome_indexes
        String output_vcf_basename
        File intervals

        Int cpu = 1
        Int memory = 2
        Int disk_size_gb = ceil((1.2 * ceil(size(vcf_to_impute, "GB"))) / 10) * 10 + 20
        String docker_img
    }

    command <<<
        gatk VariantFiltration \
            -R ~{genome_indexes[0]} \
            -V ~{vcf_to_impute[0]} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            --cluster-size 2 \
            --cluster-window-size 5 \
            -O ~{output_vcf_basename}.snp.filtered.vcf.gz \
            -L ~{intervals}

        gatk SelectVariants \
            -R ~{genome_indexes[0]} \
            -V ~{output_vcf_basename}.snp.filtered.vcf.gz \
            --exclude-filtered true \
            -O ~{output_vcf_basename}.snp.filtered.removed.vcf.gz

    >>>

    runtime {
        cpu: '~{cpu}'
        memory: "~{memory} GB"
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File filtered_snp_vcf = '~{output_vcf_basename}.snp.filtered.removed.vcf.gz'
        File filtered_snp_vcf_index = '~{output_vcf_basename}.snp.filtered.removed.vcf.gz.tbi'
    }
}

task merge_vcfs {

    input {
        String output_vcf_basename
        String output_vcf_sufix
        Array[File] vcfs
        
        # Resource
        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(2 * ceil(size(vcfs, "GB")) / 10) * 10 + 20
        
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
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        Array[File] vcf_output = [vcf, vcf_idx]
    }
}

task plink_filter {
    input {
        Array[File] input_vcf
        String output_vcf_basename

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(2 * ceil(size(input_vcf, "GB")) / 10) * 10 + 20
        
        # docker 
        String docker_img    
    }

    command <<<
        plink --vcf ~{input_vcf[0]} \
            --allow-no-sex \
            --keep-allele-order \
            --allow-extra-chr \
            --make-bed \
            --out ~{output_vcf_basename}.plink

        plink --bfile ~{output_vcf_basename}.plink \
            --allow-no-sex \
            --keep-allele-order \
            --allow-extra-chr \
            --maf 0.05 \
            --geno 0.5 \
            --make-bed \
            --out ~{output_vcf_basename}.plink.filtered
        
        plink --bfile ~{output_vcf_basename}.plink.filtered \
            --allow-no-sex \
            --keep-allele-order \
            --allow-extra-chr \
            --out ~{output_vcf_basename}.plink.filtered.transform \
            --recode vcf-iid bgz

        bcftools index -t ~{output_vcf_basename}.plink.filtered.transform.vcf.gz
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        Array[File] plink_file = ['~{output_vcf_basename}.plink.bed', '~{output_vcf_basename}.plink.bim', '~{output_vcf_basename}.plink.fam', '~{output_vcf_basename}.plink.log']
        Array[File] plink_filter_file = ['~{output_vcf_basename}.plink.filtered.bed', '~{output_vcf_basename}.plink.filtered.bim', '~{output_vcf_basename}.plink.filtered.fam', '~{output_vcf_basename}.plink.filtered.log']
        Array[File] plink_filter_vcf = ['~{output_vcf_basename}.plink.filtered.transform.vcf.gz', '~{output_vcf_basename}.plink.filtered.transform.vcf.gz.tbi']
    } 
}

task shapit4_phasing {
    input {
        String contig
        Array[File] input_vcf
        String output_vcf_basename
        Array[File] genome_indexes

        Int cpu = 4
        Int memory = 16
        Int disk_size_gb = ceil(1.2 * ceil(size(input_vcf, "GB")) / 10) * 10 + 20
        String docker_img
    }

    command <<<
        shapeit4 --input ~{input_vcf[0]} \
            --region ~{contig} \
            --output ~{output_vcf_basename}.~{contig}.phased.vcf.gz \
            --log ~{output_vcf_basename}.~{contig}.phased.log \
            --thread ~{cpu}
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File phased_vcf = '~{output_vcf_basename}.~{contig}.phased.vcf.gz'
        File phase_log = '~{output_vcf_basename}.~{contig}.phased.log'
    }
}

task reheader {
    input {
        String contig
        File vcf_to_reheader
        String output_vcf_basename
        Array[File] genome_indexes

        Int cpu = 2
        Int memory = 8
        Int disk_size_gb = ceil(3 * ceil(size(vcf_to_reheader, "GB")) / 10) * 10 + 20
        String docker_img
    }

    command <<<

        zgrep '#' ~{vcf_to_reheader} | grep -v 'contig' > newheader

        bcftools reheader -h newheader ~{vcf_to_reheader} -o ~{output_vcf_basename}.phased.reheader.vcf.gz

        bcftools index -t ~{output_vcf_basename}.phased.reheader.vcf.gz

        gatk UpdateVCFSequenceDictionary \
            -V ~{output_vcf_basename}.phased.reheader.vcf.gz \
            -R ~{genome_indexes[0]} \
            --replace \
            -O ~{output_vcf_basename}.~{contig}.phased.vcf.gz
    >>>

    runtime {
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output {
        File phased_vcf = '~{output_vcf_basename}.~{contig}.phased.vcf.gz'
        File phased_vcf_index = '~{output_vcf_basename}.~{contig}.phased.vcf.gz.tbi'
    }
}

task beagle_impute {
    input{
        String contig
        File ref_panel
        Array[File] input_vcf
        String output_vcf_basename
        Array[File] genome_indexes

        Int cpu = 4
        Int memory = 16
        Int disk_size_gb = ceil(1.2 * ceil(size(input_vcf, "GB")) / 10) * 10 + 20
        String docker_img
    }

    command <<<
        beagle -Xmx~{memory}g gt=~{input_vcf[0]} \
            out=~{output_vcf_basename}.~{contig}.beagle.impute \
            ref=~{ref_panel} \
            chrom=~{contig}
        
        bcftools index -t ~{output_vcf_basename}.~{contig}.beagle.impute.vcf.gz

        gatk UpdateVCFSequenceDictionary \
            -V ~{output_vcf_basename}.~{contig}.beagle.impute.vcf.gz \
            -R ~{genome_indexes[0]} \
            --replace \
            -O ~{output_vcf_basename}.~{contig}.beagle.imputed.vcf.gz
    >>>

    runtime{
        cpu: '~{cpu}'
        memory: '~{memory} GB'
        disk: '~{disk_size_gb} GB'
        docker: docker_img
    }

    output{
        File imputed_vcf = '~{output_vcf_basename}.~{contig}.beagle.imputed.vcf.gz'
        File imputed_vcf_index = '~{output_vcf_basename}.~{contig}.beagle.imputed.vcf.gz.tbi'
        File impute_log = '~{output_vcf_basename}.~{contig}.beagle.impute.log'
    }   
}