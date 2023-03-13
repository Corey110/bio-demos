version 1.0

workflow ChIP_seq_atac {
    meta {
        description: "ChIP-seq ATAC pipeline based on CASMAP "
        author: "Bio-OS"
        email:  ""
    }

    input {
        File REFZIP
        String sampleName
        File FASTQ1
        File? FASTQ2
        String Genome = "mm"
        String BROAD = "false"
        File BLACKLIST = "s3://bioos-wcetqgmleig4buued3380/ENCFF001TDO.merge.bed"
        String PAIR = "true"
        String IgnoreForNormalization = "chrX chrY chrM"
        String AnnotationGenome = "mm10"   
    }

    call fastqc_trim_mapping_filter{
        input:  FASTQ1 = FASTQ1,
                FASTQ2 = FASTQ2,
                REFZIP = REFZIP,
                AnnotationGenome = AnnotationGenome,
                sampleName = sampleName
    }

    call call_peaks{
        input:  BAM = fastqc_trim_mapping_filter.filter_bam,
                BAI = fastqc_trim_mapping_filter.filter_bai,
                Genome = Genome,
                BROAD = BROAD,
                sampleName = sampleName
    }

    call annotation{
        input:  peak = call_peaks.peak,
                Script = call_peaks.model,
                AnnotationGenome = AnnotationGenome,
                sampleName = sampleName
    }

    call motif{
        input:  peak = call_peaks.peak,
                AnnotationGenome = AnnotationGenome,
                sampleName = sampleName 
    }

    call normalization{
        input:  BAM = fastqc_trim_mapping_filter.filter_bam,
                BAI = fastqc_trim_mapping_filter.filter_bai,
                BLACKLIST = BLACKLIST,
                PAIR = PAIR,
                IgnoreForNormalization = IgnoreForNormalization,
                sampleName = sampleName 
    }

    call chip_qc{
        input:  BAM = fastqc_trim_mapping_filter.filter_bam,
                BAI = fastqc_trim_mapping_filter.filter_bai,
                sampleName = sampleName   
    }

    
}

task fastqc_trim_mapping_filter{
    meta {
        description: "fastqc and trim raw fq data"
    }
    input {
        File FASTQ1
        File? FASTQ2
        File REFZIP
        String AnnotationGenome

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "gzlab-cn-beijing.cr.volces.com/atac/chip-atac-qc:v1.0"
        String disk = "200 GB"

    }
    command {
        set -e
        

        # referecnce
        unzip ${REFZIP} -d .
        
        # fastqc
        fastqc -o . ${FASTQ1} ${FASTQ2}

        if [ -f ${default="foobar" FASTQ2} ]; then
            # trimming
            trim_galore --paired --fastqc --gzip ${FASTQ1} ${FASTQ2}

            # mapping
            #FASTQ_TIM1=`ls ${sampleName}*1*.gz`
            #FASTQ_TIM2=`ls ${sampleName}*2*.gz` 
            FASTQ_TIM1=`ls *1*.gz`
            FASTQ_TIM2=`ls *2*.gz` 

            bowtie2 \
                -p ${NUM_THREAD} \
                -q \
                -x REF/${AnnotationGenome} \
                -1 $FASTQ_TIM1 \
                -2 $FASTQ_TIM2 \
                |samtools view -bS - -@ ${NUM_THREAD} \
                |samtools sort - -@ ${NUM_THREAD} \
                -o ${sampleName}_sort.bam
            samtools index ${sampleName}_sort.bam
            samtools flagstat ${sampleName}_sort.bam > ${sampleName}_sort.stat

            # filter
            samtools view -q 10 \
                -b \
                -@ ${NUM_THREAD} \
                -f 2 \
                -U ${sampleName}_fail_qc.bam \
                ${sampleName}_sort.bam \
                > ${sampleName}_match.bam
        else
            # trimming
            trim_galore --fastqc --gzip ${FASTQ1}

            # mapping
            FASTQ_TIM=`ls ${sampleName}*gz`
                bowtie2 \
                -p ${NUM_THREAD} \
                -q \
                -x ${AnnotationGenome} \
                -U $FASTQ_TIM \
                |samtools view -bS - -@ ${NUM_THREAD} \
                |samtools sort - -@ ${NUM_THREAD} \
                -o ${sampleName}_sort.bam
            samtools index ${sampleName}_sort.bam
            samtools flagstat ${sampleName}_sort.bam > ${sampleName}_sort.stat

            #filter
            samtools view -q 2 \
                -b \
                -@ ${NUM_THREAD} \
                -U ${sampleName}_fail_qc.bam \
                ${sampleName}_sort.bam \
                > ${sampleName}_match.bam
        fi

        samtools index ${sampleName}_match.bam
        sambamba markdup \
            -r \
            -t ${NUM_THREAD} \
            ${sampleName}_match.bam \
            ${sampleName}_sort_filter.bam

        sambamba index ${sampleName}_sort_filter.bam
        rm -f ${sampleName}_match*
        samtools flagstat ${sampleName}_sort_filter.bam > ${sampleName}_sort_filter.stat

    }
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File filter_bam = "${sampleName}_sort_filter.bam"
        File filter_bai = "${sampleName}_sort_filter.bam.bai"
        File filter_bam_stat = "${sampleName}_sort_filter.stat"
        Array[File] zip_qcs = glob("*.zip")
        Array[File] txt_qcs = glob("*.txt")
        Array[File] html_qcs = glob("*.html")
    }
}

task call_peaks {
    meta {
        description: "call PEAKS"
    }
    input {
        File BAM
        File BAI
        String Genome = "mm"
        String BROAD = "false"

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "gzlab-cn-beijing.cr.volces.com/atac/macs2:v0.1"
        String disk = "200 GB"

    }
    command {
        set -e
        # source activate chipseq
        # source activate atac
        # source /root/miniconda3/etc/profile.d/conda.sh
        # conda activate atac 

        if [ "${BROAD}" == "false" ]; then
            macs2 callpeak -t ${BAM} -f BAM -g ${Genome} -n ${sampleName} \-q 0.01 --keep-dup all
            # 统一输出
            cp ${sampleName}_peaks.narrowPeak ${sampleName}.Peak
        elif [ "${BROAD}" == "true" ]; then
            macs2 callpeak -t ${BAM} -f BAM -g ${Genome} -n ${sampleName} --keep-dup all --broad --broad-cutoff 0.1

            cp ${sampleName}_peaks.broadPeak ${sampleName}.Peak
        else
            echo "Wrong Parameters"
        fi

        # conda deactivate
    }
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File peak = "${sampleName}.Peak"
        File model = "${sampleName}_model.r"
        File peakxls = "${sampleName}_peaks.xls"
        File summits = "${sampleName}_summits.bed" 
    }
}

task annotation{
    meta {
        description: "annotation"
    }
    input {
        File Script
        File peak
        String AnnotationGenome = "mm9"

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "r-base:4.2.2"
        String disk = "200 GB"

    }
    command <<<
        set -e

        awk -v \
            OFS="\t" \
            '{print $1,$2,$3,$4,$5}' \
            ~{peak} \
            > ~{sampleName}_tmp.bed

        Rscript ~{Script} \
            . \
            ~{sampleName}_tmp.bed \
            ~{AnnotationGenome}
    >>>
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File modelpdf = "${sampleName}_model.pdf"
    }

}

task motif {
    meta {
        description: "call motif"
    }
    input {
        File peak
        String AnnotationGenome = "mm10"

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "gzlab-cn-beijing.cr.volces.com/atac/chip_atac_motif:hg19"
        String disk = "200 GB"

    }
    command <<<
        set -e
       
        awk -v \
            OFS='\t' \
            '{print $4,$1,$2,$3,0}' \
            ~{peak} \
            > ~{sampleName}_tmp.bed

        annotatePeaks.pl \
            ~{sampleName}_tmp.bed \
            ~{AnnotationGenome} \
            > ~{sampleName}_annotated.txt

        findMotifsGenome.pl \
            ~{sampleName}_tmp.bed \
            ~{AnnotationGenome} \
            ~{sampleName}_motif/
    >>>
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File annotated = "${sampleName}_annotated.txt"
        Array[File] homerResults = glob("${sampleName}_motif/homerResults/*")
        Array[File] knownResults = glob("${sampleName}_motif/knownResults/*")
        File homerResults_html = "${sampleName}_motif/homerResults.html"
        File all_motifs = "${sampleName}_motif/homerMotifs.all.motifs"
        File motifs12 = "${sampleName}_motif/homerMotifs.motifs12"
        File motifs10 = "${sampleName}_motif/homerMotifs.motifs10"
        File motifs8 = "${sampleName}_motif/homerMotifs.motifs8"
        File knownResults_html = "${sampleName}_motif/knownResults.html"
        File knownResults_txt = "${sampleName}_motif/knownResults.txt"
        File seq_autonorm_tsv = "${sampleName}_motif/seq.autonorm.tsv"
        File motifFindingParameters = "${sampleName}_motif/motifFindingParameters.txt" 
    }

}

task normalization {
    meta {
        description: "normalization"
    }
    input {
        File BAM
        File BAI
        File BLACKLIST
        String PAIR
        String IgnoreForNormalization = "chrX chrY chrM"

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "gzlab-cn-beijing.cr.volces.com/atac/atac:v1.0"
        String disk = "200 GB"

    }
    command {
        set -e
        # source activate chipseq
        #source activate atac
        source /root/miniconda3/etc/profile.d/conda.sh
        conda activate atac 
        ## bamCoverage from deppTools
        if [ "${PAIR}" == "true" ]; then
            bamCoverage \
                --bam ${BAM} \
                -o ${sampleName}_normalize.bw \
                -bl ${BLACKLIST} \
                --binSize 50 \
                --normalizeUsing RPKM \
                --ignoreForNormalization ${IgnoreForNormalization} \
                --ignoreDuplicates \
                --extendReads \
                -of bigwig \
                -p max \
                --smoothLength 60
        elif [ "${PAIR}" == "false" ]; then
            bamCoverage \
                --bam ${BAM} \
                -o ${sampleName}_normalize.bw \
                -bl ${BLACKLIST} \
                --binSize 50 \
                --normalizeUsing RPKM \
                --ignoreForNormalization ${IgnoreForNormalization} \
                --ignoreDuplicates \
                -of bigwig \
                -p max \
                --smoothLength 60
        else
            echo "Wrong Parameters"
        fi

        conda deactivate
    }
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File bw = "${sampleName}_normalize.bw"
    }
}

task chip_qc {
    meta {
        description: "chip-qc"
    }
    input {
        File Script = "s3://bioos-wcetqgmleig4buued3380/chipseq_bamQC.py"
        File BAM
        File BAI

        String sampleName

        Int NUM_THREAD = 8
        String MEMORY = "50 GB"
        String docker = "gzlab-cn-beijing.cr.volces.com/atac/atac:v1.0"
        String disk = "200 GB"

    }
    command {
        set -e

        python ${Script} -i ${BAM} -c 8 -o ${sampleName}_sort_filter.tsv
    }
    runtime {
        docker: "${docker}" 
        cpu: "${NUM_THREAD}" 
        memory: "${MEMORY}"
        disk: "${disk}"
    }
    output {
        File tsv = "${sampleName}_sort_filter.tsv"
    }

}
