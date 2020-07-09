version 1.0

workflow simple_variant_call {
    input {
        String sample_id = "NIST7035"
        String input_prefix = "s3://aws-batch-genomics-shared/secondary-analysis/example-files/fastq"

        File reference = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
        String ref_name = basename(reference, ".fasta")
        
        # Array[Int] chromosome_ids = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
        Array[Int] chromosome_ids = [1,2]
        Array[String] chromosomes = prefix("chr", chromosome_ids)

        # WDL does not support directory globbing for inputs.
        # You need to specify each file you need explicitly.
        # See:
        # https://gatkforums.broadinstitute.org/wdl/discussion/8370/input-folders-for-wdl
        # https://github.com/broadinstitute/cromwell/issues/2269
        # https://github.com/openwdl/wdl/issues/289

        Array[File] reference_indices = [
            "${reference}.fai",
            
            "${reference}.amb",
            "${reference}.ann",
            "${reference}.bwt",
            "${reference}.pac",
            "${reference}.sa",

            "${reference}.64.alt",
            "${reference}.64.amb",
            "${reference}.64.ann",
            "${reference}.64.bwt",
            "${reference}.64.pac",
            "${reference}.64.sa"
        ]
        Array[File] fastqs = [
            "${input_prefix}/${sample_id}_R1_trim_samp-0p1.fastq.gz",
            "${input_prefix}/${sample_id}_R2_trim_samp-0p1.fastq.gz"
        ]
    }

    call bwa_mem {
        input:
            sample_id=sample_id,
            fastqs=fastqs,
            reference=reference,
            reference_indices=reference_indices
    }

    call samtools_sort {
        input:
            sample_id=sample_id,
            sam_file=bwa_mem.sam_file
    }

    call samtools_index {
        input:
            sample_id=sample_id,
            bam_file=samtools_sort.bam_file
    }

    scatter (chromosome in chromosomes) {
        call bcftools_mpileup {
            input:
                sample_id=sample_id,
                chromosome=chromosome,
                reference=reference,
                reference_indices=reference_indices,
                bam_file=samtools_sort.bam_file,
                bai_file=samtools_index.bai_file
        }

        call bcftools_call {
            input:
                sample_id=sample_id,
                chromosome=chromosome,
                mpileup_file=bcftools_mpileup.mpileup_file
        }
    }
    
    output {
        File bam_file = samtools_sort.bam_file
        File bai_file = samtools_index.bai_file
        Array[File] vcf_files = bcftools_call.vcf_file
    }
}

task bwa_mem {
    input {
        String sample_id
        File reference
        Array[File] reference_indices
        Array[File] fastqs
    }

    runtime {
        docker: "biocontainers/bwa:v0.7.15_cv4"
        memory: "64GB"
        cpu: 8
        disks: "local-disk"
    }
    command {
        bwa mem -t 16 \
            ${reference} \
            ${sep=' ' fastqs} \
            > ${sample_id}.sam
    }
    output {
        File sam_file = "${sample_id}.sam"
    }
}

task samtools_sort {
    input {
        String sample_id
        File sam_file
    }
    runtime {
        docker: "biocontainers/samtools:v1.7.0_cv4"
        memory: "32GB"
        cpu: 8
        disks: "local-disk"
    }
    command {
        samtools sort \
            -@ 16 \
            -o ${sample_id}.bam \
            ${sam_file}
    }
    output {
        File bam_file = "${sample_id}.bam"
    }
}

task samtools_index {
    input {
        String sample_id
        File bam_file
    }
    runtime {
        docker: "biocontainers/samtools:v1.7.0_cv4"
        memory: "32GB"
        cpu: 8
        disks: "local-disk"
    }
    command {
        samtools index \
            ${bam_file} \
            ${sample_id}.bam.bai
    }
    output {
        File bai_file = "${sample_id}.bam.bai"
    }
    
}

task bcftools_mpileup {
    input {
        String sample_id
        String chromosome
        File reference
        Array[File] reference_indices
        File bam_file
        File bai_file
        
        String bam_file_ = basename(bam_file)
    }
    runtime {
        docker: "biocontainers/bcftools:v1.5_cv3"
        memory: "32GB"
        cpu: 8
        disks: "local-disk"
    }
    command {
        ln -s ${bam_file}
        ln -s ${bai_file}

        ls -l

        bcftools mpileup \
            --threads 16 \
            -Oz \
            -r ${chromosome} \
            -f ${reference} \
            ${bam_file_} \
            > ${sample_id}.${chromosome}.mpileup.gz
    }
    output {
        File mpileup_file = "${sample_id}.${chromosome}.mpileup.gz"
    }
}

task bcftools_call {
    input {
        String sample_id
        String chromosome
        File mpileup_file
    }
    runtime {
        docker: "biocontainers/bcftools:v1.5_cv3"
        memory: "32GB"
        cpu: 8
        disks: "local-disk"
    }
    command {
        bcftools call \
            -m \
            --threads 16 \
            -Oz \
            -t ${chromosome} \
            -o ${sample_id}.${chromosome}.vcf.gz \
            ${mpileup_file}
    }
    output {
        File vcf_file = "${sample_id}.${chromosome}.vcf.gz"
    }
}