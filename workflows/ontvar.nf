/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ontvar_pipeline'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
include { SNIFFLES } from '../modules/nf-core/sniffles/main'
include { CUTESV   } from '../modules/nf-core/cutesv/main'
include { SEVERUS  } from '../modules/nf-core/severus/main'
include { RENAME_VCF } from '../modules/local/rename_vcf/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_SNIFFLES } from '../modules/local/rename_vcf_headers/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_CUTESV   } from '../modules/local/rename_vcf_headers/main'
include { RENAME_VCF_HEADERS as RENAME_VCF_HEADERS_SEVERUS  } from '../modules/local/rename_vcf_headers/main'
include { JASMINESV as JASMINESV_SAMPLE } from '../modules/nf-core/jasminesv/main'
include { JASMINE_HEADER_FIX } from '../modules/local/jasmine_header_fix/main'
include { JASMINESV as JASMINESV_COHORT } from '../modules/nf-core/jasminesv/main'
include { FILTER_CHR } from '../modules/local/filter_chr/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_COHORT_RAW } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_COHORT    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_PER_SAMPLE_RAW    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_ANNOTSV as ANNOTSV_PER_SAMPLE    } from '../modules/nf-core/annotsv/annotsv/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../modules/nf-core/annotsv/installannotations/main'
include { UNTAR as UNTAR_ANNOTSV } from '../modules/nf-core/untar/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLERS          } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLER_MERGED    } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_CALLER_MERGED_FILTERED  } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_ANNOTATED } from '../modules/local/summarize_sv_counts/main'
include { SUMMARIZE_SV_COUNTS as SUMMARIZE_COHORT_FILTERED  } from '../modules/local/summarize_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_RAW_CALLERS          } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_CONSENSUS            } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_FILTERED             } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_COHORT_ANNOTATED     } from '../modules/local/plot_sv_counts/main'
include { PLOT_SV_COUNTS as PLOT_COHORT_FILTERED      } from '../modules/local/plot_sv_counts/main'
include { SVDB_QUERY as SVDB_QUERY_SAMPLE } from '../modules/nf-core/svdb/query/main'
include { SVDB_QUERY as SVDB_QUERY_COHORT } from '../modules/nf-core/svdb/query/main'
include { BCFTOOLS_VIEW as CALLER_SUPPORT_FILTER } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as AF_FILTER } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as AF_FILTER_COHORT } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_SORT as SORT_VCF } from '../modules/nf-core/bcftools/sort/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONTVAR {

    take:
        ch_samplesheet // channel: samplesheet read in from --input
        ch_output_dir // channel: output directory from --outdir
        reference
        annotsv_annotations

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

   

    ch_sample_info = ch_samplesheet
    // ch_sample_info contains: [[group_id:group_id, id:sample_id, status:case|control, input_type:bam|fastq], input_path]
    //                          [   [meta], input_path    ]
    

   ch_sample_info
    .branch { meta, input_path ->
        fastq_dir:  meta.input_type == 'fastq' && file(input_path).isDirectory()
        fastq_file: meta.input_type == 'fastq'
        bam:        meta.input_type == 'bam'
    }
     .set { ch_input }
    // ──────────────────────────────────────────────────────────────────────
    // FASTQ CONCATENATION (if multiple files per sample)
    // ──────────────────────────────────────────────────────────────────────
    
    // If meta.input_type == fastq -> do concat (if dir) + alignment

    // Add meta.single_end to cat fastq input to mimic nanopore fastq structure

    ch_input_fastq_concat = ch_input.fastq_dir
        .map { meta, input_path ->
            def files = file("${input_path}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true)
            
            return [meta + [single_end:true], files ]
            }

    // Concatenate FASTQs from directories
    CAT_FASTQ(ch_input_fastq_concat)


    

    // Merge concatenated and single-file FASTQs for minimap2
    minimap2_input = ch_input.fastq_file
        .mix(CAT_FASTQ.out.reads)



    // ──────────────────────────────────────────────────────────────────────
    // ALIGNMENT (minimap2 for long-read sequencing)
    // ──────────────────────────────────────────────────────────────────────

    MINIMAP2_ALIGN(
        minimap2_input,                                                          // Input 1: [meta, reads]
        Channel.value(tuple([id: "reference"], file(params.reference))),        // Input 2: [meta2, reference]
        Channel.value(true),                                                     // Input 3: bam_format (true for BAM output)
        Channel.value('bai'),                                                    // Input 4: bam_index_extension ('bai' or 'csi')
        Channel.value(false),                                                    // Input 5: cigar_paf_format (false, not needed for BAM)
        Channel.value(true) 
    )


    ch_bam = MINIMAP2_ALIGN.out.bam.mix(ch_input.bam)



    // ──────────────────────────────────────────────────────────────────────
    // Prepare SV input
    // ──────────────────────────────────────────────────────────────────────


    ch_cases = ch_bam
        .filter { meta, bam -> meta.status == 'case' }
    
    ch_controls = ch_bam
        .filter {meta, bam -> meta.status == 'control' }
        .map { meta, bam ->
            tuple(meta.group_id, meta, bam)
        }

    sniffles_cutesv_input = ch_cases
        .map { meta, bam ->
            def bai = file("${bam}.bai")
            tuple (meta, bam, bai)
        }
    
    ch_cases_grouped = ch_cases
        .map {meta, bam ->
            tuple(meta.group_id, meta, bam)
        }

    severus_input = ch_cases_grouped
        .join(ch_controls, by: 0, remainder: true)
        .map { it -> 
                // it[0] = group_id
                // it[1] = case_meta
                // it[2] = case_bam

                def meta = it[1]
                def case_bam = it[2]
                def case_bai = file("${case_bam}.bai")

                def ctrl_bam = (it.size() == 5) ? it[4] : null

                def input_ctrl_bam = ctrl_bam ?: []
                def input_ctrl_bai = ctrl_bam ? file("${ctrl_bam}.bai") : []

                tuple(meta, case_bam, case_bai, input_ctrl_bam, input_ctrl_bai, [])
                
                }
    
    SNIFFLES(
        sniffles_cutesv_input,                                               // Input 1: [meta, bam, bai]
        Channel.value(tuple([id: "reference"], file(reference))),               // Input 2: [meta, fasta]
        Channel.value(tuple([id: "tandem"], file(params.tandem_repeats))),      // Input 3: [meta, tandem_file]
        Channel.value(true),                                                    // Input 4: vcf_output
        Channel.value(false)                                                    // Input 5: snf_output
    )

    CUTESV(
        sniffles_cutesv_input,
        Channel.value(tuple([id: "reference"], file(reference)))
    )

    SEVERUS(
        severus_input,
        Channel.value(tuple([id: "vntr"], file(params.vntr_bed)))
    )



    // ──────────────────────────────────────────────────────────────────────
    // FIX SAMPLE NAMES in VCF HEADERS
    // ──────────────────────────────────────────────────────────────────────

    RENAME_VCF_HEADERS_SNIFFLES(
        SNIFFLES.out.vcf
    )
    RENAME_VCF_HEADERS_CUTESV(
        CUTESV.out.vcf
    )


    ch_severus_outputs = SEVERUS.out.all_vcf
        .join( SEVERUS.out.somatic_vcf, remainder: true )

    ch_severus_final_out = ch_severus_outputs
        .map { meta, all_vcf, somatic_vcf ->
            def vcf_to_use = somatic_vcf ? somatic_vcf : all_vcf

            tuple( meta, vcf_to_use )
        }

    RENAME_VCF(
        ch_severus_final_out,
        'severus'
    )
    
    RENAME_VCF_HEADERS_SEVERUS(
        RENAME_VCF.out.renamed_vcf
    )

    ch_all_caller_vcfs = RENAME_VCF_HEADERS_SNIFFLES.out.mix(
        RENAME_VCF_HEADERS_CUTESV.out, RENAME_VCF_HEADERS_SEVERUS.out
    )
    .groupTuple(by:0, size:3)

    SUMMARIZE_CALLERS(
        ch_all_caller_vcfs,
        Channel.value("raw_calls")
    )

    PLOT_RAW_CALLERS(
        SUMMARIZE_CALLERS.out.json
            .map { meta, json -> tuple([id: "raw_callers_plot"], [json]) },
        Channel.value("Raw Caller SV Counts")
    ) 

    // ──────────────────────────────────────────────────────────────────────
    // Run Jasmine to merge SVs from callers per sample
    // ──────────────────────────────────────────────────────────────────────


    ch_jasminesv_input = ch_all_caller_vcfs
        .map { meta, vcf_list ->
            tuple( meta, vcf_list, [], [])
        }

    // Prepare Jasmine input channels (per-sample)
    ch_jasmine_sample_reference = Channel.value(tuple([id: "reference"], params.reference ? file(reference) : []))
    ch_jasmine_sample_fai       = Channel.value(tuple([id: "fai"], params.reference ? file("${reference}.fai") : []))
    ch_jasmine_sample_chr_norm  = Channel.value([]) // No chr norm file

    JASMINESV_SAMPLE(
        ch_jasminesv_input,
        ch_jasmine_sample_reference,
        ch_jasmine_sample_fai,
        ch_jasmine_sample_chr_norm
    )

    emit:
        multiqc_report         = ch_multiqc_files
        versions               = ch_versions    
}

/*




    jasminesv_sample_sources = jasminesv_sample_input
        .map { meta, vcf_list, bams, sample_dists -> tuple(meta.sample ?: meta.id, vcf_list) }

    jasminesv_sample_out_keyed = JASMINESV_SAMPLE.out.vcf
        .map { meta, vcf -> tuple(meta.sample ?: meta.id, tuple(meta, vcf)) }

    jasminesv_sample_out_keyed
        .join(jasminesv_sample_sources)
        .map { sample, leftVal, src_vcfs ->
            def (meta, vcf) = leftVal
            tuple(meta, vcf, src_vcfs)
        } | JASMINE_HEADER_FIX

    ch_jasmine_sample_vcfs = JASMINE_HEADER_FIX.out.vcf
    .map { meta, vcf -> tuple(meta, vcf) }

    sample_filtered = ch_jasmine_sample_vcfs | FILTER_CHR

    sample_filtered | SORT_VCF
    sample_sorted = SORT_VCF.out.vcf

    // ──────────────────────────────────────────────────────────────────────
    // Filter SVs supported by ≥2 callers
    // ──────────────────────────────────────────────────────────────────────

    bcftools_sample_input = sample_sorted
    .map { meta, vcf ->
        def v = vcf.toString()
        def updated_meta = [id: meta.sample ?: meta.id, sample: meta.sample ?: meta.id, step: "caller_support"]
        def idx = file(v + '.csi')
        if( !idx.exists() ) idx = file(v + '.tbi')
        def idx_out = idx.exists() ? idx : []
        tuple(updated_meta, file(v), idx_out)
    }

    CALLER_SUPPORT_FILTER(
        bcftools_sample_input,
        Channel.value([]), // regions
        Channel.value([]), // targets
        Channel.value([])  // samples
    )

    // Consensus summary - simple approach
    consensus_summary_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> tuple([id: "consensus_summary"], vcf_list) }

    SUMMARIZE_CALLER_MERGED(
        consensus_summary_input,
        Channel.value("consensus")
    )

    PLOT_CONSENSUS(
        SUMMARIZE_CALLER_MERGED.out.json
            .map { meta, json -> tuple([id: "consensus_plot"], [json]) },
        Channel.value("Consensus SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // SAMPLE LEVEL AF ANNOTATION + FILTERING + ANNOTSV ANNOTATION
    // ──────────────────────────────────────────────────────────────────────

    ch_per_sample_input = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> tuple(meta, vcf) }

    ch_svdb_in_occ  = Channel.value(params.svdb_in_occ ?: [])
    ch_svdb_in_frq  = Channel.value(params.svdb_in_frq ?: [])
    ch_svdb_out_occ = Channel.value(params.svdb_out_occ ?: [])
    ch_svdb_out_frq = Channel.value(params.svdb_out_frq ?: [])
    ch_svdb_dbs     = Channel.value(params.svdb_databases ? params.svdb_databases.collect { file(it) } : [])
    ch_svdb_bedpe   = Channel.value([])

    SVDB_QUERY_SAMPLE(
        ch_per_sample_input,
        ch_svdb_in_occ,
        ch_svdb_in_frq,
        ch_svdb_out_occ,
        ch_svdb_out_frq,
        ch_svdb_dbs,
        ch_svdb_bedpe
    )

    ch_per_sample_bcftools_input = SVDB_QUERY_SAMPLE.out.vcf
    .map { meta, annotated_vcf ->
        def updated_meta = [id: meta.sample ?: meta.id, sample: meta.sample ?: meta.id, step: "af_filter"]
        def idx = file(annotated_vcf.toString() + '.csi')
        if( !idx.exists() ) idx = file(annotated_vcf.toString() + '.tbi')
        def idx_out = idx.exists() ? idx : []
        tuple(updated_meta, file(annotated_vcf), idx_out)
    }

    ch_bcftools_regions = Channel.value([])
    ch_bcftools_targets = Channel.value([])
    ch_bcftools_samples = Channel.value([])

    AF_FILTER(
        ch_per_sample_bcftools_input,
        ch_bcftools_regions,
        ch_bcftools_targets,
        ch_bcftools_samples
    )

    // Filtered summary - simple approach
    filtered_summary_input = AF_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()
        .map { vcf_list -> tuple([id: "filtered_summary"], vcf_list) }

    SUMMARIZE_CALLER_MERGED_FILTERED(
        filtered_summary_input,
        Channel.value("filtered")
    )

    PLOT_FILTERED(
        SUMMARIZE_CALLER_MERGED_FILTERED.out.json
            .map { meta, json -> tuple([id: "filtered_plot"], [json]) },
        Channel.value("Filtered SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // prepare channel for AnnotSV annotations
    // ──────────────────────────────────────────────────────────────────────

    if(!annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .map { [[id:"annotsv"], it] }
            .collect()
            .set { ch_annotsv_annotations }
    } else {
        ch_annotsv_annotations_input = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        if(annotsv_annotations.endsWith(".tar.gz")){
            UNTAR_ANNOTSV(ch_annotsv_annotations_input)
            UNTAR_ANNOTSV.out.untar
                .collect()
                .set { ch_annotsv_annotations }
        } else {
            ch_annotsv_annotations = Channel.fromPath(annotsv_annotations).map{[[id:"annotsv_annotations"], it]}.collect()
        }
    }

    ch_candidate_genes      = Channel.value(tuple([id: "candidate_genes"], []))
    ch_false_positive_snv   = Channel.value(tuple([id: "false_positive_snv"], []))
    ch_gene_transcripts     = Channel.value(tuple([id: "gene_transcripts"], []))

    ANNOTSV_PER_SAMPLE_RAW(
        SVDB_QUERY_SAMPLE.out.vcf
            .map { meta, vcf -> 
                def updated_meta = meta.clone()
                updated_meta.id = "${meta.sample ?: meta.id}_annotated"
                tuple(updated_meta, vcf, [], []) 
            },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    ANNOTSV_PER_SAMPLE(
        AF_FILTER.out.vcf
            .map { meta, vcf ->
                def updated_meta = [
                    id: "${meta.sample ?: meta.id}_filtered",
                    sample: meta.sample ?: meta.id, 
                    step: "final_annotation"
                ]
                tuple(updated_meta, vcf, [], [])
            },
        ch_annotsv_annotations,
        ch_candidate_genes,
        ch_false_positive_snv,
        ch_gene_transcripts
    )

    // ──────────────────────────────────────────────────────────────────────
    // Continue to cohort-level analyses
    // ──────────────────────────────────────────────────────────────────────

    sample_consensus_vcfs = CALLER_SUPPORT_FILTER.out.vcf
        .map { meta, vcf -> vcf }
        .collect()

    jasminesv_cohort_input = sample_consensus_vcfs
        .map { vcf_list ->
        tuple([id: "cohort"], vcf_list, [], [])
    }

    ch_jasmine_cohort_reference = Channel.value(tuple([id: "reference"], params.reference ? file(reference) : []))
    ch_jasmine_cohort_fai       = Channel.value(tuple([id: "fai"], params.reference ? file("${reference}.fai") : []))
    ch_jasmine_cohort_chr_norm  = Channel.value([]) // No chr norm file

    JASMINESV_COHORT(
        jasminesv_cohort_input,
        ch_jasmine_cohort_reference,
        ch_jasmine_cohort_fai,
        ch_jasmine_cohort_chr_norm
    )

    // ──────────────────────────────────────────────────────────────────────
    // SV annotation using SVDB (cohort-level)
    // ──────────────────────────────────────────────────────────────────────

    svdb_cohort_input = JASMINESV_COHORT.out.vcf
        .map { meta, cohort_vcf ->
        tuple(meta, cohort_vcf)
    }

    ch_svdb_cohort_in_occ  = Channel.value(params.svdb_in_occ ?: [])
    ch_svdb_cohort_in_frq  = Channel.value(params.svdb_in_frq ?: [])
    ch_svdb_cohort_out_occ = Channel.value(params.svdb_out_occ ?: [])
    ch_svdb_cohort_out_frq = Channel.value(params.svdb_out_frq ?: [])
    ch_svdb_cohort_dbs     = Channel.value(params.svdb_databases ? params.svdb_databases.collect { file(it) } : [])
    ch_svdb_cohort_bedpe   = Channel.value([])

    SVDB_QUERY_COHORT(
        svdb_cohort_input,
        ch_svdb_cohort_in_occ,
        ch_svdb_cohort_in_frq,
        ch_svdb_cohort_out_occ,
        ch_svdb_cohort_out_frq,
        ch_svdb_cohort_dbs,
        ch_svdb_cohort_bedpe
    )

    ch_candidate_genes_cohort    = Channel.value(tuple([id: "candidate_genes"], []))
    ch_false_positive_snv_cohort = Channel.value(tuple([id: "false_positive_snv"], []))
    ch_gene_transcripts_cohort   = Channel.value(tuple([id: "gene_transcripts"], []))

    ANNOTSV_COHORT_RAW(
        SVDB_QUERY_COHORT.out.vcf
            .map { meta, vcf -> 
                def updated_meta = [
                    id: "${meta.id}_annotated",
                    sample: "${meta.id}_annotated",
                ]
                tuple(updated_meta, vcf, [], [])
            },
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
    )

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_ANNOTATED(
        SVDB_QUERY_COHORT.out.vcf
            .map { meta, vcf -> tuple([id: "cohort_annotated_summary"], vcf) },
        Channel.value("cohort_annotated")
    )

    PLOT_COHORT_ANNOTATED(
        SUMMARIZE_COHORT_ANNOTATED.out.json
            .map { meta, json -> tuple([id: "cohort_annotated_plot"], [json]) },
        Channel.value("Cohort Annotated SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // Filter annotated SVs based on AF
    // ──────────────────────────────────────────────────────────────────────

    bcftools_cohort_input = SVDB_QUERY_COHORT.out.vcf.map { meta, annotated_vcf ->
        tuple(meta, annotated_vcf, [])
    }

    AF_FILTER_COHORT(
        bcftools_cohort_input,
        ch_bcftools_regions,
        ch_bcftools_targets,
        ch_bcftools_samples
    )

    // Cohort summaries - one per cohort directory
    SUMMARIZE_COHORT_FILTERED(
        AF_FILTER_COHORT.out.vcf
            .map { meta, vcf -> tuple([id: "cohort_filtered_summary"], vcf) },
        Channel.value("cohort_filtered")
    )

    PLOT_COHORT_FILTERED(
        SUMMARIZE_COHORT_FILTERED.out.json
            .map { meta, json -> tuple([id: "cohort_filtered_plot"], [json]) },
        Channel.value("Cohort Filtered SV Counts")
    )

    // ──────────────────────────────────────────────────────────────────────
    // Structural variant annotation using AnnotSV
    // ──────────────────────────────────────────────────────────────────────

    annotsv_input = AF_FILTER_COHORT.out.vcf.map { meta, filtered_vcf ->
        def updated_meta = [
            id: 'cohort_filtered',
            sample: 'cohort_filtered'
        ]
        tuple(updated_meta, filtered_vcf, [], [])
    }

    ANNOTSV_COHORT(
        annotsv_input,
        ch_annotsv_annotations,
        ch_candidate_genes_cohort,
        ch_false_positive_snv_cohort,
        ch_gene_transcripts_cohort
)

    // ──────────────────────────────────────────────────────────────────────
    // Collate and save software versions
    // ──────────────────────────────────────────────────────────────────────

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SNIFFLES.out.versions)
    ch_versions = ch_versions.mix(CUTESV.out.versions)
    ch_versions = ch_versions.mix(SEVERUS_WITH_CONTROL.out.versions)
    ch_versions = ch_versions.mix(SEVERUS_NO_CONTROL.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(JASMINESV_COHORT.out.versions)
    ch_versions = ch_versions.mix(SORT_VCF.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(SVDB_QUERY_COHORT.out.versions)
    ch_versions = ch_versions.mix(CALLER_SUPPORT_FILTER.out.versions)
    ch_versions = ch_versions.mix(AF_FILTER.out.versions)
    ch_versions = ch_versions.mix(AF_FILTER_COHORT.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_PER_SAMPLE_RAW.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_PER_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(ANNOTSV_COHORT.out.versions)

    // Handle conditional AnnotSV installation versions
    if(!annotsv_annotations) {
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'ontvar_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
        multiqc_report         = MULTIQC.out.report.toList()
        versions               = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
