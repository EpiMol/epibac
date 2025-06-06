rule epibac_summary:
    input:
        counts=lambda wc: [
            f"{OUTDIR}/qc/count_reads/{sample}_counts.txt"
            for sample in get_samples().index
        ],
        amrfinder=lambda wc: [
            f"{OUTDIR}/amr_mlst/{sample}_amrfinder.tsv"
            for sample in get_samples()[get_samples()["illumina_r2"].notnull()].index
        ],
        mlst=lambda wc: [
            f"{OUTDIR}/amr_mlst/{sample}_mlst.tsv"
            for sample in get_samples()[get_samples()["illumina_r2"].notnull()].index
        ],
        resfinder=lambda wc: [
            f"{OUTDIR}/amr_mlst/resfinder/{sample}/ResFinder_results.txt"
            for sample in get_samples()[get_samples()["illumina_r2"].notnull()].index
        ],
    output:
        report_dir=directory(f"{OUTDIR}/report"),
        tsv=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.tsv",
        xlsx=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.xlsx",
    params:
        input=directory(f"{OUTDIR}/amr_mlst"),
    log:
        f"{LOGDIR}/report/summary.log",
    conda:
        "../envs/epibac_report.yml"
    threads: get_resource("summary", "threads")
    resources:
        mem_mb=get_resource("summary", "mem"),
        walltime=get_resource("summary", "walltime"),
    script:
        "../scripts/epibac_summary.py"


rule epibac_summary_gestlab:
    input:
        validated_samples_info=f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv",
        results=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.tsv",
        plasmids=f"{OUTDIR}/mge_analysis/unified_plasmid_analysis/{TAG_RUN}_plasmid_analysis_summary.tsv"
    output:
        gestlab_report_csv=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC_GESTLAB.csv",
        gestlab_report_xlsx=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC_GESTLAB.xlsx"
    log:
        f"{LOGDIR}/report/summary_gestlab.log",
    conda:
        "../envs/epibac_report.yml"
    container:
        "docker://alesanzdro/epibac_report:1.0"
    threads: get_resource("summary", "threads")
    resources:
        mem_mb=get_resource("summary", "mem"),
        walltime=get_resource("summary", "walltime"),
    shell:
        """
        python {workflow.basedir}/scripts/epibac_summary_gestlab.py \
            --samplesinfo {input.validated_samples_info} \
            --results {input.results} \
            --plasmids {input.plasmids} \
            --output_csv {output.gestlab_report_csv} \
            --output_xlsx {output.gestlab_report_xlsx} \
            --tag_run {TAG_RUN} \
            --log_file {log} \
            2> {log}
        """

rule copy_sequencing_files:
    input:
        gestlab_report=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC_GESTLAB.csv",
        report_tsv=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.tsv",
        report_xlsx=f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.xlsx",
    output:
        copy_log=f"{OUTDIR}/report/{TAG_RUN}_file_copy_log.txt",
    log:
        f"{LOGDIR}/report/copy_sequencing_files.log",
    params:
        config_file=lambda _, input: workflow.configfiles[0] if workflow.configfiles else "config.yaml"
    conda:
        "../envs/epibac_report.yml"
    threads: 1
    resources:
        mem_mb=2000,
        walltime=480,
    shell:
        """
        python {workflow.basedir}/scripts/copy_gva_files.py \
            {input.gestlab_report} \
            {input.report_tsv} \
            {input.report_xlsx} \
            {output.copy_log} \
            {OUTDIR} \
            {TAG_RUN} \
            --config-file {params.config_file} \
            2> {log}
        """