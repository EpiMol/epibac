rule epibac_amr:
    input:
        setup_db=AMRFINDER_DB_FLAG,
        fasta=f"{OUTDIR}/assembly/{{sample}}/{{sample}}.fasta",
        prokka=lambda wc: f"{OUTDIR}/annotation/{wc.sample}/{wc.sample}.faa",
        gff=lambda wc: f"{OUTDIR}/annotation/{wc.sample}/{wc.sample}.gff",
    output:
        gff=f"{OUTDIR}/annotation/{{sample}}/{{sample}}_amrfinder.gff",
        tsv=f"{OUTDIR}/amr_mlst/{{sample}}_amrfinder.tsv",
    log:
        f"{LOGDIR}/amrfinder/{{sample}}.log",
    threads: get_resource("amrfinder", "threads")
    resources:
        mem_mb=get_resource("amrfinder", "mem"),
        walltime=get_resource("amrfinder", "walltime"),
    params:
        name=lambda wc: f"{wc.sample}",
        db_dir=f"{AMRFINDER_DB_DIR}/latest",
    conda:
        "../envs/epibac_amr_annotation_plus.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation_plus:1.0"
    shell:
        r"""
        if [ ! -s {input.fasta} ]; then
            echo "[ERROR] El archivo FASTA {input.fasta} está vacío" &> {log}
            touch {output.gff}
            touch {output.tsv}
            exit 0
        fi

        # Prepara GFF para AMRFinder
        perl -pe '/^##FASTA/ && exit; s/(\W)Name=/$1OldName=/i; s/ID=([^;]+)/ID=$1;Name=$1/' \
          {input.gff} > {output.gff}

        # Ejecuta AMRFinder con ruta explícita de base de datos
        amrfinder \
          --plus \
          --threads {threads} \
          --name {params.name} \
          --database {params.db_dir} \
          -n {input.fasta} \
          -p {input.prokka} \
          -g {output.gff} \
          --coverage_min 0.7 \
          > {output.tsv} 2> {log}
        """


rule epibac_resfinder:
    input:
        setup_db=RESFINDER_DB_FLAG,
        fasta=f"{OUTDIR}/assembly/{{sample}}/{{sample}}.fasta",
    output:
        dir=directory("{}/amr_mlst/resfinder/{{sample}}".format(OUTDIR)),
        res="{}/amr_mlst/resfinder/{{sample}}/ResFinder_results.txt".format(OUTDIR),
    log:
        f"{LOGDIR}/resfinder/{{sample}}.log",
    conda:
        "../envs/epibac_amr_annotation_plus.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation_plus:1.0"
    threads: get_resource("resfinder", "threads")
    resources:
        mem_mb=get_resource("resfinder", "mem"),
        walltime=get_resource("resfinder", "walltime"),
    params:
        name=lambda wc: f"{wc.sample}",
        db_path=RESFINDER_DB_DIR,
        extra=lambda wc: config["params"]["resfinder"].get("extra", ""),
    shell:
        r"""
        # Verifica si el archivo FASTA es vacío o no
        if [ ! -s {input.fasta} ]; then
            echo "[ERROR] El archivo FASTA {input.fasta} está vacío" &> {log}
            mkdir -p {output.dir}
            touch {output.res}
            exit 0
        fi

        run_resfinder.py \
            -db_res {params.db_path} \
            -o {output.dir} \
            {params.extra} \
            -ifa {input.fasta} \
            &> {log}
        """


rule epibac_mlst:
    input:
        assembly = lambda wc: f"{OUTDIR}/assembly/{wc.sample}/{wc.sample}.fasta",
    output:
        tsv="{}/amr_mlst/{{sample}}_mlst.tsv".format(OUTDIR),
    log:
        f"{LOGDIR}/mlst/{{sample}}.log",
    conda:
        "../envs/epibac_amr_annotation.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation:1.0"
    threads: get_resource("mlst", "threads")
    resources:
        mem_mb=get_resource("mlst", "mem"),
        walltime=get_resource("mlst", "walltime"),
    params:
        name=lambda wc: f"{wc.sample}",
        blastdb_path=f"{MLST_DB_DIR}/blast/mlst.fa",
        datadir_path=MLST_DB_DIR,
    shell:
        """
        if [[ -n "${{CONDA_PREFIX:-}}" ]]; then
            perl_path=$(which perl)
            echo "Usando perl de conda: $perl_path" >&2
            export PATH=${{CONDA_PREFIX}}/bin:$PATH
        else
            echo "Entorno Singularity/Docker (sin conda)" >&2
        fi
    
        # Verifica si el archivo FASTA es vacío o no
        if [[ ! -s {input.assembly} ]]; then
            echo "[ERROR] El archivo FASTA {input.assembly} está vacío" &> {log}
            touch {output.tsv}
            exit 0
        fi

        # Ejecutar MLST (igual para ambos entornos)
        mlst \
        --exclude "ecoli" \
        --label {params.name} \
        --blastdb {params.blastdb_path} \
        --datadir {params.datadir_path} \
        {input.assembly} \
        > {output.tsv} \
        2> {log}
        """
