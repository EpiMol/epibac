rule epibac_assembly:
    input:
        r1=rules.epibac_fastp_pe.output.r1,
        r2=rules.epibac_fastp_pe.output.r2,
        nanopore=lambda wc: get_nanopore_fastq(wc) if has_nanopore_data(wc) else [],
    output:
        fasta="{}/assembly/{{sample}}/{{sample}}.fasta".format(OUTDIR),
        gfa="{}/assembly/{{sample}}/{{sample}}.gfa".format(OUTDIR),
        log_file="{}/assembly/{{sample}}/{{sample}}.log".format(OUTDIR),
    log:
        f"{LOGDIR}/unicycler/{{sample}}.log",
    conda:
        "../envs/epibac_assembly.yml"
    container:
        "docker://alesanzdro/epibac_assembly:1.0"
    threads: get_resource("unicycler", "threads")
    resources:
        mem_mb=get_resource("unicycler", "mem"),
        walltime=get_resource("unicycler", "walltime"),
    params:
        output_dir="{}/assembly/{{sample}}".format(OUTDIR),
        nanopore_flag=lambda wc: "-l " + get_nanopore_fastq(wc) if has_nanopore_data(wc) else "",
        dorado_model=lambda wc: get_dorado_model(wc) if has_nanopore_data(wc) else "",
    shell:
        r"""
        set -e
        
        # Log information about assembly type
        if [ -n "{params.nanopore_flag}" ]; then
            echo "Running hybrid assembly (Illumina + Nanopore) for sample {wildcards.sample}" >> {log}
            echo "Nanopore file: {input.nanopore}" >> {log}
            echo "Dorado model used: {params.dorado_model}" >> {log}
        else
            echo "Running Illumina-only assembly for sample {wildcards.sample}" >> {log}
        fi
        
        unicycler \
            -t {threads} \
            {config[params][unicycler][extra]} \
            -1 {input.r1} \
            -2 {input.r2} \
            {params.nanopore_flag} \
            -o {params.output_dir} \
            &> {log} || {{
                echo "Unicycler falló, probablemente debido a una entrada pequeña o artefactos. Creando archivos de salida vacíos." >> {log}
                touch {output.fasta} {output.gfa} {output.log_file}
                exit 0
            }}

        mv {params.output_dir}/assembly.fasta {output.fasta}
        mv {params.output_dir}/assembly.gfa {output.gfa}
        mv {params.output_dir}/unicycler.log {output.log_file}
        """

rule epibac_dorado_basecall:
    """
    Regla opcional para realizar basecalling con Dorado desde datos raw de nanopore.
    Esta regla se ejecuta solo si se proporciona un directorio de datos raw en lugar de FASTQ.
    """
    input:
        pod5_dir=lambda wc: config.get("params", {}).get("nanopore", {}).get("pod5_dir", ""),
    output:
        fastq="{}/nanopore/{{sample}}/{{sample}}_basecalled.fastq.gz".format(OUTDIR),
    log:
        f"{LOGDIR}/dorado/{{sample}}.log",
    conda:
        "../envs/epibac_nanopore.yml"  # Necesitaría crear este environment
    threads: get_resource("dorado", "threads", default=8)
    resources:
        mem_mb=get_resource("dorado", "mem", default=8000),
        walltime=get_resource("dorado", "walltime", default=480),
        gpu=get_resource("dorado", "gpu", default=1),  # Para uso de GPU si está disponible
    params:
        model=lambda wc: get_dorado_model(wc),
        sample_pod5_dir=lambda wc: f"{config.get('params', {}).get('nanopore', {}).get('pod5_dir', '')}/{wc.sample}",
    shell:
        r"""
        # Verificar que el directorio de POD5 existe
        if [ ! -d "{params.sample_pod5_dir}" ]; then
            echo "ERROR: No se encuentra el directorio POD5 para la muestra {wildcards.sample}: {params.sample_pod5_dir}" >> {log}
            echo "Esta regla es opcional y se salta si no hay datos raw disponibles." >> {log}
            touch {output.fastq}
            exit 0
        fi
        
        # Crear directorio de salida
        mkdir -p $(dirname {output.fastq})
        
        echo "Iniciando basecalling con Dorado para muestra {wildcards.sample}" >> {log}
        echo "Modelo: {params.model}" >> {log}
        echo "Directorio POD5: {params.sample_pod5_dir}" >> {log}
        
        # Ejecutar Dorado basecalling
        dorado basecaller \
            {params.model} \
            {params.sample_pod5_dir} \
            --verbose \
            --device cuda \
            2>> {log} | gzip > {output.fastq}
        
        echo "Basecalling completado para muestra {wildcards.sample}" >> {log}
        echo "Archivo generado: {output.fastq}" >> {log}
        """
