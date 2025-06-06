# ===== ANÁLISIS DE PLÁSMIDOS =====
rule epibac_platon:
    """Análisis de plásmidos con Platon"""
    input:
        fasta=f"{OUTDIR}/assembly/{{sample}}/{{sample}}.fasta",
        platon_db="resources/databases/platon/db",
        flag=PLATON_DB_FLAG
    output:
        tsv=f"{OUTDIR}/mge_analysis/plasmids/platon/{{sample}}/{{sample}}.tsv",
        dir=directory(f"{OUTDIR}/mge_analysis/plasmids/platon/{{sample}}")
    log:
        f"{LOGDIR}/platon/{{sample}}.log"
    conda:
        "../envs/epibac_mge.yml"
    threads: get_resource("platon", "threads")
    resources:
        mem_mb=get_resource("platon", "mem"),
        walltime=get_resource("platon", "walltime"),
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "[ERROR] El archivo FASTA {input.fasta} está vacío" &> {log}
            mkdir -p {output.dir}
            touch {output.tsv}
            exit 0
        fi

        platon {input.fasta} \
          --output {output.dir} \
          --db {input.platon_db} \
          --threads {threads} \
          --characterize \
          --mode accuracy &> {log}
        
        # Crear archivo TSV de resumen si no existe
        if [ ! -f {output.tsv} ] && [ -f {output.dir}/{wildcards.sample}.tsv ]; then
            cp {output.dir}/{wildcards.sample}.tsv {output.tsv}
        elif [ ! -f {output.tsv} ]; then
            touch {output.tsv}
        fi
        """

rule epibac_mob_suite:
    """Análisis de movilización con MOB-suite"""
    input:
        fasta=f"{OUTDIR}/assembly/{{sample}}/{{sample}}.fasta"
    output:
        results=f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{{sample}}/mobtyper_results.txt",
        dir=directory(f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{{sample}}")
    log:
        f"{LOGDIR}/mob_suite/{{sample}}.log"
    conda:
        "../envs/epibac_mge_mob_suite.yml"
    threads: get_resource("mob_suite", "threads")
    resources:
        mem_mb=get_resource("mob_suite", "mem"),
        walltime=get_resource("mob_suite", "walltime"),
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "[ERROR] El archivo FASTA {input.fasta} está vacío" &> {log}
            mkdir -p {output.dir}
            touch {output.results}
            exit 0
        fi

        # Ejecutar MOB-recon
        mob_recon \
          -i {input.fasta} \
          -o {output.dir} \
          -n {threads} \
          -s {wildcards.sample} \
          --min_length 500 \
          --max_plasmid_size 2000000 \
          --min_rep_ident 80 \
          --min_mob_ident 80 \
          --min_con_ident 80 \
          -u -c -f &> {log}
        
        # Ejecutar mob_typer en plásmidos identificados
        if ls {output.dir}/plasmid_*.fasta 1> /dev/null 2>&1; then
            echo "→ Tipificando plásmidos con MOB-typer..." >> {log}
            for PLASMID_FILE in {output.dir}/plasmid_*.fasta; do
                PLASMID_ID=$(basename "$PLASMID_FILE" .fasta)
                mob_typer \
                  --infile "$PLASMID_FILE" \
                  --out_file "{output.dir}/${{PLASMID_ID}}_mob_typer.txt" \
                  --num_threads {threads} >> {log} 2>&1
            done
            
            # Consolidar resultados
            cat {output.dir}/*_mob_typer.txt > {output.results} 2>/dev/null || touch {output.results}
        else
            echo "No se identificaron plásmidos" >> {log}
            touch {output.results}
        fi
        """

# ===== ANÁLISIS DE TRANSPOSONES =====
rule epibac_te:
    """Análisis de transposones con abricate"""
    input:
        fasta=f"{OUTDIR}/assembly/{{sample}}/{{sample}}.fasta",
        flag=TRANSPOSONS_DB_FLAG
    output:
        tsv=f"{OUTDIR}/mge_analysis/transposons/{{sample}}_abricate.tsv"
    log:
        f"{LOGDIR}/transposons/{{sample}}.log"
    conda:
        "../envs/epibac_mge.yml"
    threads: get_resource("transposons", "threads", 4)
    resources:
        mem_mb=get_resource("transposons", "mem", 4000),
        walltime=get_resource("transposons", "walltime", "01:00:00")
    shell:
        """
        if [ ! -s {input.fasta} ]; then
            echo "[ERROR] El archivo FASTA {input.fasta} está vacío" &> {log}
            touch {output.tsv}
            exit 0
        fi

        abricate \
          --threads {threads} \
          --db transposons \
          --minid 90 \
          --mincov 80 \
          {input.fasta} > {output.tsv} 2> {log}
        """


#rule epibac_mge_analysis:
#    """Análisis completo de elementos móviles genéticos"""
#    input:
#        # Plásmidos
#        platon=f"{OUTDIR}/mge_analysis/plasmids/platon/{{sample}}/{{sample}}.tsv",
#        mob_suite=f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{{sample}}/mobtyper_results.txt",
#        # Transposones
#        transposons=f"{OUTDIR}/mge_analysis/transposons/{{sample}}_abricate.tsv",
#        # Integrones
#        #integrons=f"{OUTDIR}/mge_analysis/integrons/{{sample}}/Results_Integron_Finder_{{sample}}/{{sample}}.integrons",
#        # Secuencias de inserción
#        #is_sequences=f"{OUTDIR}/mge_analysis/insertion_sequences/{{sample}}/{{sample}}/{{sample}}.fasta.sum",
#        # Profagos
#        #phages=f"{OUTDIR}/mge_analysis/phages/{{sample}}/phage.fasta",
#        #phage_vfdb=f"{OUTDIR}/mge_analysis/phages/{{sample}}/phage_abricate_vfdb.tsv",
#        #phage_resfinder=f"{OUTDIR}/mge_analysis/phages/{{sample}}/phage_abricate_resfinder.tsv"
#    output:
#        summary=f"{OUTDIR}/mge_analysis/{{sample}}_mge_summary.tsv"
#    conda:
#        "../envs/epibac_report.yml"
#    threads: get_resource("summary", "threads")
#    resources:
#        mem_mb=get_resource("summary", "mem"),
#        walltime=get_resource("summary", "walltime"),    
#    log:
#        f"{LOGDIR}/mge_summary/{{sample}}.log"
#    script:
#        "../scripts/summarize_mge_analysis.py"


# ===== ANÁLISIS UNIFICADO DE PLÁSMIDOS =====

def get_plasmid_inputs():
    """
    Obtiene los inputs de plásmidos de manera segura.
    Devuelve listas vacías si no hay muestras disponibles.
    """
    try:
        sample_ids = get_sample_index_if_exists()
        if sample_ids:
            platon_files = [
                f"{OUTDIR}/mge_analysis/plasmids/platon/{sample}/{sample}.tsv"
                for sample in sample_ids
            ]
            mob_suite_files = [
                f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{sample}/mobtyper_results.txt"
                for sample in sample_ids
            ]
            return {"platon": platon_files, "mob_suite": mob_suite_files}
        else:
            return {"platon": [], "mob_suite": []}
    except:
        return {"platon": [], "mob_suite": []}

rule epibac_unified_plasmid_analysis:
    """Análisis unificado de plásmidos combinando Platon y MOB-suite para todas las muestras"""
    input:
        # Dependencias explícitas de todos los archivos de Platon y MOB-suite
        platon_files=lambda wc: expand(
            f"{OUTDIR}/mge_analysis/plasmids/platon/{{sample}}/{{sample}}.tsv",
            sample=get_sample_index_if_exists()
        ) if get_sample_index_if_exists() else [],
        mob_suite_files=lambda wc: expand(
            f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{{sample}}/mobtyper_results.txt",
            sample=get_sample_index_if_exists()
        ) if get_sample_index_if_exists() else [],
        # También depender del archivo de validación de muestras
        validation_file=f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv"
    output:
        tsv=f"{OUTDIR}/mge_analysis/unified_plasmid_analysis/{TAG_RUN}_plasmid_analysis_summary.tsv",
        excel=f"{OUTDIR}/mge_analysis/unified_plasmid_analysis/{TAG_RUN}_Reporte_Plasmidos_Simplificado.xlsx",
        output_dir=directory(f"{OUTDIR}/mge_analysis/unified_plasmid_analysis")
    params:
        platon_dir=f"{OUTDIR}/mge_analysis/plasmids/platon",
        mob_suite_dir=f"{OUTDIR}/mge_analysis/plasmids/mob_suite",
        output_prefix=f"{TAG_RUN}",
        validation_file=f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv"
    log:
        f"{LOGDIR}/unified_plasmid_analysis/unified_analysis.log"
    conda:
        "../envs/epibac_report.yml"
    threads: get_resource("summary", "threads", 1)
    resources:
        mem_mb=get_resource("summary", "mem", 4000),
        walltime=get_resource("summary", "walltime", "02:00:00")
    shell:
        """
        # Crear directorios de salida y logs
        mkdir -p {output.output_dir}
        mkdir -p $(dirname {log})
        
        echo "Iniciando análisis unificado de plásmidos" > {log}
        echo "Fecha: $(date)" >> {log}
        echo "TAG_RUN: {params.output_prefix}" >> {log}
        
        # Verificar si existe el archivo de validación de muestras
        if [ -f {params.validation_file} ]; then
            echo "Archivo de validación encontrado: {params.validation_file}" >> {log}
        else
            echo "ADVERTENCIA: Archivo de validación no encontrado: {params.validation_file}" >> {log}
            echo "El análisis continuará buscando archivos existentes en los directorios" >> {log}
        fi
        
        # Verificar que existan directorios de entrada
        if [ ! -d {params.platon_dir} ]; then
            echo "Creando directorio Platon: {params.platon_dir}" >> {log}
            mkdir -p {params.platon_dir}
        fi
        
        if [ ! -d {params.mob_suite_dir} ]; then
            echo "Creando directorio MOB-suite: {params.mob_suite_dir}" >> {log}
            mkdir -p {params.mob_suite_dir}
        fi
        
        # Ejecutar el script unificado
        echo "Ejecutando script unificado..." >> {log}
        python {workflow.basedir}/scripts/unified_plasmid_analysis.py \
            --platon_path {params.platon_dir} \
            --mobsuite_path {params.mob_suite_dir} \
            --output_dir {output.output_dir} \
            --min_length 2000 \
            --log_file {log} 2>&1 | tee -a {log}
        
        # Verificar que se generaron los archivos de salida y renombrarlos
        base_tsv="{output.output_dir}/plasmid_analysis_summary.tsv"
        base_excel="{output.output_dir}/Reporte_Plasmidos_Simplificado.xlsx"
        
        if [ -f "$base_tsv" ]; then
            mv "$base_tsv" {output.tsv}
            echo "Archivo TSV renombrado a: {output.tsv}" >> {log}
        else
            echo "ADVERTENCIA: No se generó el archivo TSV de resumen, creando archivo vacío" >> {log}
            echo "sample_id	plasmid_id	plasmid_length	tool_source	contig_type	notes" > {output.tsv}
        fi
        
        if [ -f "$base_excel" ]; then
            mv "$base_excel" {output.excel}
            echo "Archivo Excel renombrado a: {output.excel}" >> {log}
        else
            echo "ADVERTENCIA: No se generó el archivo Excel de resumen" >> {log}
            # Crear un archivo Excel básico usando Python
            python -c "
import pandas as pd
import sys
try:
    df = pd.DataFrame({{'sample_id': [], 'plasmid_id': [], 'plasmid_length': [], 'tool_source': [], 'contig_type': [], 'notes': []}})
    df.to_excel('{output.excel}', index=False)
    print('Archivo Excel vacío creado exitosamente')
except Exception as e:
    print(f'Error creando archivo Excel: {{e}}')
    with open('{output.excel}', 'w') as f:
        f.write('')
" >> {log} 2>&1
        fi
        
        echo "Análisis unificado de plásmidos completado" >> {log}
        """

