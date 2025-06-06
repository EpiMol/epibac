rule setup_kraken2_database:
    output:
        flag=KRAKEN_DB_FLAG,
    log:
        KRAKEN_DB_LOG,
    conda:
        "../envs/epibac_qc.yml"
    container:
        "docker://alesanzdro/epibac_qc:1.0"
    params:
        db_url=KRAKEN_DB_URL,
        db_name=KRAKEN_DB_NAME,
        db_dir=KRAKEN_DB_DIR,
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p $(dirname {log})

        if [ ! -f "{params.db_dir}/taxo.k2d" ]; then
            echo "[INFO] Descargando base de datos Kraken2 desde {params.db_url}" &>> {log}
            wget --no-check-certificate -O {params.db_dir}/{params.db_name}.tar.gz {params.db_url} &>> {log}

            echo "[INFO] Descomprimiendo base de datos..." &>> {log}
            tar -xvzf {params.db_dir}/{params.db_name}.tar.gz -C {params.db_dir} &>> {log}

            echo "[INFO] Eliminando archivo .tar.gz..." &>> {log}
            rm -f {params.db_dir}/{params.db_name}.tar.gz
        else
            echo "[INFO] La base de datos ya está instalada, se omite la descarga." &>> {log}
        fi

        touch {output.flag}
        """


rule setup_prokka_database:
    """Configura todas las bases de datos de Prokka (HMM, BLAST, CM) en un directorio externo."""
    output:
        flag=PROKKA_DB_FLAG,
    log:
        PROKKA_DB_LOG,
    params:
        db_dir=PROKKA_DB_DIR,
        cm_files=["Archaea", "Bacteria", "Viruses"],
        skip=should_skip("prokka"),
    conda:
        "../envs/epibac_amr_annotation.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation:1.0"
    threads: 4
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "Omitiendo configuración de Prokka (skip_prokka=true)" > {log}
            touch {output.flag}
        else
            bash {workflow.basedir}/scripts/setup_prokka.sh {params.db_dir} {output.flag} {log}
        fi
        """


rule setup_amrfinder_database:
    output:
        flag=AMRFINDER_DB_FLAG,
        version=AMRFINDER_DB_DIR + "/VERSION.txt",
    log:
        AMRFINDER_DB_LOG,
    conda:
        "../envs/epibac_amr_annotation_plus.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation_plus:1.0"
    params:
        db_dir=AMRFINDER_DB_DIR,
    shell:
        r"""
        echo "[INFO] Descargando base de datos AMRFinder en {params.db_dir}" &>> {log}
        mkdir -p {params.db_dir}

        amrfinder_update --force_update --database {params.db_dir} &>> {log}

        # Detectar la versión descargada (directorio más reciente)
        VERSION=$(ls -1 {params.db_dir} | grep -E '^[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}}\..*' | sort -r | head -n 1)

        # Crear enlace simbólico "latest" -> versión
        ln -sfn "$VERSION" {params.db_dir}/latest

        # Guardar metadata mínima
        echo "Fecha de instalación: $(date --iso-8601=seconds)" > {output.version}
        echo "Versión: $VERSION" >> {output.version}
        echo "Ruta real: {params.db_dir}/$VERSION" >> {output.version}

        touch {output.flag}
        """


rule setup_resfinder_database:
    output:
        flag=RESFINDER_DB_FLAG,
    log:
        RESFINDER_DB_LOG,
    conda:
        "../envs/epibac_amr_annotation_plus.yml"
    container:
        "docker://alesanzdro/epibac_amr_annotation_plus:1.0"
    params:
        db_dir=RESFINDER_DB_DIR,
    shell:
        """
        mkdir -p {params.db_dir}

        if [ ! -d "{params.db_dir}/resfinder_db" ]; then
            echo "[INFO] Clonando base de datos de ResFinder..." &>> {log}
            git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git {params.db_dir} &>> {log}
        else
            echo "[INFO] Base de datos ResFinder ya existe." &>> {log}
        fi

        touch {output.flag}
        """


rule metadata_prokka_database:
    input:
        flag=PROKKA_DB_FLAG,
    output:
        version=PROKKA_DB_DIR + "/VERSION.txt",
        md5=PROKKA_DB_DIR + "/checksums.md5",
        readme=PROKKA_DB_DIR + "/README.txt",
    params:
        db_dir=PROKKA_DB_DIR,
        enable=config.get("enable_metadata", False),
    run:
        if not params.enable:
            print("Skipping metadata for Prokka DB (config[enable_metadata]=false)")
            for o in output:
                shell(f"touch {o}")
            return

            # Escribir versión detectada de HMM (si existe)
        version_info = "unknown"
        hmm_path = os.path.join(params.db_dir, "PGAP.hmm")
        if os.path.exists(hmm_path):
            version_info = shell(
                f"grep -m1 '^ACC' {hmm_path} | cut -f2", read=True
            ).strip()
        with open(output.version, "w") as f:
            f.write(f"PGAP Version: {version_info}\n")

            # Calcular checksums
        shell(
            f"cd {params.db_dir} && find . -type f -exec md5sum {{}} + > {output.md5}"
        )

        # README con contexto
        with open(output.readme, "w") as f:
            f.write(
                f"""# Prokka Custom Database Setup
Fecha de instalación: {datetime.datetime.now().isoformat()}
Ruta: {params.db_dir}
Versión HMM: {version_info}
Regla base: setup_prokka_database
Entorno: Docker/Singularity compatible
"""
            )

rule setup_transposons_database:
    """Configura la base de datos de transposones para abricate"""
    input:
        transposons_gz=f"{TRANSPOSONS_DB_DIR}/all_te.fasta.gz"
    output:
        transposons_fasta=f"{TRANSPOSONS_DB_DIR}/all_te.fasta",
        flag=TRANSPOSONS_DB_FLAG,
    log:
        TRANSPOSONS_DB_LOG,
    conda:
        "../envs/epibac_mge.yml"
    shell:
        """
        # Crear directorio para logs si no existe
        mkdir -p $(dirname {log})
        
        echo "Configurando base de datos de transposones para abricate..." &> {log}
        
        # Verificar que el archivo comprimido existe
        if [ ! -f {input.transposons_gz} ]; then
            echo "[ERROR] No se encuentra el archivo comprimido: {input.transposons_gz}" >> {log}
            exit 1
        fi
        
        # Descomprimir el archivo si no existe la versión descomprimida
        if [ ! -f {output.transposons_fasta} ]; then
            echo "Descomprimiendo {input.transposons_gz}..." >> {log}
            gunzip -c {input.transposons_gz} > {output.transposons_fasta}
            
            # Verificar que la descompresión fue correcta
            if [ ! -f {output.transposons_fasta} ] || [ ! -s {output.transposons_fasta} ]; then
                echo "[ERROR] Error al descomprimir el archivo de transposones" >> {log}
                exit 1
            fi
            
            echo "Archivo descomprimido correctamente" >> {log}
        else
            echo "El archivo {output.transposons_fasta} ya existe, omitiendo descompresión" >> {log}
        fi
        
        # Detectar el directorio de abricate en el entorno conda activo
        if [ -n "$CONDA_PREFIX" ]; then
            ABRICATE_DB_DIR="$CONDA_PREFIX/db/transposons"
        else
            echo "[ERROR] No se detectó entorno conda activo" >> {log}
            exit 1
        fi
        
        echo "Directorio de abricate detectado: $ABRICATE_DB_DIR" >> {log}
        
        # Crear directorio para la base de datos de transposons
        mkdir -p "$ABRICATE_DB_DIR"
        
        # Copiar archivo de secuencias descomprimido
        cp {output.transposons_fasta} "$ABRICATE_DB_DIR/sequences"
        
        # Configurar abricate para que reconozca la nueva base de datos
        echo "Configurando abricate..." >> {log}
        abricate --setupdb >> {log} 2>&1
        
        # Verificar que la base de datos se configuró correctamente
        echo "Verificando configuración..." >> {log}
        abricate --list >> {log} 2>&1
        
        if ! abricate --list | grep -q "transposons"; then
            echo "[WARNING] La base de datos 'transposons' no aparece en la lista de abricate" >> {log}
            echo "Esto puede ser normal si es la primera vez que se configura" >> {log}
        else
            echo "Base de datos de transposons detectada correctamente" >> {log}
        fi
        
        echo "Configuración de transposons completada" >> {log}
        
        # Crear directorio para el flag si no existe
        mkdir -p $(dirname {output.flag})
        touch {output.flag}
        """

rule setup_platon_database:
    """Descarga y configura la base de datos de Platon si no existe"""
    output:
        db_dir=directory(f"{PLATON_DB_DIR}/db"),
        flag=PLATON_DB_FLAG,
    log:
        PLATON_DB_LOG,
    conda:
        "../envs/epibac_mge.yml"
    resources:
        mem_mb=2000,
        walltime="00:30:00"
    shell:
        """
        # Crear directorio para logs si no existe
        mkdir -p $(dirname {log})
        
        # Verificar si la base de datos ya existe
        if [ -d "{PLATON_DB_DIR}/db" ]; then
            echo "Base de datos de Platon ya existe, creando flag..." &> {log}
            touch {output.flag}
            exit 0
        fi
        
        echo "Descargando base de datos de Platon..." &> {log}
        
        # Crear directorio base
        mkdir -p {PLATON_DB_DIR}
        
        # Descargar y extraer la base de datos
        wget https://zenodo.org/record/4066768/files/db.tar.gz \
          -O {PLATON_DB_DIR}/db.tar.gz >> {log} 2>&1
        
        tar -xzf {PLATON_DB_DIR}/db.tar.gz \
          -C {PLATON_DB_DIR} >> {log} 2>&1
        
        # Limpiar archivo comprimido
        rm -f {PLATON_DB_DIR}/db.tar.gz
        
        # Verificar que la base de datos se extrajo correctamente
        if [ ! -d {output.db_dir} ]; then
            echo "[ERROR] La base de datos de Platon no se pudo extraer correctamente" >> {log}
            exit 1
        fi
        
        echo "Base de datos de Platon configurada correctamente" >> {log}
        touch {output.flag}
        """