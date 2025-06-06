import glob
import re
import pandas as pd
import os
from datetime import datetime

# =================== CONSTANTES Y VARIABLES GLOBALES =================== #

# Fecha actual para reportes
DATE = datetime.now().strftime("%y%m%d")

# =================== OPCIONES DE CONFIGURACIÓN =================== #


# Procesar opciones de línea de comandos
def process_cli_args():
    """Procesa los argumentos de línea de comandos y actualiza la configuración."""
    import sys

    # Detectar skip_prokka
    if "--skip_prokka" in sys.argv:
        if "skip" not in config:
            config["skip"] = {}
        config["skip"]["prokka"] = True
        sys.argv.remove("--skip_prokka")


# Ejecutar procesamiento de argumentos
process_cli_args()


def should_skip(component):
    """
    Verifica si un componente debe ser omitido según la configuración.

    Args:
        component: Nombre del componente a verificar (ej. "prokka")

    Returns:
        bool: True si se debe omitir, False en caso contrario
    """
    # Primero verifica el formato antiguo para compatibilidad
    if f"skip_{component}" in config:
        return config[f"skip_{component}"]
    # Luego verifica el nuevo formato
    return config.get("skip", {}).get(component, False)


def check_dependencies(self):
    """Verifica que todas las dependencias necesarias estén instaladas."""
    self.logger.info("Comprobando dependencias...")

    deps = {
        "snakemake": "snakemake --version",
        "conda": "conda --version" if self.args.conda else None,
        "singularity": "singularity --version" if self.args.singularity else None,
    }

    missing = []
    for name, cmd in deps.items():
        if cmd is None:
            continue
        try:
            subprocess.run(
                cmd.split(), check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            self.logger.debug(f"✓ {name} instalado correctamente")
        except (subprocess.SubprocessError, FileNotFoundError):
            missing.append(name)
            self.logger.error(f"✗ {name} no encontrado o no funciona correctamente")

    if missing:
        self.logger.error(f"Dependencias faltantes: {', '.join(missing)}")
        self.logger.error(
            "Por favor, instala las dependencias faltantes antes de continuar"
        )
        return False

    self.logger.info("Todas las dependencias están instaladas correctamente")
    return True


# =================== CONFIGURACIÓN DE BASES DE DATOS =================== #

# ----- KRAKEN2 -----
KRAKEN_DB_URL = config["params"]["kraken2"]["db_url"]
KRAKEN_DB_NAME = os.path.basename(KRAKEN_DB_URL).replace(".tar.gz", "")
KRAKEN_DB_DIR = f"{REFDIR}/databases/kraken2/{KRAKEN_DB_NAME}"
KRAKEN_DB_FLAG = f"{KRAKEN_DB_DIR}/.installed.flag"
KRAKEN_DB_LOG = f"{REFDIR}/databases/log/{KRAKEN_DB_NAME}.log"

# ----- AMRFINDER -----
AMRFINDER_DB_DIR = f"{REFDIR}/databases/amrfinder"
AMRFINDER_DB_FLAG = f"{AMRFINDER_DB_DIR}/.installed.flag"
AMRFINDER_DB_LOG = f"{REFDIR}/databases/log/amrfinder.log"

# ----- RESFINDER -----
RESFINDER_DB_DIR = f"{REFDIR}/databases/resfinder"
RESFINDER_DB_FLAG = f"{RESFINDER_DB_DIR}/.installed.flag"
RESFINDER_DB_LOG = f"{REFDIR}/databases/log/resfinder.log"

# ----- PROKKA -----
PROKKA_DB_DIR = f"{REFDIR}/databases/prokka"
PROKKA_DB_FLAG = f"{PROKKA_DB_DIR}/.installed.flag"
PROKKA_DB_LOG = f"{REFDIR}/databases/log/prokka.log"

# ----- ABRICATE -----
TRANSPOSONS_DB_DIR = f"{REFDIR}/databases/transposons"
TRANSPOSONS_DB_FLAG = f"{TRANSPOSONS_DB_DIR}/.installed.flag"
TRANSPOSONS_DB_LOG = f"{REFDIR}/databases/log/transposons.log"

# ----- PLATON -----
PLATON_DB_DIR = f"{REFDIR}/databases/platon"
PLATON_DB_FLAG = f"{PLATON_DB_DIR}/.installed.flag"
PLATON_DB_LOG = f"{REFDIR}/databases/log/platon.log"

# ----- MLST -----
MLST_DB_DIR = f"{REFDIR}/databases/mlst"


# =================== FUNCIONES DE ACCESO A DATOS DE MUESTRAS =================== #


def get_samples_safe():
    """
    Obtiene el DataFrame de muestras de forma segura, sin lanzar errores si el archivo no existe.

    Returns:
        DataFrame: DataFrame con información de las muestras o vacío si no existe.
    """
    samples_csv = f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv"
    if os.path.exists(samples_csv):
        use_column = config.get("primary_id_column", "id")
        return pd.read_csv(samples_csv, sep=";", dtype=str).set_index(
            use_column, drop=False
        )
    else:
        return pd.DataFrame()  # Devuelve vacío, evita romper el flujo


def get_samples():
    """
    Obtiene el DataFrame de muestras. Lanza error si el archivo no existe.

    Returns:
        DataFrame: DataFrame con información de las muestras.

    Raises:
        FileNotFoundError: Si el archivo de muestras validadas no existe.
    """
    samples_csv = f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv"
    if os.path.exists(samples_csv):
        use_column = config.get("primary_id_column", "id")
        return pd.read_csv(samples_csv, sep=";", dtype=str).set_index(
            use_column, drop=False
        )
    else:
        raise FileNotFoundError(
            f"El fichero validado '{samples_csv}' aún no existe. Ejecuta primero la regla 'validate_samples'."
        )


def get_sample_index_if_exists():
    """
    Obtiene la lista de IDs de muestras si el archivo existe.

    Returns:
        list: Lista de IDs o lista vacía si el archivo no existe.
    """
    samples_csv = f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv"
    if os.path.exists(samples_csv):
        use_column = config.get("primary_id_column", "id")
        df = pd.read_csv(samples_csv, sep=";", dtype=str)
        return df[use_column].tolist()
    else:
        return []


def get_fastq(wildcards):
    """
    Obtiene las rutas de los archivos FASTQ para una muestra.

    Args:
        wildcards: Wildcards de Snakemake con el nombre de la muestra.

    Returns:
        dict: Diccionario con rutas de los archivos R1 y R2.

    Raises:
        ValueError: Si no se encuentran archivos FASTQ válidos.
    """
    samples = get_samples()
    fastqs = samples.loc[wildcards.sample, ["illumina_r1", "illumina_r2"]]
    if pd.notnull(fastqs["illumina_r1"]) and pd.notnull(fastqs["illumina_r2"]):
        return {"r1": fastqs["illumina_r1"], "r2": fastqs["illumina_r2"]}
    elif pd.notnull(fastqs["illumina_r1"]):
        return {"r1": fastqs["illumina_r1"]}
    else:
        raise ValueError(
            f"No se encontraron archivos FASTQ válidos para {wildcards.sample}"
        )


def get_filtered_samples():
    """
    Obtiene muestras que tienen archivos de validación.

    Returns:
        list: Lista de nombres de muestras validadas.
    """
    validated_samples = [
        f.split("/")[-1].split(".")[0]
        for f in glob.glob(f"{OUTDIR}/validated/*.validated")
    ]
    return validated_samples


def get_dorado_model(wildcards):
    """
    Obtiene el modelo Dorado para una muestra específica.
    
    Args:
        wildcards: Wildcards de Snakemake con el nombre de la muestra.
    
    Returns:
        str: Modelo Dorado a usar para esta muestra.
    
    Raises:
        ValueError: Si no se encuentra un modelo válido.
    """
    samples = get_samples()
    sample_row = samples.loc[wildcards.sample]
    
    # Verificar si la muestra tiene columna dorado_model y un valor específico
    if "dorado_model" in samples.columns and pd.notna(sample_row.get("dorado_model")) and str(sample_row.get("dorado_model")).strip():
        return str(sample_row.get("dorado_model")).strip()
    
    # Usar modelo global del config.yaml
    global_model = config.get("params", {}).get("nanopore", {}).get("dorado_model", None)
    if global_model:
        return global_model
    
    # Si no hay modelo disponible, lanzar error
    raise ValueError(f"No se encontró modelo Dorado para la muestra {wildcards.sample}")


def has_nanopore_data(wildcards):
    """
    Verifica si una muestra tiene datos de nanopore.
    
    Args:
        wildcards: Wildcards de Snakemake con el nombre de la muestra.
    
    Returns:
        bool: True si la muestra tiene datos de nanopore.
    """
    samples = get_samples()
    sample_row = samples.loc[wildcards.sample]
    
    nanopore_path = sample_row.get("nanopore")
    return pd.notna(nanopore_path) and str(nanopore_path).strip() != ""


def get_nanopore_fastq(wildcards):
    """
    Obtiene la ruta del archivo FASTQ de nanopore para una muestra.
    
    Args:
        wildcards: Wildcards de Snakemake con el nombre de la muestra.
    
    Returns:
        str: Ruta al archivo FASTQ de nanopore.
    
    Raises:
        ValueError: Si no se encuentran archivos nanopore válidos.
    """
    samples = get_samples()
    sample_row = samples.loc[wildcards.sample]
    
    nanopore_path = sample_row.get("nanopore")
    if pd.notna(nanopore_path) and str(nanopore_path).strip():
        return str(nanopore_path).strip()
    else:
        raise ValueError(f"No se encontró archivo nanopore para la muestra {wildcards.sample}")


# =================== UTILIDADES GENERALES =================== #


def get_resource(rule, resource, default=None):
    """
    Obtiene un recurso para una regla específica, con fallback a los valores por defecto.

    Args:
        rule: Nombre de la regla.
        resource: Tipo de recurso (threads, mem, walltime, gpu).
        default: Valor por defecto a usar si no se encuentra en la configuración.

    Returns:
        Valor del recurso.
    """
    try:
        return config["resources"][rule][resource]
    except KeyError:
        try:
            return config["resources"]["default"][resource]
        except KeyError:
            if default is not None:
                return default
            raise ValueError(f"No se encontró el recurso '{resource}' para la regla '{rule}' y no se proporcionó valor por defecto")


def sanitize_id(x):
    """
    Limpia un ID para ser utilizado como wildcard.

    Args:
        x: ID a limpiar.

    Returns:
        str: ID limpio.
    """
    return re.sub(r"[ .,-]", "_", str(x))


def get_all_inputs():
    """
    Genera la lista de archivos de salida para la regla 'all'.

    Returns:
        list: Lista de archivos esperados por la regla all.
    """
    inputs = [
        f"{LOGDIR}/samplesinfo/samplesinfo_validated.csv",
        f"{OUTDIR}/qc/multiqc.html",
        f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.tsv",
        f"{OUTDIR}/report/{TAG_RUN}_EPIBAC.xlsx",
    ]

    # Expandir archivos de platon para todas las muestras
    sample_ids = get_sample_index_if_exists()
    if sample_ids:
        inputs.extend([
            f"{OUTDIR}/mge_analysis/plasmids/platon/{sample}/{sample}.tsv"
            for sample in sample_ids
        ])
        inputs.extend([
            f"{OUTDIR}/mge_analysis/plasmids/mob_suite/{sample}/mobtyper_results.txt"
            for sample in sample_ids
        ])
        inputs.extend([
             f"{OUTDIR}/mge_analysis/transposons/{sample}_abricate.tsv"
            for sample in sample_ids
        ])
        
        # Agregar análisis unificado de plásmidos
        inputs.extend([
            f"{OUTDIR}/mge_analysis/unified_plasmid_analysis/{TAG_RUN}_plasmid_analysis_summary.tsv",
            f"{OUTDIR}/mge_analysis/unified_plasmid_analysis/{TAG_RUN}_Reporte_Plasmidos_Simplificado.xlsx"
        ])      
       

    # Si no estamos omitiendo Prokka, añadimos su metadata
    if not should_skip("prokka"):
        inputs.append(f"{REFDIR}/databases/prokka/VERSION.txt")

    # Si estamos en modo GVA, añadimos también el reporte para GESTLAB
    if config.get("mode") == "gva":
        inputs.append(f"{OUTDIR}/report/{TAG_RUN}_EPIBAC_GESTLAB.csv")
        inputs.append(f"{OUTDIR}/report/{TAG_RUN}_EPIBAC_GESTLAB.xlsx")
        inputs.append(f"{OUTDIR}/report/{TAG_RUN}_file_copy_log.txt")

    return inputs


# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",
