# Documentación EPIBAC v.1.2.1

## Índice

1.  [Link documentación y pipeline](#link-documentación-y-pipeline)
2.  [Puesta a punto EPIBAC](#puesta-a-punto-epibac)
3.  [Estructura de Carpetas](#estructura-de-carpetas)
4.  [Nombre carrera / análisis al correr EPIBAC](#nombre-carrera--análisis-al-correr-epibac)
5.  [Fichero metadatos INPUT ejemplo](#fichero-metadatos-input-ejemplo)
6.  [Lanzar EPIBAC](#lanzar-epibac)
7.  [Contacto](#contacto)

---

## 1. Link documentación y pipeline

-   **Documentación completa:** [https://alesanzdro.github.io/epibac-docs](https://alesanzdro.github.io/epibac-docs)
-   **Pipeline (GitHub):** [https://github.com/EpiMol/epibac](https://github.com/EpiMol/epibac)

---

## 2. Puesta a punto EPIBAC

Actualmente EPIBAC v.1.2.1 solo está disponible a través de **Conda**, por lo que no se podría ejecutar con Singularity.

> #### Opciones de instalación
>
> Esta guía contempla tres escenarios. En cualquier caso, desde una instalación desde cero o partiendo de una instalación previa de Conda, **obligatoriamente habrá que instalar el ambiente base de EPIBAC**:
>
> -   **Instalación desde cero:** Si no tiene Conda instalado
> -   **Instalación con Conda existente:** Si ya tiene una instalación de Conda previa
> -   **Instalación directa del ambiente EPIBAC:** Si solo necesita el ambiente de trabajo

### Opción 1: Instalación de Conda desde cero

Si no tiene Conda instalado, puede usar estos comandos para una instalación limpia:

**Descargar e instalar Miniconda**

```bash
# Descargar Miniconda
wget -q [https://repo.anaconda.com/miniconda/Miniconda3-py312_25.1.1-2-Linux-x86_64.sh](https://repo.anaconda.com/miniconda/Miniconda3-py312_25.1.1-2-Linux-x86_64.sh) -O ~/miniconda.sh

# Instalación desatendida en ~/miniconda3
bash ~/miniconda.sh -b -p ~/miniconda3

# Eliminar instalador
rm -f ~/miniconda.sh
```

**Configuración inicial de Conda**

```bash
# Cargar Conda en el shell actual
source ~/miniconda3/etc/profile.d/conda.sh

# Configurar canales en el orden correcto
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Actualizar Conda y paquetes base
conda update -n base conda -y
conda update --all -y

# Instalar Mamba para gestión más rápida de paquetes
conda install -n base -c conda-forge mamba -y
```

> #### Configurar Conda para inicio automático
>
> Para que Conda se inicie automáticamente en nuevas terminales, ejecute:
>
> ```bash
> # Configurar conda para que se inicie automáticamente
> conda init
>
> # IMPORTANTE: Después de ejecutar conda init, debe cerrar la terminal actual 
> # y abrir una nueva para que los cambios surtan efecto
> ```
>
> Este comando modificará su archivo `~/.bashrc`. Verá el prefijo `(base)` delante de su prompt: `(base) usuario@máquina:$`

**Nota:** Si está detrás de un proxy corporativo (como el de la GVA), configure el proxy:

```bash
# Configurar proxy GVA (solo si es necesario)
PROXY_URL="[http://proxy.san.gva.es:8080](http://proxy.san.gva.es:8080)"
conda config --set proxy_servers.http $PROXY_URL
conda config --set proxy_servers.https $PROXY_URL
```

> **Advertencia:** Recuerde instalar el entorno de trabajo base según se explica en **Opción 3**.

### Opción 2: Si ya tiene Conda instalado

Asegúrese de que su instalación de Conda está correctamente configurada:

```bash
# Actualizar y reorganizar canales
conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda config --set channel_priority strict && \
conda update -n base conda -y && \
conda update --all -y && \
conda install -n base -c conda-forge mamba -y
```

> **Advertencia:** Si experimenta conflictos, puede usar la opción `--override-channels`:
>
> ```bash
> conda install --override-channels -n base -c conda-forge mamba
> ```
>
> **Advertencia:** Recuerde instalar el entorno de trabajo base según se explica en **Opción 3**.

### Opción 3: Instalación directa del ambiente EPIBAC

> #### Versiones probadas y compatibles
>
> -   **Conda:** 25.5.1 o superior
> -   **Mamba:** 2.2.0 o superior
>
> Para verificar las versiones, ejecute:
>
> ```bash
> # Verificar versión de Conda
> conda -V
> 
> # Verificar versión de Mamba
> mamba --version
> ```

Cree el ambiente de trabajo para EPIBAC:

```bash
# Crear ambiente 'epibac' con todas las dependencias necesarias
mamba create -n epibac -y -c conda-forge bioconda::snakemake=9.1.1 bioconda::snakemake-wrapper-utils=0.7.2 pandas openpyxl git
```

### Activar el Entorno de Trabajo

Antes de ejecutar cualquier comando de EPIBAC, active el entorno:

```bash
conda activate epibac
```

> #### ¿Por qué usar `conda activate`?
>
> Aunque instalamos **mamba** por su velocidad, la gestión de entornos (activar/desactivar) sigue siendo responsabilidad de `conda`. Snakemake utilizará mamba automáticamente para instalar los paquetes del pipeline.

### Descargar y Acceder al Repositorio

> #### ¿Tiene una instalación previa de EPIBAC?
>
> **Cambios importantes desde la versión 1.0.0:** La estructura del pipeline ha sido modificada. **Se recomienda realizar una instalación limpia.**
>
> **Opción A: Actualizar (Solo para v1.1.x o superior)**
>
> ```bash
> cd epibac
> git pull
> ./epibac.py --version
> ```
>
> **Opción B: Nueva instalación (RECOMENDADO para v1.0.0)**
>
> ```bash
> mv epibac/ epibac_old/
> git clone [https://github.com/EpiMol/epibac.git](https://github.com/EpiMol/epibac.git)
> cd epibac
> ```

Si es una **instalación nueva**, clone el repositorio:

```bash
git clone [https://github.com/EpiMol/epibac.git](https://github.com/EpiMol/epibac.git)
cd epibac
```

> **Importante:** Lance siempre los análisis desde el directorio raíz del pipeline (`~/epibac/`) para evitar la duplicación de bases de datos.

### Instalar las Bases de Datos (solo la primera vez)

Este comando descargará y configurará todas las bases de datos necesarias.

```bash
./epibac.py setup --conda --threads 8
```

### Modificar fichero `config.yaml`

Es crucial modificar la ruta `storage_cabinet` en el fichero `config.yaml` para que EPIBAC pueda subir los resultados a la cabina de la Generalitat Valenciana.

```bash
#...
mode_config:
  gva:
    storage_cabinet: "/mnt/CabinaCS/NLSAR/deposit"
#...
```

---

## 3. Estructura de Carpetas

```bash
epibac/
├── data/
│   └── $run_name/
│       ├── fastq/
│       │   ├── sample1_R1.fastq.gz
│       │   └── sample1_R2.fastq.gz
│       └── samplesinfo_$run_name.csv
├── resources/
│   ├── databases/
│   └── setup_conda_singularity.sh
├── output/
│   └── $run_name/
├── workflow/
├── config.yaml
└── epibac.py
```

-   **data:** Contiene los datos brutos (FASTQ) y el fichero `samplesinfo`.
-   **resources:** Almacena las bases de datos y scripts de configuración.
-   **output:** Carpeta para los resultados.
-   **workflow:** Contiene la lógica del pipeline (reglas de Snakemake).
-   **config.yaml / epibac.py:** Fichero de configuración y script principal de ejecución.

---

## 4. Nombre carrera / análisis al correr EPIBAC

EPIBAC tiene dos modos de funcionamiento:

-   **Modo normal:** Menos restricciones en el formato del nombre de la carrera.
-   **Modo gva (por defecto):** Diseñado para la Vigilancia Genómica de la Comunidad Valenciana.

En modo GVA, el nombre de la carrera (`run_name`) debe seguir el formato estricto: **AAMMDD_HOSPXXX**.

-   **AAMMDD:** Año, mes y día (ej: `250319`).
-   **HOSP:** Código de cuatro letras del hospital (ej: `ALIC`).
-   **XXX:** Contador incremental (ej: `001`).

Un nombre válido sería: `250319_ALIC001`.

---

## 5. Fichero metadatos INPUT ejemplo

El fichero `samplesinfo.csv` debe tener el siguiente formato:

```bash
CODIGO_MUESTRA_ORIGEN;PETICION;FECHA_TOMA_MUESTRA;ESPECIE_SECUENCIA;MOTIVO_WGS;NUM_BROTE;CONFIRMACION;COMENTARIO_WGS;ILLUMINA_R1;ILLUMINA_R2;NANOPORE;MODELO_DORADO
234512;90000001;2025-03-11;Klebsiella pneumoniae;VIGILANCIA;;;;/FULL/PATH/data/250425_GRAL001/fastq/234512_S45_R1_001.fastq.gz;/FULL/PATH/data/250425_GRAL001/fastq/234512_S45_R2_001.fastq.gz;;
234518;90000002;2025-03-11;Escherichia coli;VIGILANCIA;;;;/FULL/PATH/data/250425_GRAL001/fastq/234518_S48_R1_001.fastq.gz;/FULL/PATH/data/250425_GRAL001/fastq/234518_S48_R2_001.fastq.gz;;
234656;90000003;2025-03-11;Pseudomonas putida;VIGILANCIA;;;;/FULL/PATH/data/250425_GRAL001/fastq/234656_S69_R1_001.fastq.gz;/FULL/PATH/data/250425_GRAL001/fastq/234656_S69_R2_001.fastq.gz;;
```

### Descripción de las Columnas

-   **CODIGO_MUESTRA_ORIGEN:** Código único de la muestra.
-   **PETICION:** Código de GestLab.
-   **FECHA_TOMA_MUESTRA:** Fecha de recolección (`AAAA-MM-DD`).
-   **ESPECIE_SECUENCIA:** Nombre científico exacto (ej: `Klebsiella pneumoniae`).
-   **MOTIVO_WGS:** `VIGILANCIA` o `BROTE`.
-   **NUM_BROTE:** Número de brote si aplica.
-   **CONFIRMACION:** Información a confirmar por FISABIO.
-   **COMENTARIO_WGS:** Comentarios del laboratorio.
-   **ILLUMINA_R1 / ILLUMINA_R2:** Ruta completa a los ficheros FASTQ.
-   **NANOPORE:** Ruta al fichero Nanopore (si existe).
-   **MODELO_DORADO:** Modelo de basecalling de Dorado.

---

## 6. Lanzar EPIBAC

### Definir el nombre de la carrera

```bash
RUN="250319_ALIC001"
```

### Creación de plantilla (Opcional)

Crea una plantilla `samplesinfo.csv` a partir de un directorio de ficheros FASTQ.

```bash
./epibac.py samplesinfo --run_name <span class="math-inline">\{RUN\} \-\-platform illumina \-\-fastq data/</span>{RUN}/fastq
```

### Validar la plantilla de la carrera

Una vez rellenada, valide la plantilla. Este paso es **obligatorio**.

```bash
./epibac.py validate --samples data/<span class="math-inline">\{RUN\}/samplesinfo\_</span>{RUN}.csv --outdir output/${RUN}
```

### Lanzar análisis

Se recomienda lanzar EPIBAC desde el directorio raíz del proyecto.

```bash
./epibac.py run --conda --threads 24 --samples data/<span class="math-inline">\{RUN\}/samplesinfo\_</span>{RUN}.csv --outdir output/${RUN} --run_name ${RUN} --resume
```

-   `--threads`: Número de hilos a usar.
-   `--samples`: Ruta al fichero de metadatos validado.
-   `--outdir`: Directorio de salida.
-   `--run_name`: Nombre de la carrera.
-   `--resume`: Reanuda ejecuciones interrumpidas.

### Ejemplo práctico: Análisis con 3 muestras de prueba

Este ejemplo descarga y analiza 3 muestras desde Zenodo.

```bash
# Definimos variables
RUN="250425_GRAL001"
EPIBAC_PATH=~/epibac
DATA_DIR="${EPIBAC_PATH}/data"
ZIP_URL="https://zenodo.org/records/15633357/files/250425_GRAL001.zip"
ZIP_FILE="${DATA_DIR}/${RUN}.zip"
SAMPLESINFO="${DATA_DIR}/${RUN}/samplesinfo_${RUN}.csv"

# Creamos el directorio si no existe
mkdir -p ${DATA_DIR}

# Descargamos el set de ejemplo
wget -O ${ZIP_FILE} ${ZIP_URL}

# Descomprimimos dentro de data/
unzip -o ${ZIP_FILE} -d ${DATA_DIR}

# Eliminamos fichero descargado
rm ${ZIP_FILE}

# Entramos al directorio del pipeline
cd ${EPIBAC_PATH}

# Activamos el entorno conda
conda activate epibac

# Corregimos el path en el archivo samples_info
sed -i "s|/FULL/PATH|$(pwd)|g" "${SAMPLESINFO}"

# Mostramos el contenido del archivo de metadatos
cat ${SAMPLESINFO}

# Validamos el archivo samplesinfo
./epibac.py validate --samples "${SAMPLESINFO}" --outdir "output/${RUN}"

# Si todo está bien, deberías ver:
# "The sample file has been successfully validated!"

# Ejecutamos el análisis de prueba
./epibac.py run --conda --threads 24 --samples "${SAMPLESINFO}" --outdir "output/${RUN}" --run_name "${RUN}" --resume
```

### Borrar archivos temporales (Opcional)

```bash
./epibac.py clean
```

---

## 7. Contacto

Para cualquier duda contactar con:

-   [vigilancia.genomica@fisabio.es](mailto:vigilancia.genomica@fisabio.es)
-   [fernando.gonzalez@fisabio.es](mailto:fernando.gonzalez@fisabio.es)
-   [alejandro.sanz@fisabio.es](mailto:alejandro.sanz@fisabio.es)
