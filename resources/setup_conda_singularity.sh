#!/usr/bin/env bash
#
# setup_env.sh
#
# Multifunction script to:
#   - Uninstall and install Singularity v3.8.7 (option 'singularity').
#   - COMPLETELY Uninstall any conda/Miniconda/Anaconda installation
#     and reinstall the latest version of Miniconda + basic configuration (option 'conda').
#
# Usage:
#   1) sudo ./setup_env.sh singularity
#   2) sudo ./setup_env.sh conda
#   3) sudo ./setup_env.sh conda http://proxy.ejemplo.com:8080
#
# WARNING!: In "conda" mode, all environments and the conda installation are deleted.
#             Make a backup if necessary.
#
# Tested on Debian/Ubuntu-like systems with 'apt'.
# -----------------------------------------------------------------------------
# Reformatted the script using: shfmt -i 4 -ci -w setup_conda_singularity.sh
# -----------------------------------------------------------------------------

# Color and format definitions
COLOR_RESET="\033[0m"
BOLD="\033[1m"
QUESTION="${BOLD}\033[33m"      # Bold golden/yellow
INFO="${BOLD}\033[90m"          # Bold gray
ERROR="${BOLD}\033[31m"         # Bold red
DANGER="${BOLD}\033[31m"        # Bold red
WARNING="${BOLD}\033[38;5;208m" # Bold orange (requires 256-color support)
TIP="${BOLD}\033[38;5;208m"     # Bold orange (requires 256-color support)
OK="${BOLD}\033[32m"            # Bold green
END="${BOLD}\033[32m"           # Bold green

# Constantes para las etiquetas de configuración - MOVER AL PRINCIPIO DEL SCRIPT
EPIBAC_TAG_START="# >>> EPIBAC_CONFIG >>>"
EPIBAC_TAG_END="# <<< EPIBAC_CONFIG <<<"
CONDA_TAG_START="# >>> conda initialize >>>"
CONDA_TAG_END="# <<< conda initialize <<<"

# Function to use in echo with colors
format_message() {
    local prefix="$1"
    local message="$2"
    local color="$3"

    echo -e "${color}[${prefix}]${COLOR_RESET} ${message}"
}

# Function to clean duplicates in PATH
clean_path() {
    local old_path="$PATH"
    PATH=$(echo "$old_path" | awk -v RS=':' '!a[$1]++ {if (NR > 1) printf ":"; printf $1}')
    export PATH
    format_message "INFO" "Cleaned duplicate entries in PATH." "$INFO"
}

# Call clean_path at the beginning of the script to ensure a clean PATH
clean_path

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 {singularity|conda|uninstall} [proxy_url]"
    echo "Examples:"
    echo "  $0 singularity"
    echo "  $0 conda"
    echo "  $0 conda http://proxy.san.gva.es:8080"
    echo "  $0 uninstall       # Remove previous installations and configurations"
    exit 1
fi

MODE="$1"
PROXY_URL=""

# If there is a second parameter, we assume it is the proxy
if [[ $# -eq 2 ]]; then
    PROXY_URL="$2"
    format_message "INFO" "Proxy will be configured: $PROXY_URL" "$INFO"
fi

# Detect the user calling sudo to edit their .bashrc, .zshrc files
REAL_USER="${SUDO_USER:-root}"
REAL_USER_HOME="$(eval echo ~"$REAL_USER")"

# Inicializar variables de instalación al principio del script
CONDA_INSTALLED=false
SINGULARITY_INSTALLED=false
GO_INSTALLED=false

# Inicializar las variables de rastreo al principio del script
GO_INSTALLED_BY_SCRIPT=false
SINGULARITY_INSTALLED_BY_SCRIPT=false
CONDA_INSTALLED_BY_SCRIPT=false

# Inicializar variables para el entorno Python
SNAKE_ENV_ALIAS_REQUESTED=false
SNAKE_ENV_PATH=""

echo -e "${BOLD}======================================${COLOR_RESET}"
echo -e "${BOLD} Script setup_env.sh - Option: $MODE ${COLOR_RESET}"
echo -e "${BOLD} Real User:       $REAL_USER${COLOR_RESET}"
echo -e "${BOLD} HOME Directory:    $REAL_USER_HOME${COLOR_RESET}"
echo -e "${BOLD}======================================${COLOR_RESET}"

# Function to ask questions with answer validation
ask_question() {
    local prompt="$1"  # Question message
    local default="$2" # Default value (Y or N)

    while true; do
        # Color the question
        read -rp "$(echo -e "${QUESTION}[QUESTION]${COLOR_RESET} $prompt")" answer

        # If it's empty, use the default value
        if [[ -z "$answer" ]]; then
            answer="$default"
        fi

        # Normalize answer to lowercase
        answer=$(echo "$answer" | tr '[:upper:]' '[:lower:]')

        case "$answer" in
            y | s | si | yes) return 0 ;; # Success = Yes
            n | no) return 1 ;;           # Failure = No
            q | quit | exit)
                format_message "INFO" "Exiting script by user request." "$INFO"
                exit 0
                ;;
            *) format_message "INFO" "Please answer (y)es, (n)o or (q) to quit." "$INFO" ;;
        esac
    done
}

# Mejorar la función para ejecutar comandos como el usuario real sin pedir contraseña
run_as_real_user() {
    local cmd="$1"
    # Si estamos ejecutando como root, usamos su para ejecutar como el usuario real
    if [[ $EUID -eq 0 ]]; then
        su -l "$REAL_USER" -c "$cmd"
    else
        # Si no somos root, simplemente ejecutamos el comando
        bash -c "$cmd"
    fi
}

###############################
# FUNCTIONS: UNINSTALLATION   #
###############################

# ----- Check disk space -----
check_disk_space() {
    local required_space_mb=30000 # ~5GB for the complete installation
    local available_space_mb=$(df -m /usr/local | awk 'NR==2 {print $4}')

    format_message "INFO" "Available space: ${available_space_mb}MB, required for complete BASE installation: ${required_space_mb}MB" "$INFO"

    if [[ $available_space_mb -lt $required_space_mb ]]; then
        format_message "ERROR" "Insufficient disk space. At least ${required_space_mb}MB are needed" "$ERROR"
        return 1
    fi

    return 0
}

# Use the function before installation
echo
if ! check_disk_space; then
    if ask_question "Continue despite the space warning? (y/N): " "n"; then
        format_message "INFO" "Continuing with limited space..." "$INFO"
    else
        format_message "INFO" "Installation cancelled." "$INFO"
        exit 1
    fi
fi

# ----- Uninstall Go -----
remove_go() {
    echo
    format_message "INFO" "Removing Go..." "$INFO"

    # 1. Remove apt packages
    format_message "INFO" "Searching for packages 'golang-go', 'golang-*' in apt..." "$INFO"
    GO_PKGS="$(dpkg -l | grep -E 'golang-go|golang-[0-9\.]+|golang-doc' || true)"
    if [[ -n "$GO_PKGS" ]]; then
        format_message "INFO" "Go packages found installed:" "$INFO"
        echo "$GO_PKGS"
        format_message "INFO" "Proceeding to uninstall them..." "$INFO"
        apt remove --purge -y golang-go golang-doc || true
        apt autoremove --purge -y || true
    else
        format_message "INFO" "No Go packages installed with apt (or not found)." "$INFO"
    fi

    # 2. Delete /usr/local/go if it exists (manual installations)
    if [[ -d /usr/local/go ]]; then
        format_message "INFO" "Deleting /usr/local/go..." "$INFO"
        rm -rf /usr/local/go
    fi

    # 3. Delete ~/go (GOPATH) directory if the user wishes
    if [[ -d "$REAL_USER_HOME/go" ]]; then
        echo
        if ask_question "Also delete the folder '$REAL_USER_HOME/go' (GOPATH)? (y/n): " "n"; then
            rm -rf "$REAL_USER_HOME/go"
            format_message "OK" "$REAL_USER_HOME/go was deleted." "$OK"
        else
            format_message "INFO" "$REAL_USER_HOME/go is kept." "$INFO"
        fi
    fi

    # 4. Remove exports in .bashrc and .zshrc
    echo
    format_message "INFO" "Cleaning Go references in .bashrc / .zshrc of user $REAL_USER..." "$INFO"

    for SHELL_RC in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
        [[ -f "$SHELL_RC" ]] || continue
        # Backup
        cp "$SHELL_RC" "$SHELL_RC.bak_$(date +%Y%m%d%H%M%S)"

        # Use multiple simple patterns instead of a complex one
        sed -i '/usr\/local\/go\/bin/d' "$SHELL_RC"
        sed -i '/GOROOT/d' "$SHELL_RC"
        sed -i '/GOPATH/d' "$SHELL_RC"
    done

    format_message "INFO" "Go references removed from shell configuration files." "$INFO"
    format_message "INFO" "(Restart or open a new shell for it to take effect)." "$INFO"
}

# ----- Install Go (v1.24.1) -----
install_go() {
    echo
    format_message "INFO" "Installing Go 1.24.1..." "$INFO"

    # Make sure we have wget
    apt-get update -y
    apt-get install -y wget

    # We create a temporary working directory
    WORK_DIR="/tmp/go-install.$$"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # Download Go with integrity verification
    format_message "INFO" "Downloading Go 1.24.1..." "$INFO"
    wget -q https://go.dev/dl/go1.24.1.linux-amd64.tar.gz

    # SHA256 verification (optional but recommended)
    # You can get the current value from https://go.dev/dl/
    EXPECTED_SHA256="cb2396bae64183cdccf81a9a6df0aea3bce9511fc21469fb89a0c00470088073"
    ACTUAL_SHA256=$(sha256sum go1.24.1.linux-amd64.tar.gz | cut -d' ' -f1)
    if [[ "$EXPECTED_SHA256" != "$ACTUAL_SHA256" ]]; then
        format_message "ERROR" "Integrity verification failed for go1.24.1.linux-amd64.tar.gz" "$ERROR"
        exit 1
    fi

    # Extract in /usr/local
    format_message "INFO" "Extracting Go in /usr/local..." "$INFO"
    rm -rf /usr/local/go
    tar -C /usr/local -xzf go1.24.1.linux-amd64.tar.gz

    # No añadimos configuración directamente a .bashrc/.zshrc
    # Solo verificamos la instalación
    export PATH=/usr/local/go/bin:$PATH
    if command -v go &>/dev/null; then
        format_message "OK" "Go 1.24.1 installed correctly:" "$OK"
        go version
    else
        format_message "ERROR" "Could not verify Go installation" "$ERROR"
    fi

    # Cleaning
    cd /
    rm -rf "$WORK_DIR"

    # Marcar explícitamente que nosotros lo instalamos
    GO_INSTALLED_BY_SCRIPT=true
}

# ----- Uninstall Singularity -----
remove_singularity() {
    echo
    format_message "INFO" "Proceeding to remove Singularity..." "$INFO"

    # Show the current version if available
    if command -v singularity &>/dev/null; then
        CURRENT_VERSION=$(singularity --version 2>/dev/null || echo "Unknown version")
        format_message "INFO" "Detected version: $CURRENT_VERSION" "$INFO"
    fi

    # 1. Remove singular packages from apt (more specific search)
    format_message "INFO" "Searching for 'singularity' packages in apt..." "$INFO"
    SING_PKGS="$(dpkg -l | grep -E 'singularity|apptainer' || true)"
    if [[ -n "$SING_PKGS" ]]; then
        format_message "INFO" "Installed packages found:" "$INFO"
        echo "$SING_PKGS"
        format_message "INFO" "Uninstalling with apt..." "$INFO"
        apt remove --purge -y singularity-container apptainer || true
        apt autoremove --purge -y || true
    else
        format_message "INFO" "No 'singularity' packages installed with apt." "$INFO"
    fi

    # 2. Remove local binaries
    for BIN in singularity apptainer; do
        if command -v $BIN &>/dev/null; then
            BIN_PATH="$(which $BIN)"
            format_message "INFO" "Deleting binary $BIN_PATH..." "$INFO"
            rm -f "$BIN_PATH"
        fi
    done

    # 3. Remove configuration and library directories
    for DIR in /usr/local/libexec/singularity /usr/local/etc/singularity \
        /usr/local/libexec/apptainer /usr/local/etc/apptainer; do
        if [[ -d "$DIR" ]]; then
            format_message "INFO" "Deleting $DIR..." "$INFO"
            rm -rf "$DIR"
        fi
    done

    format_message "INFO" "Singularity has been removed (or was not found)." "$INFO"
}
# ----- COMPLETELY Uninstall conda/Miniconda/Anaconda -----
remove_conda() {
    echo
    echo -e "${BOLD}\033[38;5;208m!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo -e "${BOLD} WARNING: ALL conda environments will be deleted,       \033[38;5;208m"
    echo -e "${BOLD} the Miniconda/Anaconda installation and the references in  \033[38;5;208m"
    echo -e "${BOLD} .bashrc / .zshrc.                                          \033[38;5;208m"
    echo -e "${BOLD}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!${COLOR_RESET}"
    echo
    if ask_question "Do you really want to delete ALL conda installations? (y/N): " "n"; then
        format_message "INFO" "Proceeding to total conda removal..." "$INFO"

        # 1. Try to deactivate conda (in case it is active)
        #    Avoids some locks, but may not be essential.
        if command -v conda &>/dev/null; then
            conda deactivate || true
        fi

        # 2. Delete typical conda/Anaconda/Miniconda directories
        #    (Some people install it in ~/anaconda3 or ~/miniconda3)
        # Before deleting important directories
        for cdir in "$REAL_USER_HOME/miniconda3" "$REAL_USER_HOME/anaconda3" "$REAL_USER_HOME/conda" "$REAL_USER_HOME/.conda"; do
            if [[ -d "$cdir" && "$cdir" != "/" && "$cdir" != "$REAL_USER_HOME" ]]; then
                format_message "INFO" "Deleting directory $cdir..." "$INFO"
                rm -rf "$cdir"
            elif [[ -d "$cdir" ]]; then
                format_message "DANGER" "$cdir will not be deleted for security reasons" "$DANGER"
            fi
        done

        # Also delete the .condarc file if it exists
        if [[ -f "$REAL_USER_HOME/.condarc" ]]; then
            format_message "INFO" "Deleting configuration file $REAL_USER_HOME/.condarc..." "$INFO"
            # Create backup first
            cp "$REAL_USER_HOME/.condarc" "$REAL_USER_HOME/.condarc.bak_$(date +%Y%m%d%H%M%S)"
            rm -f "$REAL_USER_HOME/.condarc"
        fi

        # 3. Delete lines in .bashrc / .zshrc that do 'conda init' or similar
        for SHELL_RC in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
            [[ -f "$SHELL_RC" ]] || continue
            cp "$SHELL_RC" "$SHELL_RC.bak_$(date +%Y%m%d%H%M%S)"

            # Delete the complete conda initialization block
            sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$SHELL_RC"

            # Also delete loose lines that contain references to conda for security
            sed -i '/conda\.sh\|conda activate\|conda init\|\.conda\|anaconda3\|miniconda3/d' "$SHELL_RC"

            format_message "INFO" "The file $SHELL_RC has been cleaned of references to conda." "$INFO"
        done

        format_message "INFO" "References to conda have been removed from the configuration files." "$INFO"
        format_message "INFO" "For the changes to take effect, close and reopen the terminal (or do a manual 'source')." "$INFO"
    else
        format_message "INFO" "Conda removal aborted." "$INFO"
        return 1 # Importante: devolver un código para que la función que llama sepa que se canceló
    fi
}

###############################
# FUNCTIONS: INSTALLATIONS    #
###############################

# ----- Install Singularity (v3.8.7) -----
install_singularity() {
    echo
    format_message "INFO" "Installing Singularity v3.8.7..." "$INFO"

    # Dependencies
    apt-get update -y
    apt-get install -y build-essential libseccomp-dev pkg-config squashfs-tools cryptsetup runc git

    # Use PID to create a unique temporary directory name
    WORK_DIR="/tmp/singularity-install.$$"
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    # Clone and compile
    git clone https://github.com/hpcng/singularity.git
    cd singularity

    # Verify that the version exists before checking out
    SINGULARITY_VERSION="3.8.7"
    format_message "INFO" "Verifying that v${SINGULARITY_VERSION} exists in the repository..." "$INFO"
    if ! git tag | grep -q "v${SINGULARITY_VERSION}"; then
        format_message "ERROR" "Version ${SINGULARITY_VERSION} not found" "$ERROR"
        exit 1
    fi

    git checkout "v${SINGULARITY_VERSION}"

    ./mconfig
    make -C ./builddir
    make -C ./builddir install

    # Verify the installation
    if command -v singularity &>/dev/null; then
        format_message "OK" "Singularity v${SINGULARITY_VERSION} installed correctly:" "$OK"
        singularity --version
    else
        format_message "ERROR" "Could not verify Singularity installation" "$ERROR"
    fi

    # Cleaning
    cd /
    rm -rf "$WORK_DIR"

    # Marcar explícitamente que nosotros lo instalamos
    SINGULARITY_INSTALLED_BY_SCRIPT=true
}

# ----- Install virtual Python environment with Snakemake -----
install_python_venv() {
    echo
    format_message "INFO" "Installing Python 3.12 and creating virtual environment..." "$INFO"

    # Install Python 3.12 and necessary tools
    apt-get update -y
    apt-get install -y python3.12 python3.12-venv python3.12-dev

    # Create the directory for the virtual environment if it does not exist
    ENV_DIR="$REAL_USER_HOME/snake_env"

    # Create environment as real user (not as root)
    sudo -u "$REAL_USER" python3.12 -m venv "$ENV_DIR"

    # Install packages inside the virtual environment
    sudo -u "$REAL_USER" bash -c "source $ENV_DIR/bin/activate && \
    pip install --upgrade pip && \
    pip install snakemake==9.1.1 snakemake-wrapper-utils==0.7.2 pandas openpyxl gitpython"

    # Verify installation
    if sudo -u "$REAL_USER" bash -c "source $ENV_DIR/bin/activate && snakemake --version"; then
        format_message "OK" "Virtual environment created successfully in $ENV_DIR" "$OK"
    else
        format_message "ERROR" "There was a problem creating the virtual environment" "$ERROR"
        return 1
    fi
    echo
    # No añadimos alias directamente a .bashrc/.zshrc
    if ask_question "Add alias 'snake_env' to activate the environment? (y/N): " "n"; then
        echo
        # Solo marcamos que queremos añadir el alias, sin modificar ficheros
        SNAKE_ENV_ALIAS_REQUESTED=true
        SNAKE_ENV_PATH="$ENV_DIR"
        format_message "INFO" "Alias will be added to EPIBAC configuration block" "$INFO"
    else
        format_message "INFO" "Alias will not be added. To activate the environment use:" "$INFO"
        echo "       source $ENV_DIR/bin/activate"
    fi

    echo
    format_message "INFO" "You can also add this environment to VS Code:" "$INFO"
    echo "  1. Open VS Code"
    echo "  2. Press Ctrl+Shift+P"
    echo "  3. Type 'Python: Select Interpreter'"
    echo "  4. Select the environment in $ENV_DIR/bin/python3.12"

    return 0
}

# ----- Install Miniconda (py312_25.1.1-2) -----
install_conda() {
    echo
    format_message "INFO" "Installing version py312_25.1.1-2 of Miniconda3 for Linux x86_64..." "$INFO"

    # Make sure we have wget
    apt-get update -y
    apt-get install -y wget

    # We download to the real user's home (not /root)
    cd "$REAL_USER_HOME"
    # Usar run_as_real_user en lugar de sudo -u
    run_as_real_user "wget -q https://repo.anaconda.com/miniconda/Miniconda3-py312_25.1.1-2-Linux-x86_64.sh -O $REAL_USER_HOME/miniconda.sh"

    # Verify SHA256 checksum
    EXPECTED_SHA256="4766d85b5f7d235ce250e998ebb5a8a8210cbd4f2b0fea4d2177b3ed9ea87884"
    ACTUAL_SHA256=$(sha256sum "$REAL_USER_HOME/miniconda.sh" | cut -d' ' -f1)
    if [[ "$EXPECTED_SHA256" != "$ACTUAL_SHA256" ]]; then
        format_message "ERROR" "Integrity verification failed for miniconda.sh" "$ERROR"
        exit 1
    fi

    # We give execution permissions
    chmod +x "$REAL_USER_HOME/miniconda.sh"
    # Asegurarse de que el propietario sea el usuario real
    chown "$REAL_USER:$(id -gn $REAL_USER)" "$REAL_USER_HOME/miniconda.sh"

    # Check if directory already exists and remove it if necessary
    if [[ -d "$REAL_USER_HOME/miniconda3" ]]; then
        format_message "WARNING" "Directory $REAL_USER_HOME/miniconda3 already exists" "$WARNING"
        format_message "INFO" "Removing existing directory before installation..." "$INFO"
        rm -rf "$REAL_USER_HOME/miniconda3"
    fi

    # Unattended installation in ~/miniconda3
    run_as_real_user "bash $REAL_USER_HOME/miniconda.sh -b -p $REAL_USER_HOME/miniconda3"

    # Delete installer
    run_as_real_user "rm -f $REAL_USER_HOME/miniconda.sh"

    format_message "INFO" "Miniconda installed in $REAL_USER_HOME/miniconda3." "$INFO"

    # Initialize conda only for the shells being used
    format_message "INFO" "Detecting current shell and initializing conda..." "$INFO"
    
    # Detectar el shell actual del usuario real
    USER_SHELL=$(getent passwd "$REAL_USER" | cut -d: -f7)
    format_message "INFO" "Detected shell for $REAL_USER: $USER_SHELL" "$INFO"
    
    # En lugar de usar conda init, añadimos manualmente la configuración necesaria
    CONDA_INIT_BLOCK="# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup=\"\$('$REAL_USER_HOME/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)\"
if [ \$? -eq 0 ]; then
    eval \"\$__conda_setup\"
else
    if [ -f \"$REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh\" ]; then
        . \"$REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh\"
    else
        export PATH=\"$REAL_USER_HOME/miniconda3/bin:\$PATH\"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<"

ZSH_CONDA_INIT_BLOCK="# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup=\"\$('$REAL_USER_HOME/miniconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)\"
if [ \$? -eq 0 ]; then
    eval \"\$__conda_setup\"
else
    if [ -f \"$REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh\" ]; then
        . \"$REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh\"
    else
        export PATH=\"$REAL_USER_HOME/miniconda3/bin:\$PATH\"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<"

# Inicializar conda para el shell actual del usuario
if [[ "$USER_SHELL" == *"zsh"* ]]; then
    format_message "INFO" "Adding conda initialization to zsh configuration..." "$INFO"
    if [[ -f "$REAL_USER_HOME/.zshrc" ]]; then
        # Primero eliminar cualquier inicialización existente de conda
        if grep -q "conda initialize" "$REAL_USER_HOME/.zshrc"; then
            # Guardar backup
            cp "$REAL_USER_HOME/.zshrc" "$REAL_USER_HOME/.zshrc.bak_$(date +%Y%m%d%H%M%S)"
            # Eliminar bloque existente
            sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$REAL_USER_HOME/.zshrc"
        fi
        # Añadir el nuevo bloque al final
        echo -e "\n$ZSH_CONDA_INIT_BLOCK" >> "$REAL_USER_HOME/.zshrc"
        chown "$REAL_USER:$(id -gn $REAL_USER)" "$REAL_USER_HOME/.zshrc"
        format_message "OK" "Conda initialization added to ~/.zshrc" "$OK"
    fi
else
    format_message "INFO" "Adding conda initialization to bash configuration..." "$INFO"
    if [[ -f "$REAL_USER_HOME/.bashrc" ]]; then
        # Primero eliminar cualquier inicialización existente de conda
        if grep -q "conda initialize" "$REAL_USER_HOME/.bashrc"; then
            # Guardar backup
            cp "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.bashrc.bak_$(date +%Y%m%d%H%M%S)"
            # Eliminar bloque existente
            sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$REAL_USER_HOME/.bashrc"
        fi
        # Añadir el nuevo bloque al final
        echo -e "\n$CONDA_INIT_BLOCK" >> "$REAL_USER_HOME/.bashrc"
        chown "$REAL_USER:$(id -gn $REAL_USER)" "$REAL_USER_HOME/.bashrc"
        format_message "OK" "Conda initialization added to ~/.bashrc" "$OK"
    fi
fi

# Opcionalmente preguntar por otros shells
if [[ "$USER_SHELL" != *"zsh"* ]] && [[ -f "$REAL_USER_HOME/.zshrc" ]]; then
    if ask_question "zsh configuration detected. Initialize conda for zsh as well? (y/N): " "n"; then
        # Primero eliminar cualquier inicialización existente de conda
        if grep -q "conda initialize" "$REAL_USER_HOME/.zshrc"; then
            # Guardar backup
            cp "$REAL_USER_HOME/.zshrc" "$REAL_USER_HOME/.zshrc.bak_$(date +%Y%m%d%H%M%S)"
            # Eliminar bloque existente
            sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$REAL_USER_HOME/.zshrc"
        fi
        # Añadir el nuevo bloque al final
        echo -e "\n$ZSH_CONDA_INIT_BLOCK" >> "$REAL_USER_HOME/.zshrc"
        chown "$REAL_USER:$(id -gn $REAL_USER)" "$REAL_USER_HOME/.zshrc"
        format_message "OK" "Conda initialization added to ~/.zshrc" "$OK"
    fi
elif [[ "$USER_SHELL" != *"bash"* ]] && [[ -f "$REAL_USER_HOME/.bashrc" ]]; then
    if ask_question "bash configuration detected. Initialize conda for bash as well? (y/N): " "n"; then
        # Primero eliminar cualquier inicialización existente de conda
        if grep -q "conda initialize" "$REAL_USER_HOME/.bashrc"; then
            # Guardar backup
            cp "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.bashrc.bak_$(date +%Y%m%d%H%M%S)"
            # Eliminar bloque existente
            sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$REAL_USER_HOME/.bashrc"
        fi
        # Añadir el nuevo bloque al final
        echo -e "\n$CONDA_INIT_BLOCK" >> "$REAL_USER_HOME/.bashrc"
        chown "$REAL_USER:$(id -gn $REAL_USER)" "$REAL_USER_HOME/.bashrc"
        format_message "OK" "Conda initialization added to ~/.bashrc" "$OK"
    fi
fi

    # Configure proxy if specified
    if [[ -n "$PROXY_URL" ]]; then
        format_message "INFO" "Configuring proxy for conda: $PROXY_URL" "$INFO"
        run_as_real_user "source $REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh && conda config --set proxy_servers.http $PROXY_URL && conda config --set proxy_servers.https $PROXY_URL"
        format_message "OK" "Proxy configured for conda." "$OK"
    fi

    echo
    format_message "INFO" "Conda has been initialized in the shell of user $REAL_USER." "$INFO"
    
    # Cargar conda temporalmente en la sesión actual para que los comandos siguientes funcionen
    export PATH="$REAL_USER_HOME/miniconda3/bin:$PATH"
    source "$REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
    
    format_message "INFO" "Conda loaded temporarily for current session." "$INFO"
    format_message "INFO" "For permanent access, restart your terminal or run: source ~/.bashrc" "$INFO"

    # Definir las variables de instalación explícitamente antes de usarlas
    CONDA_INSTALLED=true
    if ! command -v singularity &>/dev/null; then
        SINGULARITY_INSTALLED=false
    fi
    if ! command -v go &>/dev/null && [[ ! -d /usr/local/go ]]; then
        GO_INSTALLED=false
    fi

    # Installation verification
    if command -v conda &>/dev/null; then
        format_message "OK" "Conda installation verified:" "$OK"
        conda --version
    else
        format_message "WARNING" "Could not verify conda installation in this session." "$WARNING"
        format_message "TIP" "Restart your terminal and run 'conda --version' manually" "$TIP"
    fi

    # Actualizar la configuración EPIBAC con la nueva instalación de conda
    format_message "INFO" "Updating shell configuration with EPIBAC block containing conda..." "$INFO"
    update_epibac_shell_config "$CONDA_INSTALLED" "$SINGULARITY_INSTALLED" "$GO_INSTALLED"

    # Marcar explícitamente que nosotros lo instalamos
    CONDA_INSTALLED_BY_SCRIPT=true
}

###############################
# FUNCIONES: CONFIGURACIÓN EPIBAC #
###############################

# Constantes para las etiquetas de configuración
EPIBAC_TAG_START="# >>> EPIBAC_CONFIG >>>"
EPIBAC_TAG_END="# <<< EPIBAC_CONFIG <<<"
CONDA_TAG_START="# >>> conda initialize >>>"
CONDA_TAG_END="# <<< conda initialize <<<"

# Funciones para gestionar la configuración en .bashrc/.zshrc
# Función para verificar si existe un bloque de configuración
has_config_block() {
    local file="$1"
    local start_tag="$2"
    
    if [[ ! -f "$file" ]]; then
        return 1
    fi
    
    if grep -q "$start_tag" "$file"; then
        return 0  # El bloque existe
    else
        return 1  # El bloque no existe
    fi
}

# Función para extraer contenido de un bloque de configuración
extract_config_block() {
    local file="$1"
    local start_tag="$2"
    local end_tag="$3"
    
    if [[ ! -f "$file" ]] || ! has_config_block "$file" "$start_tag"; then
        echo ""  # Devolver cadena vacía si no existe el bloque
        return 1
    fi
    
    # Extraer contenido entre las etiquetas (incluyendo las etiquetas)
    sed -n "/$start_tag/,/$end_tag/p" "$file"
    return 0
}

# Función para eliminar un bloque de configuración
remove_config_block() {
    local file="$1"
    local start_tag="$2"
    local end_tag="$3"
    
    if [[ ! -f "$file" ]] || ! has_config_block "$file" "$start_tag"; then
        return 0  # No hay nada que eliminar
    fi
    
    # Crear backup antes de modificar
    cp "$file" "$file.bak_$(date +%Y%m%d%H%M%S)"
    
    # Eliminar el bloque completo (incluyendo las etiquetas)
    sed -i "/$start_tag/,/$end_tag/d" "$file"
    format_message "INFO" "Removed configuration block ($start_tag) from $file" "$INFO"
    return 0
}

# Función para añadir un bloque de configuración
add_config_block() {
    local file="$1"
    local start_tag="$2"
    local end_tag="$3"
    local content="$4"
    
    if [[ ! -f "$file" ]]; then
        format_message "ERROR" "File $file does not exist" "$ERROR"
        return 1
    fi
    
    # Crear backup antes de modificar
    cp "$file" "$file.bak_$(date +%Y%m%d%H%M%S)"
    
    # Si ya existe el bloque, eliminarlo primero
    if has_config_block "$file" "$start_tag"; then
        remove_config_block "$file" "$start_tag" "$end_tag"
    fi
    
    # Añadir el nuevo bloque al final del archivo
    {
        echo -e "\n$start_tag"
        echo -e "$content"
        echo -e "$end_tag"
    } >> "$file"
    
    format_message "INFO" "Added configuration block ($start_tag) to $file" "$INFO"
    return 0
}

# Función para encontrar la instalación de conda
find_conda_installation() {
    # Primero verificar si está disponible en PATH
    if command -v conda &>/dev/null; then
        # Obtener el directorio base de conda
        local conda_path=$(command -v conda)
        echo "$(dirname "$(dirname "$conda_path")")"
        return 0
    fi
    
    # Verificar ubicaciones comunes
    local common_locations=(
        "$REAL_USER_HOME/miniconda3"
        "$REAL_USER_HOME/anaconda3"
        "$REAL_USER_HOME/conda"
        "/opt/conda"
        "/opt/miniconda3"
        "/opt/anaconda3"
        "/opt/apps/conda"
        "/opt/apps/conda2"
        "/opt/apps/miniconda3"
        "/opt/apps/anaconda3"
    )
    
    for loc in "${common_locations[@]}"; do
        if [[ -f "$loc/bin/conda" ]]; then
            echo "$loc"
            return 0
        fi
    done
    
    # Buscar en los archivos .bashrc/.zshrc para buscar rutas de conda
    for rc_file in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
        if [[ -f "$rc_file" ]]; then
            # Buscar líneas tipo __conda_setup="$('/path/to/conda' 'shell.bash'...
            local conda_path=$(grep -o "'\([^']*\)/conda' 'shell" "$rc_file" 2>/dev/null | sed "s|' 'shell||g" | sed "s|'||g")
            if [[ -n "$conda_path" ]] && [[ -f "$conda_path" ]]; then
                echo "$(dirname "$(dirname "$conda_path")")"
                return 0
            fi
            
            # Buscar líneas tipo export PATH="/path/to/conda/bin:$PATH"
            local conda_bin_path=$(grep -o "PATH=\"\([^\"]*\)/conda[^/]*/bin" "$rc_file" 2>/dev/null | sed "s|PATH=\"||g" | sed "s|/bin||g")
            if [[ -n "$conda_bin_path" ]] && [[ -d "$conda_bin_path/bin" ]]; then
                echo "$conda_bin_path"
                return 0
            fi
        fi
    done
    
    # No se encontró conda
    return 1
}

# Función para generar el bloque de configuración EPIBAC
generate_epibac_config() {
    local install_conda=$1
    local install_singularity=$2
    local install_go=$3
    local conda_path=$4
    
    local config=""
    
    # Añadir sección de utilidad general (siempre se incluye)
    config+="# EPIBAC environment configuration\n"
    config+="# Helper function for clean PATH\n"
    config+="clean_path_dups() {\n"
    config+="    echo \"\$1\" | awk -v RS=':' '!a[\$0]++ {if (NR > 1) printf \":\"; printf \$0}'\n"
    config+="}\n"
    
    # Sección de Go SOLO si LO INSTALAMOS NOSOTROS
    if [[ "$GO_INSTALLED_BY_SCRIPT" == "true" ]]; then
        config+="\n# Go configuration\n"
        config+="export GOPATH=\"\$HOME/go\"\n"
        config+="export PATH=\"/usr/local/go/bin:\$PATH:\$GOPATH/bin\"\n"
    fi
    
    # Sección de Singularity SOLO si LO INSTALAMOS NOSOTROS
    if [[ "$SINGULARITY_INSTALLED_BY_SCRIPT" == "true" ]]; then
        config+="\n# Singularity configuration\n"
        config+="alias sing_version='singularity --version'\n"
        config+="alias sing_info='singularity --version && which singularity'\n"
    fi
    
    # Añadir alias de snake_env si fue solicitado
    if [[ "$SNAKE_ENV_ALIAS_REQUESTED" == "true" && -n "$SNAKE_ENV_PATH" ]]; then
        config+="\n# Python virtual environment configuration\n"
        config+="alias snake_env='source $SNAKE_ENV_PATH/bin/activate'\n"
    fi
    
    # Sección de Conda SOLO si LO INSTALAMOS NOSOTROS
    if [[ "$CONDA_INSTALLED_BY_SCRIPT" == "true" && -n "$conda_path" && -d "$conda_path" ]]; then
        config+="\n# Conda configuration\n"
        config+="conda_info() {\n"
        config+="    conda info --envs\n"
        config+="    echo -e \"\\nActive environment: \$(conda info | grep 'active env' | cut -d':' -f2 | xargs)\"\n"
        config+="}\n\n"
        config+="# Auto-deactivate conda on shell exit (better cleanup)\n"
        config+="trap \"conda deactivate 2>/dev/null || true\" EXIT\n\n"
        
        # Asegurarse de que conda esté en el PATH sin importar cómo se inicializa
        config+="# Ensure conda is always available in PATH\n"
        config+="[[ \":\$PATH:\" != *\":$conda_path/bin:\"* ]] && export PATH=\"$conda_path/bin:\$PATH\"\n\n"
        
        # Añadir la inicialización de conda - solo si existe el archivo
        if [[ -f "$conda_path/etc/profile.d/conda.sh" ]]; then
            config+="# Inicialización para cargar conda en la shell\n"
            config+="if [ -f \"$conda_path/etc/profile.d/conda.sh\" ]; then\n"
            config+="    . \"$conda_path/etc/profile.d/conda.sh\"\n"
            config+="fi\n"
        fi
    fi
    
    echo -e "$config"
}

# Función para actualizar o crear el bloque EPIBAC en los archivos de shell
update_epibac_shell_config() {
    local install_conda=$1
    local install_singularity=$2
    local install_go=$3
    
    # Asegurarse de que las variables estén definidas
    install_conda=${install_conda:-false}
    install_singularity=${install_singularity:-false}
    install_go=${install_go:-false}
    
    # Debug: Mostrar qué componentes se incluirán en el bloque
    format_message "INFO" "Updating EPIBAC configuration with:" "$INFO"
    [[ "$install_conda" == "true" ]] && format_message "INFO" "- Conda configuration" "$INFO"
    [[ "$install_singularity" == "true" ]] && format_message "INFO" "- Singularity configuration" "$INFO"
    [[ "$install_go" == "true" ]] && format_message "INFO" "- Go configuration" "$INFO"
    
    # Si no se instala nada nuevo, no tocar la configuración
    if [[ "$install_conda" != "true" && "$install_singularity" != "true" && "$install_go" != "true" ]]; then
        format_message "INFO" "No installations detected, shell configuration will remain unchanged" "$INFO"
        return 0
    fi
    
    # Encontrar la ruta de conda si está instalado
    local conda_path=""
    if [[ "$install_conda" == "true" ]]; then
        if ! conda_path=$(find_conda_installation); then
            format_message "WARNING" "Conda installation not found, using default path" "$WARNING"
            conda_path="$REAL_USER_HOME/miniconda3"
        fi
        format_message "INFO" "Using conda installation at: $conda_path" "$INFO"
    fi
    
    # Generar la configuración EPIBAC
    local epibac_config=$(generate_epibac_config "$install_conda" "$install_singularity" "$install_go" "$conda_path")
    
    # Actualizar los archivos de shell
    for shell_rc in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
        [[ -f "$shell_rc" ]] || continue
        
        format_message "INFO" "Updating EPIBAC configuration in: $(basename "$shell_rc")" "$INFO"
        
        # IMPORTANTE: NO eliminar el bloque conda initialize, es generado por conda init
        # Si hay un bloque EPIBAC, lo actualizamos pero preservamos el bloque conda
        if has_config_block "$shell_rc" "$EPIBAC_TAG_START"; then
            remove_config_block "$shell_rc" "$EPIBAC_TAG_START" "$EPIBAC_TAG_END"
        fi
        
        # Añadir/actualizar el bloque EPIBAC
        add_config_block "$shell_rc" "$EPIBAC_TAG_START" "$EPIBAC_TAG_END" "$epibac_config"
        
        # Si conda no está en el PATH del usuario en su .bashrc o .zshrc, asegurarse de que se agregue fuera del bloque EPIBAC
        if [[ "$install_conda" == "true" ]] && ! grep -q "miniconda3/bin" "$shell_rc"; then
            echo -e "\n# Ensure conda is in PATH (added by setup_conda_singularity.sh)" >> "$shell_rc"
            echo "export PATH=\"$REAL_USER_HOME/miniconda3/bin:\$PATH\"" >> "$shell_rc"
            format_message "INFO" "Added conda to PATH in $(basename "$shell_rc")" "$INFO"
        fi
    done
    
    format_message "OK" "EPIBAC configuration updated in shell files" "$OK"
    format_message "INFO" "Changes will take effect after restarting the terminal or running 'source ~/.bashrc'" "$INFO"
    return 0
}

# Función para limpiar la configuración EPIBAC
clean_epibac_shell_config() {
    for shell_rc in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
        [[ -f "$shell_rc" ]] || continue
        
        if has_config_block "$shell_rc" "$EPIBAC_TAG_START"; then
            format_message "INFO" "Removing EPIBAC configuration from $(basename "$shell_rc")" "$INFO"
            remove_config_block "$shell_rc" "$EPIBAC_TAG_START" "$EPIBAC_TAG_END"
        fi
    done
    
    format_message "INFO" "EPIBAC configuration removed from shell files" "$INFO"
    return 0
}

###############################
# FUNCTION: UNINSTALL EVERYTHING
###############################
uninstall_all() {
    echo
    format_message "INFO" "Starting uninstallation process..." "$INFO"
    
    # Check for conda installation
    CONDA_DETECTED=false
    if command -v conda &>/dev/null || [[ -d "$REAL_USER_HOME/miniconda3" ]] || [[ -d "$REAL_USER_HOME/anaconda3" ]]; then
        CONDA_DETECTED=true
        CONDA_PATH=""
        
        if command -v conda &>/dev/null; then
            CONDA_PATH=$(command -v conda)
        elif [[ -f "$REAL_USER_HOME/miniconda3/bin/conda" ]]; then
            CONDA_PATH="$REAL_USER_HOME/miniconda3/bin/conda"
        elif [[ -f "$REAL_USER_HOME/anaconda3/bin/conda" ]]; then
            CONDA_PATH="$REAL_USER_HOME/anaconda3/bin/conda"
        fi
        
        format_message "INFO" "Conda installation detected: $CONDA_PATH" "$INFO"
        if [[ -n "$CONDA_PATH" ]]; then
            # Try to get version info
            sudo -u "$REAL_USER" "$CONDA_PATH" --version || true
        fi
        
        # Ask to uninstall
        echo
        if ask_question "Do you want to uninstall conda and all environments? (y/N): " "n"; then
            remove_conda
            format_message "OK" "Conda has been uninstalled" "$OK"
        else
            format_message "INFO" "Conda will be kept" "$INFO"
        fi
    else
        format_message "INFO" "No conda installation detected" "$INFO"
    fi
    
    # Check for Go installation
    GO_DETECTED=false
    if command -v go &>/dev/null || [[ -d /usr/local/go ]]; then
        GO_DETECTED=true
        GO_PATH=""
        
        if command -v go &>/dev/null; then
            GO_PATH=$(command -v go)
            format_message "INFO" "Go detected in PATH: $GO_PATH" "$INFO"
        elif [[ -d /usr/local/go/bin ]]; then
            GO_PATH="/usr/local/go/bin/go"
            format_message "INFO" "Go detected in: /usr/local/go" "$INFO"
        fi
        
        if [[ -n "$GO_PATH" ]] && [[ -x "$GO_PATH" ]]; then
            "$GO_PATH" version || true
        else
            PATH="$PATH:/usr/local/go/bin" go version || true
        fi
        
        # Ask to uninstall
        echo
        if ask_question "Do you want to uninstall Go? (y/N): " "n"; then
            remove_go
            format_message "OK" "Go has been uninstalled" "$OK"
        else
            format_message "INFO" "Go will be kept" "$INFO"
        fi
    else
        format_message "INFO" "No Go installation detected" "$INFO"
    fi
    
    # Check for Singularity installation
    SINGULARITY_DETECTED=false
    if command -v singularity &>/dev/null; then
        SINGULARITY_DETECTED=true
        format_message "INFO" "Singularity detected: $(which singularity)" "$INFO"
        singularity --version || true
        
        # Ask to uninstall
        echo
        if ask_question "Do you want to uninstall Singularity? (y/N): " "n"; then
            remove_singularity
            format_message "OK" "Singularity has been uninstalled" "$OK"
        else
            format_message "INFO" "Singularity will be kept" "$INFO"
        fi
    else
        format_message "INFO" "No Singularity installation detected" "$INFO"
    fi
    
    # Clean up EPIBAC configuration from shell files
    echo
    if ask_question "Do you want to remove EPIBAC configuration from shell files (.bashrc/.zshrc)? (y/N): " "n"; then
        clean_epibac_shell_config
        format_message "OK" "EPIBAC configuration has been removed from shell files" "$OK"
    else
        format_message "INFO" "EPIBAC configuration will be kept in shell files" "$INFO"
    fi
    
    echo
    format_message "INFO" "Uninstallation process completed" "$INFO"
    format_message "TIP" "You may need to restart your terminal for all changes to take effect" "$TIP"
}

###############################
# MAIN FLOW LOGIC             #
###############################
case "$MODE" in
    singularity)
        format_message "INFO" "Mode: Singularity" "$INFO"

        # 1. Check if Go is installed and offer to uninstall it
        GO_AVAILABLE=false
        if command -v go &>/dev/null || [[ -d /usr/local/go ]] || [[ -d /opt/apps/go ]]; then
            GO_AVAILABLE=true
            GO_PATH=""

            # Determine the path to Go by searching in common locations
            if command -v go &>/dev/null; then
                GO_PATH=$(command -v go)
                format_message "INFO" "Go detected in PATH: $GO_PATH" "$INFO"
            elif [[ -d /usr/local/go/bin ]]; then
                GO_PATH="/usr/local/go/bin/go"
                format_message "INFO" "Go detected in: /usr/local/go" "$INFO"
            elif [[ -d /opt/apps/go/bin ]]; then
                GO_PATH="/opt/apps/go/bin/go"
                format_message "INFO" "Go detected in: /opt/apps/go" "$INFO"
            fi

            # Show version with expanded PATH
            if [[ -n "$GO_PATH" ]] && [[ -x "$GO_PATH" ]]; then
                "$GO_PATH" version || true
            else
                # Try with manually expanded PATH
                PATH="$PATH:/usr/local/go/bin:/opt/apps/go/bin" go version || true
            fi
            echo
            if ask_question "Uninstall previous Go? (y/N): " "n"; then
                remove_go
                GO_AVAILABLE=false
            else
                format_message "INFO" "Current Go installation kept." "$INFO"
            fi
        else
            format_message "INFO" "Go not detected in the system." "$INFO"
        fi
        echo
        # 2. Offer to install Go 1.24.1 if not available
        if ! $GO_AVAILABLE; then
            if ask_question "Install Go 1.24.1? (y/N): " "n"; then
                install_go
                GO_AVAILABLE=true
            else
                format_message "INFO" "Go will not be installed." "$INFO"
            fi
        fi

        # 3. Check if Singularity is present and offer to uninstall
        SINGULARITY_AVAILABLE=false
        if command -v singularity &>/dev/null; then
            SINGULARITY_AVAILABLE=true
            format_message "INFO" "Singularity detected: $(which singularity)" "$INFO"
            singularity --version || true
            echo
            if ask_question "Uninstall previous Singularity? (y/N): " "n"; then
                remove_singularity
                SINGULARITY_AVAILABLE=false
            else
                format_message "INFO" "Current Singularity installation kept." "$INFO"
            fi
        else
            format_message "INFO" "Singularity not detected in the system." "$INFO"
        fi
        echo
        # 4. Offer to install version 3.8.7
        if ! $SINGULARITY_AVAILABLE; then
            if ask_question "Install Singularity v3.8.7? (y/N): " "n"; then
                # Check if Go is available (required to compile Singularity)
                if ! $GO_AVAILABLE; then
                    format_message "WARNING" "Go is required to install Singularity." "$WARNING"
                    if ask_question "Install Go 1.24.1 first? (Y/n): " "y"; then
                        install_go
                        GO_AVAILABLE=true
                    else
                        format_message "INFO" "Singularity cannot be installed without Go." "$INFO"
                        # 'continue' cannot be used here since we are not in a loop
                        exit 1
                    fi
                fi

                install_singularity
                SINGULARITY_AVAILABLE=true
            else
                format_message "INFO" "Singularity will not be installed." "$INFO"
            fi
        fi
        echo
        # 5. Offer to install Python environment for Snakemake
        if ask_question "Create Python virtual environment with Snakemake? (y/N): " "n"; then
            # Captura el código de retorno explícitamente
            install_python_venv
            python_env_status=$?
            
            if [ $python_env_status -ne 0 ]; then
                format_message "WARNING" "There might have been issues with the Python environment installation" "$WARNING"
            fi
        else
            format_message "INFO" "Python environment will not be installed." "$INFO"
        fi
        ;;

    conda)
        format_message "INFO" "Mode: conda" "$INFO"

        # Variable to track if conda is available
        CONDA_AVAILABLE=false

        # 1. First, check if installation directories exist, even if they might be broken
        CONDA_DIRS_EXIST=false
        for CONDA_DIR in "$REAL_USER_HOME/miniconda3" "$REAL_USER_HOME/anaconda3" "$REAL_USER_HOME/conda"; do
            if [[ -d "$CONDA_DIR" ]]; then
                CONDA_DIRS_EXIST=true
                format_message "WARNING" "Found existing conda directory: $CONDA_DIR" "$WARNING"
            fi
        done

        # 2. Check if conda is installed and working correctly
        if command -v conda &>/dev/null || [[ -f "$REAL_USER_HOME/miniconda3/bin/conda" ]] || [[ -f "$REAL_USER_HOME/anaconda3/bin/conda" ]]; then
            # Determine the path to conda
            CONDA_PATH=""
            if command -v conda &>/dev/null; then
                CONDA_PATH=$(command -v conda)
            elif [[ -f "$REAL_USER_HOME/miniconda3/bin/conda" ]]; then
                CONDA_PATH="$REAL_USER_HOME/miniconda3/bin/conda"
            elif [[ -f "$REAL_USER_HOME/anaconda3/bin/conda" ]]; then
                CONDA_PATH="$REAL_USER_HOME/anaconda3/bin/conda"
            fi

            format_message "INFO" "Conda detected in: $CONDA_PATH" "$INFO"
            sudo -u "$REAL_USER" "$CONDA_PATH" --version || true
            echo
            
            # Simplificar la decisión a una pregunta directa
            if ask_question "Completely uninstall conda and reinstall? (Not necessary for configuring channels or environments) (y/N): " "n"; then
                format_message "WARNING" "Preparing to uninstall conda completely..." "$WARNING"
                if remove_conda; then
                    CONDA_AVAILABLE=false
                else
                    format_message "ERROR" "Conda uninstallation canceled." "$ERROR"
                    exit 1
                fi
            else
                format_message "INFO" "Keeping existing conda installation." "$INFO"
                CONDA_AVAILABLE=true
                CONDA_INSTALLED_BY_SCRIPT=true  # Tratar la instalación existente como nuestra
            fi
        elif $CONDA_DIRS_EXIST; then
            # 3. If directories were found but conda doesn't work, it's likely a broken installation
            format_message "ERROR" "Found conda directories but conda is not working correctly!" "$ERROR"
            format_message "WARNING" "This might be an incomplete or corrupted installation." "$WARNING"

            if ask_question "Remove existing conda directories before proceeding? (Y/n): " "y"; then
                remove_conda
                CONDA_AVAILABLE=false
            else
                format_message "ERROR" "Cannot continue installation with existing directories." "$ERROR"
                format_message "TIP" "Remove '$REAL_USER_HOME/miniconda3' manually or use this script to uninstall it." "$TIP"
                exit 1
            fi
        else
            format_message "INFO" "Conda not detected in the system." "$INFO"
        fi
        
        # 3. Offer to install conda only if it's not available
        if ! $CONDA_AVAILABLE; then
            echo
            if ask_question "Install Miniconda3: Conda 25.1.1 - Python 3.12.9? (y/N): " "n"; then
                install_conda
                CONDA_AVAILABLE=true
                CONDA_INSTALLED_BY_SCRIPT=true
            else
                format_message "INFO" "Miniconda will not be installed." "$INFO"
            fi
        fi

        # 4. Channel configuration and mamba now as a separate option
        if $CONDA_AVAILABLE; then
            echo
            if ask_question "Configure channels (conda-forge, bioconda) and install mamba? (y/N): " "n"; then
                format_message "INFO" "Configuring conda..." "$INFO"

                # Usar run_as_real_user en lugar de sudo -u
                run_as_real_user "source $REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh && conda config --remove-key channels || true && \
                    conda config --add channels defaults && \
                    conda config --add channels bioconda && \
                    conda config --add channels conda-forge && \
                    conda config --set channel_priority strict && \
                    conda update -n base conda -y && \
                    conda update --all -y && \
                    conda install -n base -c conda-forge mamba -y"

                format_message "OK" "Channels configured and mamba installed." "$OK"

                # Variable to track if mamba is available
                MAMBA_AVAILABLE=true
            else
                format_message "INFO" "Channels and mamba will not be configured." "$INFO"
                MAMBA_AVAILABLE=false
            fi
            
            # 5. Creation of the snake environment as a separate option
            echo
            if ask_question "Create the 'snake' environment with Snakemake? (y/N): " "n"; then
                format_message "INFO" "Creating 'snake' environment with Snakemake..." "$INFO"

                # If mamba is available, use it, otherwise use conda
                # In the section where you define INSTALL_CMD, add "conda" to the packages
                if $MAMBA_AVAILABLE; then
                    INSTALL_CMD="mamba create -n snake -y -c conda-forge bioconda::snakemake=9.1.1 \
                                                bioconda::snakemake-wrapper-utils=0.7.2 \
                                                pandas openpyxl git"
                else
                    INSTALL_CMD="conda create -n snake -y -c conda-forge bioconda::snakemake=9.1.1 \
                                                bioconda::snakemake-wrapper-utils=0.7.2 \
                                                pandas openpyxl git"
                fi

                # Execute the installation
                run_as_real_user "source $REAL_USER_HOME/miniconda3/etc/profile.d/conda.sh && $INSTALL_CMD && conda activate snake && snakemake --version"

                format_message "OK" "'snake' environment created successfully." "$OK"

                # Ensure that conda is always available
                echo
                format_message "INFO" "Ensuring conda availability in all environments..." "$INFO"
                for SHELL_RC in "$REAL_USER_HOME/.bashrc" "$REAL_USER_HOME/.zshrc"; do
                    [[ -f "$SHELL_RC" ]] || continue

                    # Check if a line already exists that adds miniconda3/bin to the PATH
                    if ! grep -q "miniconda3/bin" "$SHELL_RC"; then
                        echo "# Ensure access to conda in all environments" >>"$SHELL_RC"
                        echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >>"$SHELL_RC"
                        format_message "INFO" "miniconda3/bin has been added to PATH in $SHELL_RC" "$INFO"
                    fi
                done
                format_message "INFO" "To activate it use: conda activate snake" "$INFO"
            else
                format_message "INFO" "The snake environment will not be created." "$INFO"
            fi
        else
            format_message "WARNING" "Conda must be installed to configure channels, mamba, or create environments." "$WARNING"
        fi
        ;;
        
    uninstall)
        format_message "INFO" "Mode: Uninstall" "$INFO"
        uninstall_all
        ;;

    *)
        format_message "ERROR" "Invalid mode: $MODE" "$ERROR"
        echo "Usage: $0 {singularity|conda|uninstall}"
        exit 1
        ;;
esac

# Inicializar las variables de instalación según el estado actual
if command -v conda &>/dev/null; then
    CONDA_INSTALLED=true
else
    CONDA_INSTALLED=false
fi

if command -v singularity &>/dev/null; then
    SINGULARITY_INSTALLED=true
else
    SINGULARITY_INSTALLED=false
fi

if command -v go &>/dev/null || [[ -d /usr/local/go ]]; then
    GO_INSTALLED=true
else
    GO_INSTALLED=false
fi

format_message "INFO" "Checking if shell configuration needs to be updated..." "$INFO"
update_epibac_shell_config "$CONDA_INSTALLED" "$SINGULARITY_INSTALLED" "$GO_INSTALLED"

echo
echo -e "${BOLD}\033[32m====================================================================${COLOR_RESET}"
format_message "END" "Script 'setup_conda_singularity.sh' finished successfully." "$END"
echo -e "${BOLD}\033[32m====================================================================${COLOR_RESET}"
exit 0