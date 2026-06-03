import os
from pathlib import Path

def collect_tfa_sequences(base_dir: str, output_file: str = "mlst.fa"):
    base_path = Path(base_dir)
    with open(output_file, "w") as outfile:
        for species_dir in sorted(base_path.iterdir()):
            if species_dir.is_dir():
                species_name = species_dir.name
                for tfa_file in sorted(species_dir.glob("*.tfa")):
                    gene_name = tfa_file.stem
                    with open(tfa_file, "r") as infile:
                        for line in infile:
                            if line.startswith(">"):
                                header = line.strip().lstrip(">")
                                new_header = f">{species_name}.{header}"
                                outfile.write(new_header + "\n")
                            else:
                                outfile.write(line)
    print(f"Fichero '{output_file}' generado correctamente.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combina archivos .tfa de múltiples especies en un único mlst.fa")
    parser.add_argument("input_dir", help="Directorio raíz que contiene carpetas por especie con archivos .tfa")
    parser.add_argument("-o", "--output", default="mlst.fa", help="Nombre del fichero de salida (por defecto: mlst.fa)")

    args = parser.parse_args()
    collect_tfa_sequences(args.input_dir, args.output)
