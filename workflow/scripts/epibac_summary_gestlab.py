#!/usr/bin/env python3
"""
Script para generar el reporte GESTLAB con formato específico y salida en Excel y CSV.

Autor: GitHub Copilot
Fecha: 6 de junio de 2025
"""

import pandas as pd
import os
import sys
import re
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils.dataframe import dataframe_to_rows


def setup_logging():
    """Configurar logging básico."""
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def get_species_conversion_dict():
    """Diccionario para convertir nombres de especies abreviados a nombres completos."""
    return {
        'abaumannii': 'Acinetobacter baumannii',
        'acinetobacter_baumannii': 'Acinetobacter baumannii',
        'ecoli': 'Escherichia coli',
        'escherichia_coli': 'Escherichia coli',
        'kpneumoniae': 'Klebsiella pneumoniae',
        'klebsiella_pneumoniae': 'Klebsiella pneumoniae',
        'senterica': 'Salmonella enterica',
        'salmonella_enterica': 'Salmonella enterica',
        'paeruginosa': 'Pseudomonas aeruginosa',
        'pseudomonas_aeruginosa': 'Pseudomonas aeruginosa',
        'saureus': 'Staphylococcus aureus',
        'staphylococcus_aureus': 'Staphylococcus aureus',
        'efaecium': 'Enterococcus faecium',
        'enterococcus_faecium': 'Enterococcus faecium',
        'efaecalis': 'Enterococcus faecalis',
        'enterococcus_faecalis': 'Enterococcus faecalis',
        'sagalactiae': 'Streptococcus agalactiae',
        'streptococcus_agalactiae': 'Streptococcus agalactiae',
        'spneumoniae': 'Streptococcus pneumoniae',
        'streptococcus_pneumoniae': 'Streptococcus pneumoniae',
        'spyogenes': 'Streptococcus pyogenes',
        'streptococcus_pyogenes': 'Streptococcus pyogenes',
        'lmonocytogenes': 'Listeria monocytogenes',
        'listeria_monocytogenes': 'Listeria monocytogenes',
        'campylobacter_jejuni': 'Campylobacter jejuni',
        'cjejuni': 'Campylobacter jejuni',
        'campylobacter_coli': 'Campylobacter coli',
        'ccoli': 'Campylobacter coli'
    }


def convert_species_name(species_name):
    """Convertir nombre de especie abreviado a nombre completo."""
    if pd.isna(species_name) or not species_name or str(species_name).strip() == "" or str(species_name).strip() == "-":
        return "ND"
    
    species_dict = get_species_conversion_dict()
    species_lower = str(species_name).lower().strip()
    
    # Casos especiales
    if species_lower == "-" or species_lower == "":
        return "ND"
    
    # Buscar coincidencia exacta
    if species_lower in species_dict:
        return species_dict[species_lower]
    
    # Buscar coincidencia parcial
    for key, value in species_dict.items():
        if key in species_lower or species_lower in key:
            return value
    
    # Si no hay coincidencia, devolver "ND"
    return "ND"


def resumen_plasmidos(row, max_rep=8):
    """Generar resumen de plásmidos usando los datos de una fila del DataFrame."""
    try:
        partes = [
            f"{int(row['MOB_Suite_Num_Predicted_Plasmids_in_Sample'])} plasmids",
            f"size={int(row['MOB_Suite_Total_Size_Predicted_Plasmids_in_Sample'])}bp"
        ]
        if pd.notna(row['Unique_rep_type']):
            reps = [r.strip() for r in row['Unique_rep_type'].split(';') if r.strip()]
            partes.append(f"rep={', '.join(reps[:max_rep])}{'...' if len(reps) > max_rep else ''}")
        if pd.notna(row['Unique_predicted_mobility']):
            partes.append(f"mobility={row['Unique_predicted_mobility']}")
        return " | ".join(partes)
    except Exception as e:
        return "ND"


def generate_plasmid_summary(sample_id, plasmids_file, logger):
    """Generar resumen de plásmidos para una muestra específica."""
    try:
        logger.info(f"=== Procesando plásmidos para muestra: {sample_id} ===")
        
        # Verificar que el archivo existe
        if not os.path.exists(plasmids_file):
            logger.warning(f"Archivo de plásmidos no encontrado: {plasmids_file}")
            return "ND"
        
        # Leer archivo de plásmidos
        df_plasmids = pd.read_csv(plasmids_file, sep="\t")
        logger.info(f"Archivo de plásmidos cargado: {len(df_plasmids)} filas")
        logger.info(f"Columnas: {list(df_plasmids.columns)}")
        logger.info(f"Sample_IDs disponibles: {df_plasmids['Sample_ID'].tolist()}")
        
        # Buscar la muestra específica
        sample_data = df_plasmids[df_plasmids['Sample_ID'] == str(sample_id)]
        logger.info(f"Filas encontradas para '{sample_id}': {len(sample_data)}")
        
        if len(sample_data) == 0:
            logger.warning(f"No se encontraron datos de plásmidos para '{sample_id}'")
            return "ND"
        
        # Obtener la primera fila (debería ser única por muestra)
        row = sample_data.iloc[0]
        logger.info(f"Datos encontrados: {dict(row)}")
        
        # Verificar si tiene plásmidos
        num_plasmids = row.get('MOB_Suite_Num_Predicted_Plasmids_in_Sample', 0)
        if pd.isna(num_plasmids) or num_plasmids == 0:
            logger.info(f"Muestra '{sample_id}' no tiene plásmidos predichos")
            return "ND"
        
        # Usar la función resumen_plasmidos directamente
        summary = resumen_plasmidos(row)
        logger.info(f"Resumen generado para '{sample_id}': {summary}")
        return summary
        
    except Exception as e:
        logger.error(f"Error procesando plásmidos para '{sample_id}': {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        return "ND"


def process_path_nlsar(file_path):
    """
    Procesa una ruta de archivo para que comience desde NLSAR.
    
    Ejemplo:
    /home/asanzc/CabinaCS/NLSAR/deposit/CVA_EPIM/illumina/250425_EPIM199/fastq/EPI00579_S43_R1_001.fastq.gz
    -> NLSAR/deposit/CVA_EPIM/illumina/250425_EPIM199/fastq/EPI00579_S43_R1_001.fastq.gz
    """
    if pd.isna(file_path) or not file_path:
        return file_path
    
    # Buscar la posición de "NLSAR" en la ruta
    nlsar_pos = file_path.find("NLSAR")
    if nlsar_pos != -1:
        return file_path[nlsar_pos:]
    else:
        # Si no contiene NLSAR, devolver la ruta original
        return file_path


def determine_seq_method(row):
    """Determinar el tipo de secuenciación (ILLUMINA, NANOPORE o HYBRID)."""
    has_illumina = pd.notna(row.get("ILLUMINA_R1", "")) and row.get("ILLUMINA_R1", "") != ""
    has_nanopore = pd.notna(row.get("NANOPORE", "")) and row.get("NANOPORE", "") != ""

    if has_illumina and has_nanopore:
        return "HYBRID"
    elif has_illumina:
        return "ILLUMINA"
    elif has_nanopore:
        return "NANOPORE"
    else:
        return pd.NA


def build_path(row, config, tag_run):
    """Construir la ruta de FICHERO_LECTURAS_WGS."""
    if pd.isna(row.get("CARRERA", "")):
        return pd.NA
    
    parts = row["CARRERA"].split("_")
    if len(parts) < 2:
        return pd.NA
    
    # Obtener los parámetros específicos del modo GVA
    storage_cabinet = (
        config.get("mode_config", {})
        .get("gva", {})
        .get("storage_cabinet", "/home/asanzc/CABINAX/CabinaCS/NLSAR/deposit")
    )
    hosp = parts[1][:4]  # Primeros 4 caracteres del segundo segmento
    
    # Para simplificar, usamos la plataforma primaria para la ruta mostrada
    platform = "illumina" if row["OBS_MET_WGS"] in ["ILLUMINA", "HYBRID"] else "nanopore"
    
    # Construir ruta de destino base - usar forward slash (/) para Linux
    base_path = f"{storage_cabinet}/CVA_{hosp}/{platform}/{row['CARRERA']}/fastq"
    
    # Añadir el nombre del archivo si existe
    if platform == "illumina" and pd.notna(row.get("ILLUMINA_R1", "")) and row.get("ILLUMINA_R1", ""):
        # Usar solo el nombre del archivo, no la ruta completa
        filename = os.path.basename(row["ILLUMINA_R1"])
        full_path = f"{base_path}/{filename}"
    elif platform == "nanopore" and pd.notna(row.get("NANOPORE", "")) and row.get("NANOPORE", ""):
        # Usar solo el nombre del archivo, no la ruta completa
        filename = os.path.basename(row["NANOPORE"])
        full_path = f"{base_path}/{filename}"
    else:
        full_path = base_path
    
    # Procesar la ruta para que comience desde NLSAR
    return process_path_nlsar(full_path)


def clean_dataframe_for_excel(df, logger):
    """Limpiar DataFrame para compatibilidad con Excel."""
    try:
        # Crear una copia del dataframe para no modificar el original
        df_clean = df.copy()
        
        # Convertir valores pandas NA/NaN/None a "ND"
        df_clean = df_clean.fillna('ND')
        
        # Convertir cualquier objeto dtype a string para evitar problemas
        for col in df_clean.columns:
            if df_clean[col].dtype == 'object':
                df_clean[col] = df_clean[col].astype(str)
                # Reemplazar 'nan', '<NA>', 'None' por 'ND'
                df_clean[col] = df_clean[col].replace(['nan', '<NA>', 'None', 'NaN', ''], 'ND')
        
        logger.info("DataFrame limpiado para compatibilidad con Excel - campos vacíos rellenados con 'ND'")
        return df_clean
        
    except Exception as e:
        logger.error(f"Error limpiando DataFrame: {e}")
        return df


def create_excel_with_formatting(df, output_path, logger):
    """Crear archivo Excel con formato específico."""
    try:
        # Limpiar DataFrame antes de escribir a Excel
        df_clean = clean_dataframe_for_excel(df, logger)
        
        # Crear workbook y worksheet
        wb = Workbook()
        ws = wb.active
        ws.title = "GESTLAB_Report"
        
        # Agregar datos al worksheet
        for r in dataframe_to_rows(df_clean, index=False, header=True):
            ws.append(r)
        
        # Formatear la primera fila (encabezados)
        header_font = Font(bold=True, color="FFFFFF")
        header_fill = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")  # Azul clarito
        
        for cell in ws[1]:
            cell.font = header_font
            cell.fill = header_fill
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        
        # Ajustar ancho de columnas automáticamente
        for column in ws.columns:
            max_length = 0
            column_letter = column[0].column_letter
            
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            
            # Establecer ancho mínimo y máximo
            adjusted_width = min(max(max_length + 2, 10), 50)
            ws.column_dimensions[column_letter].width = adjusted_width
        
        # Aplicar wrap text a todas las celdas
        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = Alignment(wrap_text=True, vertical="top")
        
        # Guardar archivo
        wb.save(output_path)
        logger.info(f"Archivo Excel guardado exitosamente en: {output_path}")
        
    except Exception as e:
        logger.error(f"Error creando archivo Excel: {e}")
        raise


def main():
    """Función principal del script."""
    logger = setup_logging()
    
    # Obtener argumentos del snakemake
    try:
        validated_samples_info = snakemake.input.validated_samples_info
        results_file = snakemake.input.results
        plasmids_file = snakemake.input.plasmids
        output_csv = snakemake.output.gestlab_report_csv
        output_excel = snakemake.output.gestlab_report_xlsx
        tag_run = snakemake.params.tag_run
        outdir = snakemake.params.outdir
        config = snakemake.config
    except NameError:
        logger.error("Este script debe ejecutarse desde Snakemake")
        sys.exit(1)
    
    logger.info("=== Iniciando procesamiento de reporte GESTLAB ===")
    logger.info(f"TAG_RUN: {tag_run}")
    
    try:
        # Leer archivos de entrada
        logger.info("Leyendo archivos de entrada...")
        df_meta = pd.read_csv(validated_samples_info, sep=";")
        df_results = pd.read_csv(results_file, sep="\t")
        
        # Detectar si Sample coincide con id o con id2
        if df_results["Sample"].isin(df_meta["id"]).all():
            merge_col = "id"
        elif df_results["Sample"].isin(df_meta["id2"]).all():
            merge_col = "id2"
        else:
            raise ValueError(
                "No se puede encontrar coincidencia entre 'Sample' y 'id' o 'id2'."
            )
        
        logger.info(f"Usando columna '{merge_col}' para el merge")
        
        # Realizar merge
        df_results = df_results.rename(columns={"Sample": merge_col})
        df_merged = df_meta.merge(df_results, on=merge_col, how="left")
        
        # Eliminar columnas no deseadas
        columns_to_remove = ["SCOPE_core", "GENE_resfinder", "PHENO_resfinder", "MLST"]
        for col in columns_to_remove:
            if col in df_merged.columns:
                df_merged.drop(columns=col, inplace=True)
                logger.info(f"Eliminada columna: {col}")
        
        # Añadir CARRERA
        df_merged["CARRERA"] = tag_run
        
        # Renombrar columnas según el mapa GVA
        rename_map_gva = {
            "CODIGO_MUESTRA_ORIGEN": "id",
            "PETICION": "id2",
            "FECHA_TOMA_MUESTRA": "collection_date",
            "ESPECIE_SECUENCIA": "organism",
            "MOTIVO_WGS": "relevance",
            "ILLUMINA_R1": "illumina_r1",
            "ILLUMINA_R2": "illumina_r2",
            "NANOPORE": "nanopore",
            "ID_WGS": "Scheme_mlst",  # Cambiado de ID_WS a ID_WGS
            "ST_WGS": "ST",
            "MLST_WGS": "MLST",
            "R_Geno_WGS": "AMR",
            "PHENO_WGS": "PHENO_resfinder",
            "V_WGS": "VIRULENCE",
            "CONFIRMACION": "confirmation_note",
            "NUM_BROTE": "outbreak_id",
            "COMENTARIO_WGS": "comment",
        }
        
        df_merged.rename(
            columns={v: k for k, v in rename_map_gva.items()}, inplace=True
        )
        
        # Convertir nombres de especies (aplicar a ID_WGS y ESPECIE_SECUENCIA)
        if "ID_WGS" in df_merged.columns:
            df_merged["ID_WGS"] = df_merged["ID_WGS"].apply(convert_species_name)
            logger.info("Nombres de especies convertidos en ID_WGS")
        
        if "ESPECIE_SECUENCIA" in df_merged.columns:
            df_merged["ESPECIE_SECUENCIA"] = df_merged["ESPECIE_SECUENCIA"].apply(convert_species_name)
            logger.info("Nombres de especies convertidos en ESPECIE_SECUENCIA")
        
        # Generar resumen de plásmidos
        logger.info("Generando resumen de plásmidos...")
        logger.info(f"Columnas disponibles después del merge: {list(df_merged.columns)}")
        
        # Determinar qué columna usar para el ID de muestra
        # En el modo GVA, CODIGO_MUESTRA_ORIGEN se mapea a 'id' después del rename
        # y debe coincidir con Sample_ID en el archivo de plásmidos
        id_column = None
        if "CODIGO_MUESTRA_ORIGEN" in df_merged.columns:
            id_column = "CODIGO_MUESTRA_ORIGEN"
        elif "id" in df_merged.columns:
            id_column = "id"
        elif merge_col in df_merged.columns:
            id_column = merge_col
        
        if id_column:
            logger.info(f"Usando columna '{id_column}' para buscar plásmidos")
            logger.info(f"Valores en {id_column}: {df_merged[id_column].tolist()}")
            
            df_merged["PLASMIDOS_WGS"] = df_merged[id_column].apply(
                lambda sample_id: generate_plasmid_summary(sample_id, plasmids_file, logger)
            )
            logger.info(f"Resumen de plásmidos generado usando columna '{id_column}'")
        else:
            logger.warning("No se encontró columna de ID para generar resumen de plásmidos")
            df_merged["PLASMIDOS_WGS"] = "ND"
        
        # Añadir columnas extra
        df_merged["MO_EST_WGS"] = "BACTERIA"
        df_merged["TIPO_ANALISIS_WGS"] = "VERIFICACION"
        
        # Determinar el tipo de secuenciación
        df_merged["OBS_MET_WGS"] = df_merged.apply(determine_seq_method, axis=1)
        
        # Generar FICHERO_LECTURAS_WGS
        df_merged["FICHERO_LECTURAS_WGS"] = df_merged.apply(
            lambda row: build_path(row, config, tag_run), axis=1
        )
        
        # Añadir las nuevas columnas requeridas
        df_merged["ENVIAR_REDLABRA"] = "ND"
        df_merged["VERIFICACION_WGS"] = "ND"
        
        # Generar ruta para FICHERO_ANALISIS_WGS
        analysis_file_path = f"epibac/output/{tag_run}/report/{tag_run}_EPIBAC.tsv"
        df_merged["FICHERO_ANALISIS_WGS"] = analysis_file_path
        
        # Eliminar las columnas especificadas del archivo final
        columns_to_drop_final = ["ILLUMINA_R1", "ILLUMINA_R2", "NANOPORE", "dorado_model"]
        for col in columns_to_drop_final:
            if col in df_merged.columns:
                df_merged.drop(columns=col, inplace=True)
                logger.info(f"Eliminada columna del resultado final: {col}")
        
        # Rellenar campos vacíos con "ND" para ambos archivos de salida
        df_merged = df_merged.fillna('ND')
        for col in df_merged.columns:
            if df_merged[col].dtype == 'object':
                df_merged[col] = df_merged[col].astype(str)
                df_merged[col] = df_merged[col].replace(['nan', '<NA>', 'None', 'NaN', ''], 'ND')
        
        logger.info("Campos vacíos rellenados con 'ND' en el DataFrame final")
        
        # Guardar archivo CSV
        df_merged.to_csv(output_csv, sep=";", index=False)
        logger.info(f"Archivo CSV guardado en: {output_csv}")
        
        # Crear archivo Excel con formato
        create_excel_with_formatting(df_merged, output_excel, logger)
        
        # Estadísticas
        logger.info("=== Estadísticas finales ===")
        logger.info(f"Total de muestras procesadas: {len(df_merged)}")
        logger.info(f"Columnas en el archivo final: {len(df_merged.columns)}")
        logger.info(f"Columnas: {list(df_merged.columns)}")
        
        logger.info("=== Procesamiento completado exitosamente ===")
        
    except Exception as e:
        logger.error(f"Error durante el procesamiento: {e}")
        raise


if __name__ == "__main__":
    main()
