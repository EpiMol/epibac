#!/usr/bin/env python3
"""
Script unificado para análisis completo de plásmidos y elementos móviles genéticos (MGE)
Combina funcionalidades de Platon y MOB-suite para generar resúmenes por muestra.

Autor: GitHub Copilot
Fecha: 6 de junio de 2025
"""

import pandas as pd
import numpy as np
import json
import glob
import os
import sys
import logging
import re
from pathlib import Path
from typing import Dict, List, Any, Optional


def setup_logging(log_file: Optional[str] = None) -> logging.Logger:
    """Configurar logging para el script."""
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)


def safe_get(data_dict: Dict[str, Any], key_list: List[str], default: Any = "N/A") -> Any:
    """
    Intenta obtener un valor de un diccionario usando una lista de claves posibles.
    Devuelve el valor de la primera clave encontrada o el valor por defecto.
    """
    if not isinstance(data_dict, dict):
        return default
    
    for key in key_list:
        if key in data_dict and data_dict[key] is not None:
            return data_dict[key]
    
    return default


def format_hit_details(hit_object: Dict[str, Any]) -> str:
    """
    Formatea una cadena de detalles para un hit individual,
    incluyendo dinámicamente la información disponible.
    """
    if not isinstance(hit_object, dict):
        return "N/A"
    
    details = []
    
    # Claves comunes y su orden preferido de aparición
    # Prioridad para nombre/tipo
    primary_id = safe_get(hit_object, ["gene", "type", "name", "hmm-id"], None)
    if primary_id:
        details.append(str(primary_id))

    accession = safe_get(hit_object, ["accession"], None)
    if accession and accession not in str(primary_id):
        details.append(f"Acc:{accession}")
    
    sub_details = []
    
    # Cobertura e Identidad (si existen y son numéricos)
    try:
        coverage = float(safe_get(hit_object, ["coverage"], None))
        if coverage is not None:
            sub_details.append(f"Cov:{coverage:.1f}%")
    except (ValueError, TypeError, AttributeError):
        pass
    
    try:
        identity = float(safe_get(hit_object, ["identity"], None))
        if identity is not None:
            sub_details.append(f"Id:{identity:.1f}%")
    except (ValueError, TypeError, AttributeError):
        pass

    # Bitscore y E-value (si existen)
    bitscore = safe_get(hit_object, ["bitscore"], None)
    if bitscore is not None:
        sub_details.append(f"Bits:{bitscore}")
    
    evalue = safe_get(hit_object, ["evalue"], None)
    if evalue is not None:
        sub_details.append(f"E:{evalue}")

    if sub_details:
        details.append(f"[{','.join(sub_details)}]")
        
    return " ".join(details) if details else "N/A"


def extract_array_hits_info(hits_array: List[Dict[str, Any]]) -> str:
    """
    Extrae y formatea información de un array de 'hits' de Platon,
    creando una cadena única separada por punto y coma.
    """
    if not hits_array or not isinstance(hits_array, list):
        return "N/A"
    
    formatted_hits = [format_hit_details(hit) for hit in hits_array if isinstance(hit, dict)]
    valid_formatted_hits = [fh for fh in formatted_hits if fh and fh != "N/A"]
            
    return "; ".join(valid_formatted_hits) if valid_formatted_hits else "N/A"


def extract_plasmid_db_hits_info(hits_array: List[Dict[str, Any]]) -> str:
    """Extrae y formatea información específica del array 'plasmid_hits'."""
    if not hits_array or not isinstance(hits_array, list):
        return "N/A"

    formatted_hits = []
    for hit in hits_array:
        if isinstance(hit, dict):
            name = safe_get(hit, ["name"], "unknown")
            coverage = safe_get(hit, ["coverage"], None)
            identity = safe_get(hit, ["identity"], None)
            
            hit_str = str(name)
            if coverage is not None and identity is not None:
                hit_str += f" [Cov:{coverage:.1f}%,Id:{identity:.1f}%]"
            
            formatted_hits.append(hit_str)
            
    return "; ".join(formatted_hits) if formatted_hits else "N/A"


def extract_platon_data_from_json(json_filepath: str) -> List[Dict[str, Any]]:
    """
    Extrae datos completos de un archivo JSON de Platon para una muestra.
    
    Args:
        json_filepath (str): Ruta al archivo JSON de Platon
        
    Returns:
        List[Dict]: Lista de diccionarios con información de cada contig
    """
    sample_id = os.path.basename(json_filepath).replace('.json', '')
    contigs_data = []
    
    try:
        with open(json_filepath, 'r') as f:
            platon_data = json.load(f)
        
        if 'records' not in platon_data:
            logging.warning(f"No se encontró 'records' en {json_filepath}")
            return contigs_data
        
        for record in platon_data['records']:
            contig_info = {
                'Sample_ID': sample_id,
                'Contig_ID': safe_get(record, ['id'], 'unknown'),
                'Length': safe_get(record, ['length'], 0),
                'Is_Circular': 'yes' if safe_get(record, ['circular'], False) else 'no',
                'Sequence_Coverage': safe_get(record, ['coverage'], 'N/A'),
                'Protein_Score': safe_get(record, ['rds'], 'N/A'),
                'Inc_Types_Info': safe_get(record, ['inc_types'], []),
                'Replication_Hits_Info': extract_array_hits_info(safe_get(record, ['replication'], [])),
                'Mobilization_Hits_Info': extract_array_hits_info(safe_get(record, ['mobilization'], [])),
                'OriT_Hits_Info': extract_array_hits_info(safe_get(record, ['orit'], [])),
                'Conjugation_Hits_Info': extract_array_hits_info(safe_get(record, ['conjugation'], [])),
                'AMR_Platon_Hits_Info': extract_array_hits_info(safe_get(record, ['amr'], [])),
                'PlasmidDB_Hits_Info': extract_plasmid_db_hits_info(safe_get(record, ['plasmid_hits'], [])),
                'rRNA_Hits_Info': extract_array_hits_info(safe_get(record, ['rrna'], []))
            }
            
            # Convertir inc_types de lista a string
            if isinstance(contig_info['Inc_Types_Info'], list):
                contig_info['Inc_Types_Info'] = "; ".join(contig_info['Inc_Types_Info']) if contig_info['Inc_Types_Info'] else "N/A"
            
            contigs_data.append(contig_info)
            
    except Exception as e:
        logging.error(f"Error procesando {json_filepath}: {e}")
    
    return contigs_data


def safe_split_and_unique(series: pd.Series, separator: str = '; ') -> str:
    """
    Toma una serie de pandas, donde cada celda puede contener una cadena
    con múltiples valores separados por 'separator'.
    Devuelve una cadena única de valores únicos, ordenados y separados por 'separator'.
    Maneja 'N/A' y valores vacíos.
    """
    all_items = set()
    for item_list_str in series.dropna():
        if isinstance(item_list_str, str) and item_list_str != "N/A":
            items = item_list_str.split(separator)
            for item in items:
                cleaned_item = item.strip()
                if cleaned_item and cleaned_item != "N/A":
                    all_items.add(cleaned_item)
    
    return separator.join(sorted(list(all_items))) if all_items else "N/A"


def extract_primary_ids(text_series: pd.Series, item_separator: str = '; ') -> str:
    """
    Toma una serie de Pandas con strings formateados (ej. 'ID1 [detalles]; ID2 [detalles]').
    Extrae solo los identificadores principales (la parte antes del primer '[' o '(').
    Devuelve una cadena única de IDs únicos, ordenados y separados por 'item_separator'.
    """
    if not isinstance(text_series, pd.Series):
        return "N/A"
    
    primary_ids = set()
    for item_str in text_series.dropna():
        if isinstance(item_str, str) and item_str != "N/A":
            items = item_str.split(item_separator)
            for item in items:
                item_stripped = item.strip()
                if item_stripped and item_stripped != "N/A":
                    # Extraer parte principal antes de '[', '(', o ' ['
                    primary_part = re.split(r'[\[\(]', item_stripped)[0].strip()
                    if primary_part:
                        primary_ids.add(primary_part)
                        
    return item_separator.join(sorted(list(primary_ids))) if primary_ids else "N/A"


def filter_platon_plasmids(df_detailed: pd.DataFrame, min_length: int = 2000) -> pd.DataFrame:
    """
    Filtra contigs de Platon para identificar plásmidos potenciales.
    
    Args:
        df_detailed (pd.DataFrame): DataFrame con datos detallados de Platon
        min_length (int): Longitud mínima para considerar un contig
        
    Returns:
        pd.DataFrame: DataFrame filtrado con plásmidos potenciales
    """
    if df_detailed.empty:
        return df_detailed
    
    # Convertir 'Length' a numérico
    df_detailed['Length'] = pd.to_numeric(df_detailed['Length'], errors='coerce')
    df_detailed = df_detailed.dropna(subset=['Length'])

    # Criterio para considerar un contig como plasmídico
    is_plasmid_criteria = (
        (df_detailed['Is_Circular'] == 'yes') |
        (df_detailed['Inc_Types_Info'].notna() & (df_detailed['Inc_Types_Info'] != 'N/A')) |
        (df_detailed['Replication_Hits_Info'].notna() & (df_detailed['Replication_Hits_Info'] != 'N/A')) |
        (df_detailed['PlasmidDB_Hits_Info'].notna() & (df_detailed['PlasmidDB_Hits_Info'] != 'N/A'))
    )
    
    df_filtered = df_detailed[
        (df_detailed['Length'] >= min_length) & is_plasmid_criteria
    ].copy()

    return df_filtered


def summarize_platon_per_sample(df_filtered_plasmids: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Crea resumen de plásmidos de Platon por muestra.
    
    Args:
        df_filtered_plasmids (pd.DataFrame): DataFrame con plásmidos filtrados
        
    Returns:
        List[Dict]: Lista con resúmenes por muestra
    """
    summary_list = []

    for sample_id, group in df_filtered_plasmids.groupby('Sample_ID'):
        summary = {
            'Sample_ID': sample_id,
            'Num_Filtered_Plasmid_Contigs': len(group),
            'Total_Filtered_Plasmid_Length': int(group['Length'].sum()),
            'Circular_Filtered_Plasmids_Count': len(group[group['Is_Circular'] == 'yes']),
            'Unique_Inc_Types_on_Plasmids': safe_split_and_unique(group['Inc_Types_Info']),
            'Unique_Replication_Hits_on_Plasmids': extract_primary_ids(group['Replication_Hits_Info']),
            'Unique_Mobilization_Hits_on_Plasmids': extract_primary_ids(group['Mobilization_Hits_Info']),
            'Unique_PlasmidDB_Hits_on_Plasmids': extract_primary_ids(group['PlasmidDB_Hits_Info'])
        }
        summary_list.append(summary)

    return summary_list


def parse_mob_suite_file(filepath: str) -> pd.DataFrame:
    """
    Parsea un archivo de MOB-suite manejando múltiples headers.
    
    Args:
        filepath (str): Ruta al archivo de MOB-suite
        
    Returns:
        pd.DataFrame: DataFrame con datos de MOB-suite
    """
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return pd.DataFrame()
    
    try:
        # Leer todas las líneas del archivo
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Encontrar todas las líneas de header
        header_indices = []
        for i, line in enumerate(lines):
            if line.startswith('sample_id\t') or line.startswith('#FILE\t'):
                header_indices.append(i)
        
        if not header_indices:
            # Intentar lectura normal si no hay headers especiales
            return pd.read_csv(filepath, sep='\t')
        
        # Usar el primer header encontrado y filtrar líneas de datos
        header_line = lines[header_indices[0]]
        data_lines = []
        
        for i, line in enumerate(lines):
            # Saltear líneas de header y líneas vacías
            if i not in header_indices and line.strip():
                # Solo añadir si la línea tiene el mismo número de columnas que el header
                if len(line.split('\t')) == len(header_line.split('\t')):
                    data_lines.append(line)
        
        if not data_lines:
            return pd.DataFrame()
        
        # Crear DataFrame a partir de las líneas filtradas
        from io import StringIO
        csv_string = header_line + ''.join(data_lines)
        df = pd.read_csv(StringIO(csv_string), sep='\t')
        
        return df
        
    except Exception as e:
        logging.warning(f"Error parseando {filepath}: {e}")
        return pd.DataFrame()


def safe_join_unique_sorted(series: pd.Series, separator: str = '; ') -> str:
    """
    Toma una serie de Pandas y devuelve una cadena única de valores únicos,
    ordenados y separados por 'separator'. Maneja 'N/A', '-', y valores vacíos/nulos.
    """
    all_items = set()
    for item_val in series.dropna():
        s_item_val = str(item_val).strip()
        if s_item_val and s_item_val != 'N/A' and s_item_val != '-':
            # Split por comas si hay múltiples valores en una celda
            sub_items = s_item_val.split(',')
            for sub_item in sub_items:
                cleaned_sub_item = sub_item.strip()
                if cleaned_sub_item and cleaned_sub_item != 'N/A' and cleaned_sub_item != '-':
                    all_items.add(cleaned_sub_item)
                    
    return separator.join(sorted(list(all_items))) if all_items else "N/A"


def summarize_mob_suite_per_sample(df_mob_all: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Crea resumen de MOB-suite por muestra.
    
    Args:
        df_mob_all (pd.DataFrame): DataFrame con todos los datos de MOB-suite
        
    Returns:
        List[Dict]: Lista con resúmenes por muestra
    """
    if df_mob_all.empty:
        return []
    
    # Columnas que queremos agregar/colapsar
    columns_to_aggregate = {
        'rep_type(s)': 'Unique_rep_type',
        'rep_type_accession(s)': 'Unique_rep_type_accession',
        'relaxase_type(s)': 'Unique_relaxase_type',
        'relaxase_type_accession(s)': 'Unique_relaxase_type_accession',
        'mpf_type': 'Unique_mpf_type',
        'mpf_type_accession(s)': 'Unique_mpf_type_accession',
        'orit_type(s)': 'Unique_orit_type',
        'orit_accession(s)': 'Unique_orit_accession',
        'predicted_mobility': 'Unique_predicted_mobility',
        'mash_neighbor_identification': 'Unique_mash_neighbor_identification'
    }

    summary_list = []

    # Determinar la columna de ID de muestra
    sample_id_col = None
    possible_cols = ['Original_Sample_ID', 'Sample_ID', 'sample_id']
    for col in possible_cols:
        if col in df_mob_all.columns:
            sample_id_col = col
            break
    
    if not sample_id_col:
        logging.error("No se encontró columna de ID de muestra en datos de MOB-suite")
        return []

    # Agrupar por la muestra
    for sample_id, group in df_mob_all.groupby(sample_id_col):
        summary = {
            'Sample_ID': sample_id,
            'MOB_Suite_Num_Predicted_Plasmids_in_Sample': len(group),
            'MOB_Suite_Total_Size_Predicted_Plasmids_in_Sample': int(group.get('size', pd.Series([0])).sum())
        }
        
        # Agregar columnas específicas
        for orig_col, new_col in columns_to_aggregate.items():
            if orig_col in group.columns:
                summary[new_col] = safe_join_unique_sorted(group[orig_col])
            else:
                summary[new_col] = "N/A"
        
        summary_list.append(summary)

    return summary_list


def process_platon_data(platon_base_path: str, min_length: int = 2000) -> pd.DataFrame:
    """
    Procesa todos los archivos JSON de Platon y crea resumen por muestra.
    
    Args:
        platon_base_path (str): Directorio base con archivos JSON de Platon
        min_length (int): Longitud mínima para filtrar plásmidos
        
    Returns:
        pd.DataFrame: DataFrame con resumen de Platon por muestra
    """
    logging.info(f"Procesando datos de Platon desde: {platon_base_path}")
    
    # Buscar archivos JSON de Platon
    json_pattern = os.path.join(platon_base_path, "**", "EPI*.json")
    json_files = glob.glob(json_pattern, recursive=True)
    
    if not json_files:
        logging.warning(f"No se encontraron archivos JSON de Platon en {platon_base_path}")
        return pd.DataFrame()

    logging.info(f"Se encontraron {len(json_files)} archivos JSON de Platon")
    
    # Extraer datos de todos los archivos JSON
    all_contigs_data = []
    for json_filepath in sorted(json_files):
        sample_contigs = extract_platon_data_from_json(json_filepath)
        all_contigs_data.extend(sample_contigs)

    if not all_contigs_data:
        logging.warning("No se pudieron extraer datos de ningún archivo JSON de Platon")
        return pd.DataFrame()

    # Crear DataFrame detallado
    df_detailed = pd.DataFrame(all_contigs_data)
    logging.info(f"Total de contigs extraídos de Platon: {len(df_detailed)}")

    # Filtrar plásmidos potenciales
    df_filtered_plasmids = filter_platon_plasmids(df_detailed, min_length)
    logging.info(f"Contigs considerados plasmídicos después de filtrar: {len(df_filtered_plasmids)}")

    if df_filtered_plasmids.empty:
        logging.warning("No se encontraron plásmidos después de aplicar filtros")
        return pd.DataFrame()

    # Crear resumen por muestra
    summary_list = summarize_platon_per_sample(df_filtered_plasmids)
    df_summary = pd.DataFrame(summary_list)
    
    return df_summary


def process_mob_suite_data(mob_suite_base_path: str) -> pd.DataFrame:
    """
    Procesa todos los archivos de MOB-suite y crea resumen por muestra.
    
    Args:
        mob_suite_base_path (str): Directorio base con archivos de MOB-suite
        
    Returns:
        pd.DataFrame: DataFrame con resumen de MOB-suite por muestra
    """
    logging.info(f"Procesando datos de MOB-suite desde: {mob_suite_base_path}")
    
    # Buscar archivos de MOB-suite (mobtyper_results.txt)
    mob_pattern = os.path.join(mob_suite_base_path, "**", "mobtyper_results.txt")
    mob_files = glob.glob(mob_pattern, recursive=True)
    
    if not mob_files:
        logging.warning(f"No se encontraron archivos de MOB-suite en {mob_suite_base_path}")
        return pd.DataFrame()

    logging.info(f"Se encontraron {len(mob_files)} archivos de MOB-suite")
    
    # Procesar todos los archivos de MOB-suite
    all_plasmids_list = []
    for filepath in sorted(mob_files):
        # Extraer sample_id del path
        path_parts = filepath.split(os.sep)
        sample_id = None
        for part in reversed(path_parts):
            if part.startswith('EPI'):
                sample_id = part
                break
        
        if not sample_id:
            logging.warning(f"No se pudo extraer sample_id de {filepath}")
            continue
            
        logging.info(f"Procesando MOB-suite para muestra: {sample_id}")
        
        df_mob_sample = parse_mob_suite_file(filepath)
        if not df_mob_sample.empty:
            df_mob_sample.insert(0, 'Original_Sample_ID', sample_id)
            all_plasmids_list.append(df_mob_sample)

    if not all_plasmids_list:
        logging.warning("No se pudieron procesar datos de ningún archivo de MOB-suite")
        return pd.DataFrame()

    # Combinar todos los DataFrames
    df_mob_all = pd.concat(all_plasmids_list, ignore_index=True)
    logging.info(f"Total de plásmidos predichos por MOB-suite: {len(df_mob_all)}")

    # Crear resumen por muestra
    summary_list = summarize_mob_suite_per_sample(df_mob_all)
    df_summary = pd.DataFrame(summary_list)
    
    return df_summary


def merge_platon_mob_suite(df_platon: pd.DataFrame, df_mob_suite: pd.DataFrame) -> pd.DataFrame:
    """
    Combina los resúmenes de Platon y MOB-suite.
    
    Args:
        df_platon (pd.DataFrame): Resumen de Platon por muestra
        df_mob_suite (pd.DataFrame): Resumen de MOB-suite por muestra
        
    Returns:
        pd.DataFrame: DataFrame combinado
    """
    logging.info("Combinando datos de Platon y MOB-suite")
    
    if df_platon.empty and df_mob_suite.empty:
        logging.warning("Ambos DataFrames están vacíos")
        return pd.DataFrame()
    elif df_platon.empty:
        logging.info("Solo se utilizarán datos de MOB-suite")
        return df_mob_suite
    elif df_mob_suite.empty:
        logging.info("Solo se utilizarán datos de Platon")
        return df_platon
    else:
        logging.info("Realizando merge de datos de Platon y MOB-suite")
        df_combined = pd.merge(df_platon, df_mob_suite, on='Sample_ID', how='outer')
        
        # Llenar NaNs con N/A
        for col in df_combined.columns:
            if df_combined[col].dtype == 'object':
                df_combined[col] = df_combined[col].fillna('N/A')
            else:
                df_combined[col] = df_combined[col].fillna(0)
        
        return df_combined


def simplify_and_format_for_excel(df: pd.DataFrame) -> pd.DataFrame:
    """
    Simplifica las columnas para generar un formato más limpio para Excel.
    
    Args:
        df (pd.DataFrame): DataFrame con datos combinados
        
    Returns:
        pd.DataFrame: DataFrame simplificado para Excel
    """
    # Columnas que queremos simplificar (extraer solo IDs principales)
    columns_to_simplify = [
        "Unique_Inc_Types_on_Plasmids",
        "Unique_Replication_Hits_on_Plasmids",
        "Unique_Mobilization_Hits_on_Plasmids",
        "Unique_predicted_mobility"
    ]
    
    df_excel = df.copy()
    
    for col in columns_to_simplify:
        if col in df_excel.columns:
            # Simplificar extrayendo solo la parte principal antes de '[' o '('
            df_excel[col] = df_excel[col].apply(lambda x: simplify_cell_content(x))
    
    return df_excel


def simplify_cell_content(cell_content: Any, item_separator: str = '; ') -> str:
    """
    Simplifica el contenido de una celda extrayendo solo los identificadores principales.
    
    Args:
        cell_content (Any): Contenido de la celda
        item_separator (str): Separador entre elementos
        
    Returns:
        str: Contenido simplificado
    """
    if pd.isna(cell_content) or str(cell_content).strip() == "N/A" or str(cell_content).strip() == "":
        return "N/A"

    items = str(cell_content).split(item_separator)
    primary_ids = set()
    
    for item in items:
        item_stripped = item.strip()
        if item_stripped and item_stripped != "N/A":
            # Extraer la parte principal antes de '[', '(', o ' ['
            primary_part = re.split(r'[\[\(]', item_stripped)[0].strip()
            if primary_part:
                primary_ids.add(primary_part)
                
    return "\n".join(sorted(list(primary_ids))) if primary_ids else "N/A"


def main():
    """Función principal del script."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Análisis unificado de plásmidos Platon + MOB-suite")
    parser.add_argument("--platon_path", required=True, help="Directorio base con archivos JSON de Platon")
    parser.add_argument("--mobsuite_path", required=True, help="Directorio base con archivos de MOB-suite")
    parser.add_argument("--output_dir", default=".", help="Directorio de salida (default: directorio actual)")
    parser.add_argument("--min_length", type=int, default=2000, help="Longitud mínima para filtrar plásmidos (default: 2000)")
    parser.add_argument("--log_file", help="Archivo de log (opcional)")
    
    try:
        args = parser.parse_args()
    except SystemExit as e:
        print(f"Error en argumentos: {e}")
        return
    
    # Configurar logging
    logger = setup_logging(args.log_file)
    logger.info("=== Iniciando análisis unificado de plásmidos ===")
    logger.info(f"Directorio Platon: {args.platon_path}")
    logger.info(f"Directorio MOB-suite: {args.mobsuite_path}")
    logger.info(f"Directorio salida: {args.output_dir}")
    logger.info(f"Longitud mínima: {args.min_length}")
    
    try:
        # Crear directorio de salida si no existe
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Procesar datos de Platon
        logger.info("=== Procesando datos de Platon ===")
        df_platon = process_platon_data(args.platon_path, args.min_length)
        
        # Procesar datos de MOB-suite
        logger.info("=== Procesando datos de MOB-suite ===")
        df_mob_suite = process_mob_suite_data(args.mobsuite_path)
        
        # Combinar datos
        logger.info("=== Combinando datos ===")
        df_combined = merge_platon_mob_suite(df_platon, df_mob_suite)
        
        if df_combined.empty:
            logger.error("No se pudieron procesar datos. Saliendo.")
            sys.exit(1)
        
        # Guardar archivo TSV completo
        tsv_output = os.path.join(args.output_dir, "plasmid_analysis_summary.tsv")
        df_combined.to_csv(tsv_output, sep='\t', index=False, na_rep="N/A")
        logger.info(f"Archivo TSV completo guardado en: {tsv_output}")
        
        # Crear versión simplificada para Excel
        logger.info("=== Creando versión simplificada para Excel ===")
        df_excel = simplify_and_format_for_excel(df_combined)
        
        # Seleccionar columnas principales para Excel
        excel_columns = [
            "Sample_ID",
            "Num_Filtered_Plasmid_Contigs",
            "Total_Filtered_Plasmid_Length",
            "Circular_Filtered_Plasmids_Count",
            "Unique_Inc_Types_on_Plasmids",
            "Unique_Replication_Hits_on_Plasmids",
            "Unique_Mobilization_Hits_on_Plasmids",
            "MOB_Suite_Num_Predicted_Plasmids_in_Sample",
            "MOB_Suite_Total_Size_Predicted_Plasmids_in_Sample",
            "Unique_predicted_mobility"
        ]
        
        # Filtrar columnas que existen
        available_columns = [col for col in excel_columns if col in df_excel.columns]
        df_excel_filtered = df_excel[available_columns]
        
        # Guardar archivo Excel
        excel_output = os.path.join(args.output_dir, "Reporte_Plasmidos_Simplificado.xlsx")
        
        try:
            with pd.ExcelWriter(excel_output, engine='openpyxl') as writer:
                df_excel_filtered.to_excel(writer, sheet_name='Resumen_Plasmidos', index=False)
                
                # Configurar formato de Excel
                worksheet = writer.sheets['Resumen_Plasmidos']
                
                # Ajustar ancho de columnas
                for column_cells in worksheet.columns:
                    length = max(len(str(cell.value)) for cell in column_cells)
                    worksheet.column_dimensions[column_cells[0].column_letter].width = min(length + 2, 50)
                
                # Configurar wrapping de texto
                try:
                    from openpyxl.styles import Alignment
                    for row in worksheet.iter_rows():
                        for cell in row:
                            cell.alignment = Alignment(wrap_text=True, vertical='top')
                except ImportError:
                    logger.warning("No se pudo importar openpyxl.styles, continuando sin formato avanzado")
            
            logger.info(f"Archivo Excel simplificado guardado en: {excel_output}")
        except Exception as e:
            logger.warning(f"Error creando archivo Excel: {e}. Continuando sin Excel.")
        
        # Estadísticas finales
        logger.info("=== Estadísticas finales ===")
        logger.info(f"Total de muestras procesadas: {len(df_combined)}")
        if 'Num_Filtered_Plasmid_Contigs' in df_combined.columns:
            total_platon_plasmids = df_combined['Num_Filtered_Plasmid_Contigs'].sum()
            logger.info(f"Total plásmidos detectados por Platon: {total_platon_plasmids}")
        if 'MOB_Suite_Num_Predicted_Plasmids_in_Sample' in df_combined.columns:
            total_mob_plasmids = df_combined['MOB_Suite_Num_Predicted_Plasmids_in_Sample'].sum()
            logger.info(f"Total plásmidos predichos por MOB-suite: {total_mob_plasmids}")
        
        logger.info("=== Análisis completado exitosamente ===")
        
    except Exception as e:
        logger.error(f"Error durante el procesamiento: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
