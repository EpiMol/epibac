#!/usr/bin/env python3
"""
Script para resumir análisis de elementos móviles genéticos (MGE)
Consolida resultados de plásmidos, transposones, integrones, etc.
"""

import pandas as pd
import os
import sys
import logging
from pathlib import Path

def setup_logging(log_file):
    """Configurar logging para el script."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def parse_platon_results(platon_file):
    """
    Parsea resultados de Platon.
    
    Args:
        platon_file (str): Ruta al archivo TSV de Platon
        
    Returns:
        dict: Diccionario con resumen de plásmidos
    """
    platon_data = {
        'plasmids_detected': 0,
        'potential_plasmids': [],
        'plasmid_summary': [],
        'total_plasmid_size': 0,
        'avg_plasmid_size': 0,
        'amr_bearing_plasmids': 0
    }
    
    try:
        if os.path.exists(platon_file) and os.path.getsize(platon_file) > 0:
            df = pd.read_csv(platon_file, sep='\t')
            
            if not df.empty:
                # Filtrar contigs que tengan evidencia de ser plásmidos
                # Criterios: tienen genes de replicación, movilización, conjugación, o hits de plásmidos
                plasmid_criteria = (
                    (df.get('# Replication', 0) > 0) |
                    (df.get('# Mobilization', 0) > 0) |
                    (df.get('# Conjugation', 0) > 0) |
                    (df.get('# Plasmid Hits', 0) > 0) |
                    (df.get('Inc Type(s)', '0') != '0')
                )
                
                potential_plasmids = df[plasmid_criteria] if plasmid_criteria.any() else pd.DataFrame()
                
                if not potential_plasmids.empty:
                    platon_data['plasmids_detected'] = len(potential_plasmids)
                    platon_data['total_plasmid_size'] = potential_plasmids['Length'].sum()
                    platon_data['avg_plasmid_size'] = round(potential_plasmids['Length'].mean(), 0)
                    platon_data['amr_bearing_plasmids'] = len(potential_plasmids[potential_plasmids.get('# AMRs', 0) > 0])
                    
                    # Crear resumen formateado de plásmidos
                    for _, row in potential_plasmids.iterrows():
                        features = []
                        if row.get('# Replication', 0) > 0:
                            features.append(f"Rep({row['# Replication']})")
                        if row.get('# Mobilization', 0) > 0:
                            features.append(f"Mob({row['# Mobilization']})")
                        if row.get('# Conjugation', 0) > 0:
                            features.append(f"Conj({row['# Conjugation']})")
                        if row.get('# AMRs', 0) > 0:
                            features.append(f"AMR({row['# AMRs']})")
                        if row.get('Inc Type(s)', '0') != '0':
                            features.append(f"Inc({row['Inc Type(s)']})")
                        
                        size_kb = round(row['Length'] / 1000, 1)
                        coverage = row.get('Coverage', 'N/A')
                        feature_str = ','.join(features) if features else 'unknown'
                        
                        plasmid_summary = f"Contig_{row['ID']} [{size_kb}kb,{coverage}x,{feature_str}]"
                        platon_data['plasmid_summary'].append(plasmid_summary)
                        
    except Exception as e:
        logging.warning(f"Error parsing Platon results: {e}")
    
    return platon_data

def parse_mob_suite_results(mob_suite_dir):
    """
    Parsea resultados completos de MOB-suite incluyendo múltiples archivos.
    
    Args:
        mob_suite_dir (str): Directorio con todos los archivos de MOB-suite
        
    Returns:
        dict: Diccionario con resumen completo de tipificación de plásmidos
    """
    mob_data = {
        'mob_plasmids_detected': 0,
        'conjugative_plasmids': 0,
        'mobilizable_plasmids': 0,
        'plasmid_summary': [],
        'inc_types': set(),
        'mob_types': set(),
        'host_ranges': set(),
        'mge_elements': [],
        'total_plasmid_size': 0,
        'avg_plasmid_size': 0,
        'mash_neighbors': set(),
        'cluster_diversity': 0
    }
    
    try:
        # Parsear archivo principal mobtyper_results.txt
        mobtyper_file = os.path.join(mob_suite_dir, 'mobtyper_results.txt')
        if os.path.exists(mobtyper_file) and os.path.getsize(mobtyper_file) > 0:
            df = pd.read_csv(mobtyper_file, sep='\t')
            
            if not df.empty:
                mob_data['mob_plasmids_detected'] = len(df)
                
                # Calcular estadísticas de tamaño
                sizes = df.get('size', []).dropna()
                if len(sizes) > 0:
                    mob_data['total_plasmid_size'] = int(sizes.sum())
                    mob_data['avg_plasmid_size'] = round(sizes.mean(), 0)
                
                # Contar tipos de movilidad
                mob_data['conjugative_plasmids'] = len(df[df.get('predicted_mobility', '') == 'conjugative'])
                mob_data['mobilizable_plasmids'] = len(df[df.get('predicted_mobility', '') == 'mobilizable'])
                
                # Extraer diversidad de clusters
                primary_clusters = df.get('primary_cluster_id', []).dropna().unique()
                mob_data['cluster_diversity'] = len(primary_clusters)
                
                # Crear resumen formateado de plásmidos MOB-suite
                for _, row in df.iterrows():
                    plasmid_id = str(row.get('sample_id', 'unknown')).replace('plasmid_', '')
                    size_kb = round(row.get('size', 0) / 1000, 1)
                    mobility = str(row.get('predicted_mobility', 'unknown')).strip()
                    gc_content = round(row.get('gc', 0) * 100, 1)
                    
                    # Extraer tipos Inc/Rep
                    rep_types = str(row.get('rep_type(s)', '')).split(',') if pd.notna(row.get('rep_type(s)')) else []
                    inc_types = [rep.strip() for rep in rep_types if rep.strip() and rep.strip() != '-']
                    
                    # Extraer información de relaxasa
                    relaxase_types = str(row.get('relaxase_type(s)', '')).split(',') if pd.notna(row.get('relaxase_type(s)')) else []
                    relaxases = [rel.strip() for rel in relaxase_types if rel.strip() and rel.strip() != '-']
                    
                    # Extraer rango de hospedador
                    host_range = str(row.get('predicted_host_range_overall_name', '')).strip()
                    if host_range and host_range != '-':
                        mob_data['host_ranges'].add(host_range)
                    
                    # Extraer vecino más cercano por MASH
                    mash_neighbor = str(row.get('mash_neighbor_identification', '')).strip()
                    if mash_neighbor and mash_neighbor != '-':
                        mob_data['mash_neighbors'].add(mash_neighbor)
                    
                    # Crear resumen del plásmido con información rica
                    features = []
                    if inc_types:
                        features.extend(inc_types[:2])  # Limitar a 2 tipos Inc para brevedad
                        mob_data['inc_types'].update(inc_types)
                    if relaxases:
                        features.append(f"Rel:{relaxases[0]}")  # Usar primer tipo de relaxasa
                        mob_data['mob_types'].update(relaxases)
                    
                    # Añadir información de movilidad y GC
                    if mobility and mobility != 'unknown':
                        features.append(mobility.capitalize())
                    features.append(f"GC:{gc_content}%")
                    
                    feature_str = ','.join(features) if features else 'untyped'
                    plasmid_summary = f"{plasmid_id} [{size_kb}kb,{feature_str}]"
                    mob_data['plasmid_summary'].append(plasmid_summary)
        
        # Parsear archivo mge.report.txt para elementos móviles adicionales
        mge_file = os.path.join(mob_suite_dir, 'mge.report.txt')
        if os.path.exists(mge_file) and os.path.getsize(mge_file) > 0:
            try:
                mge_df = pd.read_csv(mge_file, sep='\t')
                if not mge_df.empty:
                    # Extraer información de elementos móviles
                    for _, mge_row in mge_df.iterrows():
                        mge_type = str(mge_row.get('mge_type', 'unknown')).strip()
                        mge_subtype = str(mge_row.get('mge_subtype', '')).strip()
                        mge_length = mge_row.get('mge_length', 0)
                        molecule_type = str(mge_row.get('molecule_type', 'unknown')).strip()
                        pident = mge_row.get('pident', 0)
                        
                        if mge_type and mge_type != 'unknown':
                            mge_element = f"{mge_type}({mge_subtype}) [{mge_length}bp,{pident:.1f}%,{molecule_type}]"
                            mob_data['mge_elements'].append(mge_element)
            except Exception as e:
                logging.warning(f"Error parsing MGE report: {e}")
                
    except Exception as e:
        logging.warning(f"Error parsing MOB-suite results: {e}")
    
    # Convertir sets a listas para serialización
    mob_data['inc_types'] = list(mob_data['inc_types'])
    mob_data['mob_types'] = list(mob_data['mob_types'])
    mob_data['host_ranges'] = list(mob_data['host_ranges'])
    mob_data['mash_neighbors'] = list(mob_data['mash_neighbors'])
    
    return mob_data

def parse_transposon_results(transposon_file):
    """
    Parsea resultados de análisis de transposones.
    
    Args:
        transposon_file (str): Ruta al archivo TSV de abricate para transposones
        
    Returns:
        dict: Diccionario con resumen de transposones
    """
    te_data = {
        'transposons_detected': 0,
        'transposon_summary': [],
        'te_families': set(),
        'avg_coverage': 0,
        'avg_identity': 0,
        'high_quality_hits': 0  # >90% coverage y >95% identity
    }
    
    try:
        if os.path.exists(transposon_file) and os.path.getsize(transposon_file) > 0:
            df = pd.read_csv(transposon_file, sep='\t')
            
            if not df.empty:
                te_data['transposons_detected'] = len(df)
                
                # Calcular estadísticas
                coverages = df.get('%COVERAGE', [])
                identities = df.get('%IDENTITY', [])
                
                if len(coverages) > 0:
                    te_data['avg_coverage'] = round(sum(coverages) / len(coverages), 1)
                if len(identities) > 0:
                    te_data['avg_identity'] = round(sum(identities) / len(identities), 1)
                
                # Contar hits de alta calidad
                high_quality = df[
                    (df.get('%COVERAGE', 0) >= 90) & 
                    (df.get('%IDENTITY', 0) >= 95)
                ]
                te_data['high_quality_hits'] = len(high_quality)
                
                # Crear resumen formateado: "Nombre [tamaño_pb, %identidad]"
                for _, row in df.iterrows():
                    gene_name = str(row.get('GENE', 'unknown')).strip()
                    coverage = row.get('%COVERAGE', 0)
                    identity = row.get('%IDENTITY', 0)
                    start = row.get('START', 0)
                    end = row.get('END', 0)
                    
                    # Calcular tamaño del hit
                    hit_size = abs(end - start) if start and end else 0
                    
                    # Formatear según especificación: "Tn1 [1231pb,98%]"
                    transposon_summary = f"{gene_name} [{hit_size}pb,{identity:.1f}%]"
                    te_data['transposon_summary'].append(transposon_summary)
                    
                    # Añadir familia de transposón
                    if gene_name and gene_name != 'unknown':
                        te_data['te_families'].add(gene_name)
                        
    except Exception as e:
        logging.warning(f"Error parsing transposon results: {e}")
    
    # Convertir set a lista para serialización
    te_data['te_families'] = list(te_data['te_families'])
    
    return te_data

def create_summary_table(sample_id, platon_data, mob_data, te_data):
    """
    Crea tabla resumen consolidada de MGE con formato optimizado.
    
    Args:
        sample_id (str): ID de la muestra
        platon_data (dict): Datos de Platon
        mob_data (dict): Datos de MOB-suite
        te_data (dict): Datos de transposones
        
    Returns:
        pd.DataFrame: DataFrame con resumen consolidado
    """
    summary_data = []
    
    # Resumen de plásmidos Platon
    platon_details = '; '.join(platon_data['plasmid_summary']) if platon_data['plasmid_summary'] else 'None'
    summary_data.append({
        'sample_id': sample_id,
        'mge_type': 'PLASMIDS',
        'tool': 'Platon',
        'elements_detected': platon_data['plasmids_detected'],
        'details': platon_details,
        'additional_info': f"Total size: {platon_data.get('total_plasmid_size', 0)}bp; AMR-bearing: {platon_data.get('amr_bearing_plasmids', 0)}"
    })
    
    # Resumen mejorado de plásmidos MOB-suite
    mob_details = '; '.join(mob_data['plasmid_summary']) if mob_data['plasmid_summary'] else 'None'
    mob_info_parts = [
        f"Conjugative: {mob_data.get('conjugative_plasmids', 0)}",
        f"Mobilizable: {mob_data.get('mobilizable_plasmids', 0)}",
        f"Total size: {mob_data.get('total_plasmid_size', 0)}bp"
    ]
    
    if mob_data.get('inc_types'):
        mob_info_parts.append(f"Inc types: {', '.join(mob_data['inc_types'][:3])}")
    if mob_data.get('mash_neighbors'):
        mob_info_parts.append(f"MASH neighbors: {', '.join(list(mob_data['mash_neighbors'])[:2])}")
    if mob_data.get('cluster_diversity', 0) > 0:
        mob_info_parts.append(f"Cluster diversity: {mob_data['cluster_diversity']}")
    
    mob_info = '; '.join(mob_info_parts)
    
    summary_data.append({
        'sample_id': sample_id,
        'mge_type': 'PLASMIDS',
        'tool': 'MOB-suite',
        'elements_detected': mob_data['mob_plasmids_detected'],
        'details': mob_details,
        'additional_info': mob_info
    })
    
    # Resumen de elementos móviles adicionales de MOB-suite
    if mob_data.get('mge_elements'):
        mge_details = '; '.join(mob_data['mge_elements'][:5])  # Limitar a 5 elementos
        summary_data.append({
            'sample_id': sample_id,
            'mge_type': 'INSERTION_SEQUENCES',
            'tool': 'MOB-suite',
            'elements_detected': len(mob_data['mge_elements']),
            'details': mge_details,
            'additional_info': f"Total MGE elements detected in plasmids and chromosome"
        })
    
    # Resumen de transposones con formato optimizado
    te_details = '; '.join(te_data['transposon_summary']) if te_data['transposon_summary'] else 'None'
    te_info = f"Families: {len(te_data.get('te_families', []))}; Avg identity: {te_data.get('avg_identity', 0)}%; High quality: {te_data.get('high_quality_hits', 0)}"
    
    summary_data.append({
        'sample_id': sample_id,
        'mge_type': 'TRANSPOSONS',
        'tool': 'Abricate',
        'elements_detected': te_data['transposons_detected'],
        'details': te_details,
        'additional_info': te_info
    })
    
    return pd.DataFrame(summary_data)

def create_detailed_summary(sample_id, platon_data, mob_data, te_data):
    """
    Crea resumen detallado con estadísticas adicionales optimizadas.
    
    Args:
        sample_id (str): ID de la muestra
        platon_data (dict): Datos de Platon
        mob_data (dict): Datos de MOB-suite  
        te_data (dict): Datos de transposones
        
    Returns:
        pd.DataFrame: DataFrame con resumen detallado
    """
    detailed_data = []
    
    # Estadísticas de plásmidos Platon
    if platon_data['plasmids_detected'] > 0:
        detailed_data.extend([
            {
                'sample_id': sample_id,
                'category': 'plasmid_statistics',
                'metric': 'total_plasmids',
                'value': platon_data['plasmids_detected'],
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'plasmid_statistics',
                'metric': 'total_plasmid_size',
                'value': platon_data.get('total_plasmid_size', 0),
                'unit': 'bp'
            },
            {
                'sample_id': sample_id,
                'category': 'plasmid_statistics',
                'metric': 'average_plasmid_size',
                'value': platon_data.get('avg_plasmid_size', 0),
                'unit': 'bp'
            },
            {
                'sample_id': sample_id,
                'category': 'plasmid_statistics',
                'metric': 'amr_bearing_plasmids',
                'value': platon_data.get('amr_bearing_plasmids', 0),
                'unit': 'count'
            }
        ])
    
    # Estadísticas mejoradas de MOB-suite
    if mob_data['mob_plasmids_detected'] > 0:
        detailed_data.extend([
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'typed_plasmids',
                'value': mob_data['mob_plasmids_detected'],
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'conjugative_plasmids',
                'value': mob_data.get('conjugative_plasmids', 0),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'mobilizable_plasmids',
                'value': mob_data.get('mobilizable_plasmids', 0),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'total_mob_plasmid_size',
                'value': mob_data.get('total_plasmid_size', 0),
                'unit': 'bp'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'avg_mob_plasmid_size',
                'value': mob_data.get('avg_plasmid_size', 0),
                'unit': 'bp'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'inc_types_detected',
                'value': len(mob_data.get('inc_types', [])),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'relaxase_types_detected',
                'value': len(mob_data.get('mob_types', [])),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'host_ranges_predicted',
                'value': len(mob_data.get('host_ranges', [])),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'cluster_diversity',
                'value': mob_data.get('cluster_diversity', 0),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'mobilization_statistics',
                'metric': 'mash_neighbors_identified',
                'value': len(mob_data.get('mash_neighbors', [])),
                'unit': 'count'
            }
        ])
    
    # Estadísticas de elementos móviles adicionales de MOB-suite
    if mob_data.get('mge_elements'):
        detailed_data.append({
            'sample_id': sample_id,
            'category': 'insertion_sequence_statistics',
            'metric': 'total_mge_elements',
            'value': len(mob_data['mge_elements']),
            'unit': 'count'
        })
    
    # Estadísticas de transposones optimizadas
    if te_data['transposons_detected'] > 0:
        detailed_data.extend([
            {
                'sample_id': sample_id,
                'category': 'transposon_statistics',
                'metric': 'total_transposons',
                'value': te_data['transposons_detected'],
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'transposon_statistics',
                'metric': 'transposon_families',
                'value': len(te_data.get('te_families', [])),
                'unit': 'count'
            },
            {
                'sample_id': sample_id,
                'category': 'transposon_statistics',
                'metric': 'average_identity',
                'value': te_data.get('avg_identity', 0),
                'unit': 'percent'
            },
            {
                'sample_id': sample_id,
                'category': 'transposon_statistics',
                'metric': 'average_coverage',
                'value': te_data.get('avg_coverage', 0),
                'unit': 'percent'
            },
            {
                'sample_id': sample_id,
                'category': 'transposon_statistics',
                'metric': 'high_quality_hits',
                'value': te_data.get('high_quality_hits', 0),
                'unit': 'count'
            }
        ])
    
    return pd.DataFrame(detailed_data)

def main():
    """Función principal del script."""
    # Obtener argumentos de Snakemake
    try:
        sample_id = snakemake.wildcards.sample
        log_file = snakemake.log[0]
        
        # Archivos de entrada
        platon_file = snakemake.input.platon
        mob_suite_dir = os.path.dirname(snakemake.input.mob_suite)  # Usar directorio en lugar del archivo
        transposon_file = snakemake.input.transposons
        
        # Archivo de salida
        output_file = snakemake.output.summary
        
    except NameError:
        # Para testing manual sin Snakemake
        import sys
        if len(sys.argv) >= 6:
            sample_id = sys.argv[1]
            platon_file = sys.argv[2]
            mob_suite_dir = sys.argv[3]  # Ahora espera directorio
            transposon_file = sys.argv[4]
            output_file = sys.argv[5]
            log_file = sys.argv[6] if len(sys.argv) > 6 else "mge_analysis.log"
        else:
            print("Usage: python script.py sample_id platon_file mob_suite_dir transposon_file output_file [log_file]")
            sys.exit(1)
    
    # Configurar logging
    logger = setup_logging(log_file)
    logger.info(f"Iniciando resumen de análisis MGE para muestra: {sample_id}")
    
    logger.info(f"Procesando archivos de entrada:")
    logger.info(f"  - Platon: {platon_file}")
    logger.info(f"  - MOB-suite directory: {mob_suite_dir}")
    logger.info(f"  - Transposones: {transposon_file}")
    
    try:
        # Parsear resultados
        logger.info("Parseando resultados de Platon...")
        platon_data = parse_platon_results(platon_file)
        
        logger.info("Parseando resultados de MOB-suite...")
        mob_data = parse_mob_suite_results(mob_suite_dir)  # Ahora usa directorio
        
        logger.info("Parseando resultados de transposones...")
        te_data = parse_transposon_results(transposon_file)
        
        # Crear tabla resumen
        logger.info("Creando tabla resumen...")
        summary_df = create_summary_table(sample_id, platon_data, mob_data, te_data)
        
        # Crear directorio de salida si no existe
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Guardar resumen
        summary_df.to_csv(output_file, sep='\t', index=False)
        
        # Log de estadísticas mejoradas
        logger.info(f"Resumen completado:")
        logger.info(f"  - Plásmidos detectados (Platon): {platon_data['plasmids_detected']}")
        logger.info(f"  - Plásmidos tipificados (MOB-suite): {mob_data['mob_plasmids_detected']}")
        logger.info(f"    * Conjugativos: {mob_data.get('conjugative_plasmids', 0)}")
        logger.info(f"    * Movilizables: {mob_data.get('mobilizable_plasmids', 0)}")
        logger.info(f"    * Tipos Inc detectados: {len(mob_data.get('inc_types', []))}")
        logger.info(f"    * Elementos móviles adicionales: {len(mob_data.get('mge_elements', []))}")
        logger.info(f"  - Transposones detectados: {te_data['transposons_detected']}")
        logger.info(f"    * Familias únicas: {len(te_data.get('te_families', []))}")
        logger.info(f"    * Hits de alta calidad: {te_data.get('high_quality_hits', 0)}")
        logger.info(f"  - Archivo de salida: {output_file}")
        
        # Crear también archivo detallado (opcional)
        detailed_file = output_file.replace('.tsv', '_detailed.tsv')
        detailed_df = create_detailed_summary(sample_id, platon_data, mob_data, te_data)
        detailed_df.to_csv(detailed_file, sep='\t', index=False)
        logger.info(f"  - Archivo detallado: {detailed_file}")
        
    except Exception as e:
        logger.error(f"Error durante el procesamiento: {e}")
        # Crear archivo vacío para que Snakemake no falle
        with open(output_file, 'w') as f:
            f.write("sample_id\tmge_type\ttool\telements_detected\tdetails\tadditional_info\n")
            f.write(f"{sample_id}\tERROR\tNA\t0\tProcessing failed: {e}\t\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
