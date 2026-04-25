#!/usr/bin/env python3
"""
Mitophen Database Parser
Extracts variant-phenotype associations from mitophen SQL dump

Creates a JSON file mapping mitochondrial variants to disease phenotypes (HPO terms)
"""

import gzip
import re
import json
import argparse
from pathlib import Path
from collections import defaultdict


def parse_sql_insert(line, table_name):
    """Parse SQL INSERT statement and extract values"""
    # Pattern: INSERT INTO `table` VALUES (val1,val2,...),(val3,val4,...);
    pattern = r"INSERT INTO `" + table_name + r"` VALUES\s*(.+);"
    match = re.search(pattern, line)
    if not match:
        return []

    values_str = match.group(1)
    rows = []

    # Parse each row (val1,val2,...)
    # Handle nested parentheses and quoted strings
    pos = 0
    while pos < len(values_str):
        if values_str[pos] == '(':
            # Find matching closing parenthesis
            depth = 1
            start = pos + 1
            pos += 1
            while pos < len(values_str) and depth > 0:
                if values_str[pos] == '(':
                    depth += 1
                elif values_str[pos] == ')':
                    depth -= 1
                pos += 1
            row_str = values_str[start:pos-1]
            rows.append(parse_row_values(row_str))
        else:
            pos += 1

    return rows


def parse_row_values(row_str):
    """Parse individual row values, handling quotes and NULL"""
    values = []
    pos = 0
    current_val = ""

    while pos < len(row_str):
        char = row_str[pos]

        if char == "'":
            # Quoted string - find end quote
            pos += 1
            quoted_val = ""
            while pos < len(row_str):
                if row_str[pos] == "'" and (pos + 1 >= len(row_str) or row_str[pos + 1] != "'"):
                    break
                if row_str[pos] == "'" and row_str[pos + 1] == "'":
                    quoted_val += "'"
                    pos += 2
                else:
                    quoted_val += row_str[pos]
                    pos += 1
            values.append(quoted_val)
            pos += 1  # Skip closing quote
            # Skip comma after quoted value if present
            if pos < len(row_str) and row_str[pos] == ',':
                pos += 1
            current_val = ""
        elif char == ',':
            if current_val.strip():
                values.append(convert_value(current_val.strip()))
            current_val = ""
            pos += 1
        elif char == 'N' and pos + 3 < len(row_str) and row_str[pos:pos+4] == 'NULL':
            values.append(None)
            pos += 4
            # Skip comma after NULL if present
            if pos < len(row_str) and row_str[pos] == ',':
                pos += 1
            current_val = ""
        else:
            current_val += char
            pos += 1

    if current_val.strip():
        values.append(convert_value(current_val.strip()))

    return values


def convert_value(val_str):
    """Convert string value to appropriate type"""
    val_str = val_str.strip()
    if val_str == '':
        return None
    try:
        if '.' in val_str:
            return float(val_str)
        return int(val_str)
    except ValueError:
        return val_str


def parse_mitophen_sql(sql_path):
    """Parse mitophen SQL dump and extract all tables"""
    tables = {
        'variant': {},
        'sample': [],
        'patient': {},
        'patient_hpo': [],
        'hpo': {}
    }

    # Open gzipped SQL file
    opener = gzip.open if sql_path.endswith('.gz') else open
    with opener(sql_path, 'rt', encoding='utf-8') as f:
        current_table = None

        for line in f:
            # Detect which table we're processing
            if 'INSERT INTO `variant`' in line:
                rows = parse_sql_insert(line, 'variant')
                for row in rows:
                    if len(row) >= 3:
                        tables['variant'][row[0]] = {
                            'id': row[0],
                            'pos': row[1],
                            'name': row[2]
                        }

            elif 'INSERT INTO `hpo`' in line:
                rows = parse_sql_insert(line, 'hpo')
                for row in rows:
                    if len(row) >= 2:
                        tables['hpo'][row[0]] = row[1]

            elif 'INSERT INTO `patient`' in line:
                rows = parse_sql_insert(line, 'patient')
                for row in rows:
                    if len(row) >= 1:
                        tables['patient'][row[0]] = {
                            'id': row[0],
                            'family_id': row[1] if len(row) > 1 else None,
                            'gender': row[2] if len(row) > 2 else None,
                            'affected': row[3] if len(row) > 3 else None,
                            'age_at_onset': row[4] if len(row) > 4 else None
                        }

            elif 'INSERT INTO `sample`' in line:
                rows = parse_sql_insert(line, 'sample')
                for row in rows:
                    if len(row) >= 4:
                        tables['sample'].append({
                            'id': row[0],
                            'patient_id': row[1],
                            'variant_id': row[2],
                            'tissue_id': row[3],
                            'heteroplasmy': row[4] if len(row) > 4 else None,
                            'variant_presence': row[5] if len(row) > 5 else None
                        })

            elif 'INSERT INTO `patient_hpo`' in line:
                rows = parse_sql_insert(line, 'patient_hpo')
                for row in rows:
                    if len(row) >= 2:
                        tables['patient_hpo'].append({
                            'patient_id': row[0],
                            'hpo_code': row[1],
                            'minimal_set': row[2] if len(row) > 2 else None,
                            'lead_term': row[3] if len(row) > 3 else None
                        })

    return tables


def build_variant_phenotype_map(tables):
    """Build mapping from variant to phenotypes"""
    # First build patient -> variant mapping
    patient_to_variants = defaultdict(set)
    for sample in tables['sample']:
        if sample['patient_id'] and sample['variant_id']:
            patient_to_variants[sample['patient_id']].add(sample['variant_id'])

    # Then build variant -> phenotypes mapping
    variant_to_hpo = defaultdict(set)
    for ph in tables['patient_hpo']:
        patient_id = ph['patient_id']
        hpo_code = ph['hpo_code']

        # Get variants for this patient
        for variant_id in patient_to_variants.get(patient_id, set()):
            variant_to_hpo[variant_id].add(hpo_code)

    # Build final mapping with HPO names
    variant_map = {}
    for variant_id, hpo_codes in variant_to_hpo.items():
        variant_info = tables['variant'].get(variant_id, {})
        phenotypes = []
        for code in sorted(hpo_codes):
            name = tables['hpo'].get(code, code)
            phenotypes.append({
                'code': code,
                'name': name
            })

        variant_map[variant_id] = {
            'id': variant_id,
            'pos': variant_info.get('pos', 0),
            'name': variant_info.get('name', variant_id),
            'phenotypes': phenotypes,
            'phenotype_count': len(phenotypes)
        }

    return variant_map


def build_position_variant_map(variant_map):
    """Build mapping from position to variants for quick lookup"""
    pos_map = defaultdict(list)
    for variant_id, info in variant_map.items():
        pos = info.get('pos', 0)
        pos_map[pos].append({
            'variant_id': variant_id,
            'name': info['name'],
            'phenotypes': info['phenotypes']
        })
    return dict(pos_map)


def main():
    parser = argparse.ArgumentParser(description='Parse mitophen SQL dump')
    parser.add_argument('-i', '--input', required=True, help='Input SQL file (gzipped)')
    parser.add_argument('-o', '--output', required=True, help='Output JSON file')

    args = parser.parse_args()

    sql_path = args.input
    output_path = args.output

    print(f"Parsing mitophen SQL: {sql_path}")
    tables = parse_mitophen_sql(sql_path)

    print(f"Found {len(tables['variant'])} variants")
    print(f"Found {len(tables['hpo'])} HPO terms")
    print(f"Found {len(tables['patient'])} patients")
    print(f"Found {len(tables['sample'])} samples")
    print(f"Found {len(tables['patient_hpo'])} patient-HPO associations")

    # Build variant-phenotype mapping
    print("\nBuilding variant-phenotype mapping...")
    variant_map = build_variant_phenotype_map(tables)
    pos_map = build_position_variant_map(variant_map)

    print(f"Built mapping for {len(variant_map)} variants with phenotype data")

    # Save output
    output_data = {
        'variants': variant_map,
        'position_map': pos_map,
        'hpo_dict': tables['hpo'],
        'metadata': {
            'source': 'mitophen',
            'version': '1.7',
            'total_variants': len(tables['variant']),
            'variants_with_phenotypes': len(variant_map),
            'total_hpo_terms': len(tables['hpo'])
        }
    }

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)

    print(f"\nOutput saved to: {output_path}")

    # Print some examples
    print("\n=== Example variants with phenotypes ===")
    for variant_id in list(variant_map.keys())[:5]:
        info = variant_map[variant_id]
        print(f"\n{info['name']} (pos: {info['pos']})")
        print(f"  Phenotypes ({info['phenotype_count']}):")
        for ph in info['phenotypes'][:5]:
            print(f"    - {ph['code']}: {ph['name']}")


if __name__ == '__main__':
    main()