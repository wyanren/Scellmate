#!/usr/bin/env python3

import sys
import re
import argparse

def main(species_file, qa_file, output_file):
    qa_data = {}
    with open(qa_file, 'r') as f:
        header_line = f.readline().strip()
        header = re.split(r'\t+|\s{2,}', header_line)
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = re.split(r'\t+|\s{2,}', line)
            if len(fields) < 7:
                print(f'Warning: failed to parse line: {line}', file=sys.stderr)
                continue
            bin_id = fields[0]
            try:
                completeness = float(fields[5])
                contamination = float(fields[6])
            except ValueError:
                print(f'Warning: failed to parse Completeness or Contamination for {bin_id}', file=sys.stderr)
                continue
            qa_data[bin_id] = {
                'line': line,
                'Completeness': completeness,
                'Contamination': contamination,
                'fields': fields,
            }

    species_clusters = {}
    with open(species_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                print(f'Warning: invalid line in species file: {line}', file=sys.stderr)
                continue
            cluster_id, species_name = parts
            species_clusters.setdefault(species_name, []).append(cluster_id)

    with open(output_file, 'w') as out_f:
        out_f.write('\t'.join(header + ['Species', 'Grade']) + '\n')

        for species_name, cluster_ids in species_clusters.items():
            genomes = []
            for cluster_id in cluster_ids:
                if cluster_id in qa_data:
                    completeness = qa_data[cluster_id]['Completeness']
                    contamination = qa_data[cluster_id]['Contamination']
                    if completeness > 90 and contamination < 5:
                        grade = 'High'
                    elif completeness < 50 or contamination > 10:
                        grade = 'Low'
                    else:
                        grade = 'Medium'
                    score = completeness * (1 - contamination * 5)
                    genomes.append({
                        'cluster_id': cluster_id,
                        'grade': grade,
                        'score': score,
                        'qa_line': qa_data[cluster_id]['line'],
                        'fields': qa_data[cluster_id]['fields'],
                    })
                else:
                    print(f'Warning: Cluster ID {cluster_id} was not found in qa.txt', file=sys.stderr)

            if not genomes:
                print(f'Warning: no available genome data for species {species_name}', file=sys.stderr)
                continue

            grades = [g['grade'] for g in genomes]
            if 'High' in grades:
                best_grade = 'High'
            elif 'Medium' in grades:
                best_grade = 'Medium'
            else:
                best_grade = 'Low'

            best_genomes = [g for g in genomes if g['grade'] == best_grade]

            best_genomes.sort(key=lambda x: x['score'], reverse=True)
            representative = best_genomes[0]

            qa_fields = representative['fields']
            output_line = '\t'.join(qa_fields + [species_name, representative['grade']])
            out_f.write(output_line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Select representative genomes based on Completeness and Contamination.')
    parser.add_argument('-s', '--species_file', default='species_to_nocluster_deduplicated.tsv', help='Input species-to-cluster file')
    parser.add_argument('-q', '--qa_file', default='qa.txt', help='Input qa.txt file')
    parser.add_argument('-o', '--output_file', default='qa-representative.txt', help='Output representative genome file')
    args = parser.parse_args()

    main(args.species_file, args.qa_file, args.output_file)

