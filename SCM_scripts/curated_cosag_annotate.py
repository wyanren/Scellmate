#!/usr/bin/env python3

import csv
import argparse
from collections import defaultdict

def parse_taxonomy(taxonomy_str):
    taxonomy = {}
    for taxon in taxonomy_str.strip().split(';'):
        if '__' in taxon:
            prefix, name = taxon.split('__', 1)
            taxonomy[prefix + '__'] = name
    return taxonomy

def read_checkm_results(checkm_file):
    checkm_data = {}
    with open(checkm_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        bin_id_idx = headers.index('Bin Id')
        completeness_idx = headers.index('Completeness')
        contamination_idx = headers.index('Contamination')
        for row in reader:
            if not row:
                continue
            bin_id = row[bin_id_idx]
            completeness = float(row[completeness_idx])
            contamination = float(row[contamination_idx])
            checkm_data[bin_id] = {
                'completeness': completeness,
                'contamination': contamination
            }
    return checkm_data

def main(input_file, output_file, checkm_file, max_CoSAG_cont):
    ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    higher_ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__']

    rank_names = {
        'd__': 'Domain',
        'p__': 'Phylum',
        'c__': 'Class',
        'o__': 'Order',
        'f__': 'Family',
        'g__': 'Genus',
        's__': 'Species'
    }

    checkm_data = read_checkm_results(checkm_file)

    taxon_counts = {rank: defaultdict(int) for rank in ranks}

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader)
        for row in reader:
            taxonomy_str = row[1]
            taxonomy = parse_taxonomy(taxonomy_str)
            for rank in ranks:
                taxon = taxonomy.get(rank, '')
                if taxon:
                    taxon_counts[rank][taxon] += 1

    taxon_sp_counts = defaultdict(lambda: defaultdict(int))

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        writer.writerow(['bin_name', 'taxonomy', 'annotate', 'level'])
        next(reader)

        for row in reader:
            bin_name = row[0]
            taxonomy_str = row[1]
            taxonomy = parse_taxonomy(taxonomy_str)

            completeness = checkm_data.get(bin_name, {}).get('completeness', 0)
            contamination = checkm_data.get(bin_name, {}).get('contamination', 100)

            annotate = ''
            level = ''

            if 's__' in taxonomy and taxonomy['s__'] and contamination < max_CoSAG_cont:
                annotate = taxonomy['s__']
                level = 'Species'

            elif completeness > 50 and contamination < max_CoSAG_cont:
                for rank in ['g__', 'f__', 'o__', 'c__', 'p__']:
                    taxon = taxonomy.get(rank, '')
                    if taxon:
                        taxon_sp_counts[rank][taxon] += 1
                        sp_number = taxon_sp_counts[rank][taxon]
                        annotate = f"{taxon} sp. {sp_number}"
                        level = 'Species'
                        break
                else:
                    taxon_sp_counts['Unclassified']['Unclassified'] += 1
                    sp_number = taxon_sp_counts['Unclassified']['Unclassified']
                    annotate = f"Unclassified sp. {sp_number}"
                    level = 'Unclassified'

            else:
                for rank in reversed(higher_ranks):
                    taxon = taxonomy.get(rank, '')
                    if taxon:
                        count = taxon_counts[rank][taxon]
                        rank_name = rank_names[rank]
                        if count == 1:
                            annotate = taxon
                            level = rank_name
                        else:
                            annotate = f"{taxon} spp."
                            level = rank_name
                        break
                else:
                    annotate = 'Unclassified'
                    level = 'Unclassified'

            writer.writerow([bin_name, taxonomy_str, annotate, level])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', default='gtdbtk_2.summary.tsv')
    parser.add_argument('-o', '--output_file', default='annotated_bins.tsv')
    parser.add_argument('-q', '--checkm_file', default='qa_2.txt')
    parser.add_argument('--max_CoSAG_cont', type=float, default=10)
    args = parser.parse_args()

    main(args.input_file, args.output_file, args.checkm_file, args.max_CoSAG_cont)
