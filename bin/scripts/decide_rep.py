#!/usr/bin/env python3

import sys
import re
import argparse

def main(species_file, qa_file, output_file):
    # 读取 qa.txt，建立 Bin Id 到数据的映射
    qa_data = {}
    with open(qa_file, 'r') as f:
        header_line = f.readline().strip()
        # 使用制表符或多个空格作为分隔符
        header = re.split(r'\t+|\s{2,}', header_line)
        for line in f:
            line = line.strip()
            if not line:
                continue
            # 使用正则表达式分割行
            fields = re.split(r'\t+|\s{2,}', line)
            if len(fields) < 7:
                print(f'警告：无法解析行：{line}', file=sys.stderr)
                continue  # 跳过不完整的行
            bin_id = fields[0]
            try:
                completeness = float(fields[5])
                contamination = float(fields[6])
            except ValueError:
                print(f'警告：无法解析 {bin_id} 的 Completeness 或 Contamination', file=sys.stderr)
                continue
            qa_data[bin_id] = {
                'line': line,
                'Completeness': completeness,
                'Contamination': contamination,
                'fields': fields,
            }

    # 读取 species_to_nocluster_deduplicated.tsv，建立物种到聚类ID的映射
    species_clusters = {}
    with open(species_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                print(f'警告：物种文件中的无效行：{line}', file=sys.stderr)
                continue
            cluster_id, species_name = parts
            species_clusters.setdefault(species_name, []).append(cluster_id)

    # 打开输出文件
    with open(output_file, 'w') as out_f:
        # 写入标题行，添加 Species 和 Grade 列
        out_f.write('\t'.join(header + ['Species', 'Grade']) + '\n')

        # 对于每个物种，选择代表性基因组
        for species_name, cluster_ids in species_clusters.items():
            genomes = []
            for cluster_id in cluster_ids:
                if cluster_id in qa_data:
                    completeness = qa_data[cluster_id]['Completeness']
                    contamination = qa_data[cluster_id]['Contamination']
                    # 评定等级
                    if completeness > 90 and contamination < 5:
                        grade = 'High'
                    elif completeness < 50 or contamination > 10:
                        grade = 'Low'
                    else:
                        grade = 'Medium'
                    # 计算分数
                    score = completeness * (1 - contamination * 5)
                    genomes.append({
                        'cluster_id': cluster_id,
                        'grade': grade,
                        'score': score,
                        'qa_line': qa_data[cluster_id]['line'],
                        'fields': qa_data[cluster_id]['fields'],
                    })
                else:
                    print(f'警告：在 qa.txt 中找不到 Cluster ID {cluster_id}', file=sys.stderr)

            if not genomes:
                print(f'警告：物种 {species_name} 没有可用的基因组数据', file=sys.stderr)
                continue

            # 确定最高可用等级
            grades = [g['grade'] for g in genomes]
            if 'High' in grades:
                best_grade = 'High'
            elif 'Medium' in grades:
                best_grade = 'Medium'
            else:
                best_grade = 'Low'

            # 筛选出具有最高等级的基因组
            best_genomes = [g for g in genomes if g['grade'] == best_grade]

            # 选择分数最高的基因组作为代表
            best_genomes.sort(key=lambda x: x['score'], reverse=True)
            representative = best_genomes[0]

            # 将代表性基因组的信息写入输出文件，添加物种和等级信息
            qa_fields = representative['fields']
            output_line = '\t'.join(qa_fields + [species_name, representative['grade']])
            out_f.write(output_line + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='根据 Completeness 和 Contamination 选择代表性基因组。')
    parser.add_argument('-s', '--species_file', default='species_to_nocluster_deduplicated.tsv', help='输入的物种到聚类文件')
    parser.add_argument('-q', '--qa_file', default='qa.txt', help='输入的 qa.txt 文件')
    parser.add_argument('-o', '--output_file', default='qa-representative.txt', help='输出的代表性基因组文件')
    args = parser.parse_args()

    main(args.species_file, args.qa_file, args.output_file)

