import pandas as pd
import json
import os
import argparse

# 解析命令行参数
parser = argparse.ArgumentParser(description='根据 cluster 分配合并 JSON 数据。')
parser.add_argument('-c', '--cluster', type=str, required=True, help='包含 cluster 分配的输入 CSV 文件')
parser.add_argument('-j', '--json', type=str, required=True, help='原始 JSON 数据的输入文件')
parser.add_argument('-o', '--output', type=str, required=True, help='合并数据的输出 JSON 文件')
parser.add_argument('-x', '--suffix', type=str, default='.500.fasta', help='要从序列路径中删除的后缀（默认：.500.fasta）')
parser.add_argument('-t', '--txt', type=str, default='clustered_sequences.txt', help='列出被聚类序列的输出 TXT 文件')
args = parser.parse_args()

# 读取 cluster 分配
cluster_df = pd.read_csv(args.cluster)

# 找出包含多个序列的 clusters
cluster_counts = cluster_df['cluster'].value_counts()
clusters_to_merge = cluster_counts[cluster_counts > 1].index.tolist()

# 创建一个从原始 cluster ID 到新的从 1 开始的 cluster ID 的映射
cluster_mapping = {original_id: new_id+1 for new_id, original_id in enumerate(clusters_to_merge)}

# 过滤 cluster_df 以仅包含需要合并的 clusters
filtered_cluster_df = cluster_df[cluster_df['cluster'].isin(clusters_to_merge)]

# 读取原始 JSON 数据
with open(args.json, 'r') as json_file:
    round_data = json.load(json_file)

# 创建一个字典来保存每个 cluster 的合并元素
cluster_dict = {}

# 记录已合并的序列
merged_sequences = {}

# 记录已使用的 round_id
used_round_ids = set()

# 遍历过滤后的 cluster 分配并合并 JSON 元素
for _, row in filtered_cluster_df.iterrows():
    sequence_path = row['sequence']
    original_cluster_id = row['cluster']
    new_cluster_id = cluster_mapping[original_cluster_id]
    cluster_id = f"cluster_{new_cluster_id}"

    # 使用指定的后缀从序列路径中提取 round ID
    round_id = os.path.basename(sequence_path).replace(args.suffix, '')

    # 将此 round 的元素添加到对应的 cluster
    if cluster_id not in cluster_dict:
        cluster_dict[cluster_id] = []

    if round_id in round_data:
        cluster_dict[cluster_id].extend(round_data[round_id])
        used_round_ids.add(round_id)

        # 记录哪些序列被合并
        if cluster_id not in merged_sequences:
            merged_sequences[cluster_id] = []
        merged_sequences[cluster_id].append(round_id)

# 从原始数据中移除已使用的 round_id
for round_id in used_round_ids:
    del round_data[round_id]

# 将合并的数据保存到新的 JSON 文件中（仅包含合并的 clusters）
with open(args.output, 'w') as output_file:
    json.dump(cluster_dict, output_file, indent=4)

# 将被聚类的序列列表保存到 TXT 文件中
with open(args.txt, 'w') as txt_file:
    for cluster_id, sequences in merged_sequences.items():
        txt_file.write(f"{cluster_id}: {', '.join(sequences)}\n")

print(f"合并完成。结果已保存到 '{args.output}'。")
print(f"被聚类的序列列表已保存到 '{args.txt}'。")

