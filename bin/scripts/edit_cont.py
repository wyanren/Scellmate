import json

with open('contamination_map.json', 'r') as f:
    contamination_map = json.load(f)

with open('final.reformat.add_cluster.json', 'r') as f:
    original_data = json.load(f)

final_data = {}
edited_data = {}

for key in original_data:
    if key in contamination_map:
    
        original_values = original_data[key]

        values_to_remove = contamination_map[key]
        
        new_values = [v for v in original_values if v not in values_to_remove]
        
        new_key = key + '_edit'
        
        final_data[new_key] = new_values
        
        edited_data[new_key] = new_values
    else:
        
        final_data[key] = original_data[key]

with open('final.reformat.add_cluster.edit_cont.json', 'w') as f:
    json.dump(final_data, f, indent=2)

with open('edit_cont.json', 'w') as f:
    json.dump(edited_data, f, indent=2)

print("Generate 'final.reformat.add_cluster.edit_cont.json' and 'edit_cont.json'")

