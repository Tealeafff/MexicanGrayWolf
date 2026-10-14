import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_heterozygosity(ped_file):
    wolf_stats = {}

    with open(ped_file, 'r') as file:
        for line in file:
            columns = line.strip().split()
            
            if len(columns) < 7:
                continue 

            wolf_id = columns[1]
            genetics = columns[6:]

            heterozygous_count = 0
            total_valid_sites = 0
            
            for i in range(0, len(genetics), 2):
                allele1 = genetics[i]
                allele2 = genetics[i+1]
                
                if allele1 == '0' or allele2 == '0':
                    continue
                    
                total_valid_sites += 1
                
                if allele1 != allele2:
                    heterozygous_count += 1
            
            if total_valid_sites > 0:
                heterozygosity_rate = heterozygous_count / total_valid_sites
                wolf_stats[wolf_id] = heterozygosity_rate
                
    return wolf_stats


results = calculate_heterozygosity("MW.clean.ped")

df = pd.DataFrame(list(results.items()), columns=['Wolf_ID', 'Heterozygosity_Rate'])

df = df.sort_values(by='Heterozygosity_Rate', ascending=False)

print("\nTop 5 Most Genetically Diverse Wolves:")
print(df.head())

def build_breeding_heatmap(target_filename, top_wolves_df):
    print("\nCalculating Pairwise Similarity for the top candidates...")
    
    top_ids = [str(x) for x in top_wolves_df.head(15)['Wolf_ID'].tolist()]
    wolf_genetics = {}
    
    with open(target_filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) < 7: continue
                
            wolf_id = str(columns[1])
            if wolf_id in top_ids:
                wolf_genetics[wolf_id] = columns[6:]
    
    found_ids = [wid for wid in top_ids if wid in wolf_genetics]
    if not found_ids:
        print("Error: Could not find any of the top wolf IDs in the ped file.")
        return

    num_wolves = len(found_ids)
    similarity_matrix = np.zeros((num_wolves, num_wolves))
    
    for i in range(num_wolves):
        for j in range(num_wolves):
            id_1 = found_ids[i]
            id_2 = found_ids[j]
            
            matches = sum(1 for a, b in zip(wolf_genetics[id_1], wolf_genetics[id_2]) if a == b and a != '0')
            total = sum(1 for a, b in zip(wolf_genetics[id_1], wolf_genetics[id_2]) if a != '0' and b != '0')
            
            similarity_matrix[i, j] = matches / total if total > 0 else 0

    plt.figure(figsize=(10, 8))
    
    plt.imshow(similarity_matrix, cmap='RdYlGn_r', interpolation='nearest')
    plt.colorbar(label='Genetic Similarity (Red = Danger, Green = Safe)')
    
    plt.xticks(range(num_wolves), found_ids, rotation=45)
    plt.yticks(range(num_wolves), found_ids)
    
    plt.title('Mexican Gray Wolf Breeding Recommendation Matrix', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Wolf ID')
    plt.ylabel('Wolf ID')
    
    plt.tight_layout()
    plt.savefig("Breeding_Matrix.png", dpi=300)
    print("Success! Heatmap saved as Breeding_Matrix.png")
    plt.show()

build_breeding_heatmap("MW.clean.ped", df)