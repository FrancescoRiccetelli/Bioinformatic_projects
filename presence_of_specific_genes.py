import os
import csv
from collections import defaultdict

# Directory contenente i risultati
results_dir = r"C:/Users/Francesco/Desktop/Percorsoeccellenza/"
summary_file = os.path.join(results_dir, "summary_results.csv")
output_file = os.path.join(results_dir, "presence_absence_profile.csv")

# Soglie
identity_threshold = 90.0
coverage_threshold = 95.0

# Leggere i risultati dal file summary_results.csv
genome_gene_presence = defaultdict(dict)
genes = set()

with open(summary_file, mode='r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        gene = row["GENE"]
        genome = row["GENOME"]
        identity = float(row["IDENTITY"])
        coverage = float(row["COVERAGE"])

        genes.add(gene)
        if identity >= identity_threshold and coverage >= coverage_threshold:
            genome_gene_presence[genome][gene] = 1
        else:
            genome_gene_presence[genome][gene] = 0

# Ordinare i geni per la creazione delle colonne
sorted_genes = sorted(genes)

# Scrivere il profilo di presenze/assenze nel file CSV
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    header = ["GENOME"] + sorted_genes
    writer.writerow(header)

    for genome in genome_gene_presence:
        row = [genome]
        for gene in sorted_genes:
            row.append(genome_gene_presence[genome].get(gene, 0))
        writer.writerow(row)

print(f"Profilo di presenze/assenze salvato in {output_file}")
