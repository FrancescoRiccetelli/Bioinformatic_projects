import os
import csv
import subprocess
from Bio.Blast import NCBIXML

# Directory contenente i genomi e i geni
genomes_dir = "C:/Users/Francesco/Desktop/Percorsoeccellenza/genomi/"
genes_dir = "C:/Users/Francesco/Desktop/Percorsoeccellenza/geni/"
results_dir = "C:/Users/Francesco/Desktop/Percorsoeccellenza/risultati/"
summary_file = os.path.join(results_dir, "summary_results.csv")

# Creare le directory se non esistono
os.makedirs(results_dir, exist_ok=True)

# Lista dei file dei genomi e dei geni
genome_files = [f for f in os.listdir(genomes_dir) if f.endswith('.fasta')]
gene_files = [f for f in os.listdir(genes_dir) if f.endswith('.fasta')]

# Passo 1: Convertire i genomi in database BLAST
for genome_file in genome_files:
    genome_path = os.path.join(genomes_dir, genome_file)
    db_name = os.path.splitext(genome_file)[0]
    makeblastdb_cmd = f"makeblastdb -dbtype nucl -in {genome_path} -out {os.path.join(results_dir, db_name)}"
    print(f"Esecuzione comando: {makeblastdb_cmd}")  # Debug message
    subprocess.run(makeblastdb_cmd, shell=True, check=True)
    print(f"Database {db_name} creato con successo.")  # Debug message

# Passo 2: Eseguire gli allineamenti con BLAST
for genome_file in genome_files:
    db_name = os.path.join(results_dir, os.path.splitext(genome_file)[0])
    for gene_file in gene_files:
        gene_path = os.path.join(genes_dir, gene_file)
        result_file = os.path.join(results_dir, f"{os.path.splitext(genome_file)[0]}_{os.path.splitext(gene_file)[0]}.xml")
        blastn_cmd = f"blastn -query {gene_path} -db {db_name} -outfmt 5 -out {result_file}"
        print(f"Esecuzione comando: {blastn_cmd}")  # Debug message
        subprocess.run(blastn_cmd, shell=True, check=True)
        print(f"Risultato salvato in {result_file}.")  # Debug message

# Passo 3: Analizzare i risultati e salvarli in un file CSV
def parse_blast_results(result_file, csv_writer):
    with open(result_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                identity = hsp.identities / hsp.align_length
                query_coverage = hsp.align_length / blast_record.query_length
                csv_writer.writerow([os.path.splitext(os.path.basename(result_file))[0], alignment.hit_def, f"{identity*100:.2f}", f"{query_coverage*100:.2f}"])

# Creare e aprire il file CSV per la scrittura
with open(summary_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["GENE", "GENOME", "IDENTITY", "COVERAGE"])
    print(f"File CSV creato: {summary_file}")  # Debug message
    
    # Analizzare i risultati e scriverli nel file CSV
    for result_file in os.listdir(results_dir):
        if result_file.endswith('.xml'):
            print(f"Analisi del file {result_file}")  # Debug message
            parse_blast_results(os.path.join(results_dir, result_file), writer)
            print(f"Risultato del file {result_file} aggiunto al CSV.")  # Debug message
