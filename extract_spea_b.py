import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

def extract_speAB_region(input_folder):
    with open(output_fasta, 'w') as out_fasta:

        for filename in os.listdir(input_folder):
            if filename.endswith(".gbk"):
                print(filename)
                file_path = os.path.join(input_folder, filename)
                with open(file_path, 'r') as input_handle:
                    records = list(SeqIO.parse(input_handle, "genbank"))

                for record in records:
                    speA_location = None
                    speB_location = None

                    for feature in record.features:
                        if feature.type == "CDS":
                            gene_name = feature.qualifiers.get('gene', [''])[0].lower()
                            if gene_name == 'spea':
                                speA_location = feature.location
                            elif gene_name == 'speb':
                                speB_location = feature.location

                    if speA_location and speB_location:
                        print("Genes found")
                        start = min(speA_location.start, speB_location.start)# - 15
                        end = max(speA_location.end, speB_location.end)# + 15

                        if start < 0:  # Ensure start is not negative
                            start = 0

                        new_sequence = record.seq[start:end]
                        new_features = []

                        for feature in record.features:
                            if feature.location.start >= start and feature.location.end <= end:
                                new_start = feature.location.start - start
                                new_end = feature.location.end - start
                                new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
                                new_feature = SeqFeature(location=new_location, type=feature.type, qualifiers=feature.qualifiers)
                                new_features.append(new_feature)

                        if speA_location.strand == -1:
                            new_sequence = new_sequence.reverse_complement()
                            new_features = [
                                SeqFeature(
                                    location=FeatureLocation(
                                        len(new_sequence) - feature.location.end,
                                        len(new_sequence) - feature.location.start,
                                        strand=-feature.location.strand
                                    ),
                                    type=feature.type,
                                    qualifiers=feature.qualifiers
                                ) for feature in new_features
                            ]

                        new_record = record[start:end]
                        new_record.seq = new_sequence
                        new_record.features = new_features
                        new_record.annotations = record.annotations
                        new_record.id = record.id
                        new_record.description = f"Region between speA and speB"

                        output_filename = f"{os.path.splitext(filename)[0]}_speAB.gbk"
                        output_path = os.path.join(input_folder, output_filename)
                        with open(output_path, 'w') as output_handle:
                            SeqIO.write(new_record, output_handle, "genbank")
                        out_fasta.write(f">{new_record.id}_{filename}\n{new_record.seq}\n")
        #input("Fermo")

# Example usage:
input_folder = "/home/francesco/Desktop/Projects/Prosseda_05_2024/annotations/"
output_fasta = "/home/francesco/Desktop/Projects/Prosseda_05_2024/annotations/speA_speB_all.fasta"

extract_speAB_region(input_folder)

