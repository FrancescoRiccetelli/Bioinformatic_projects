from Bio import SeqIO, AlignIO, Entrez, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
import itertools

def download_sequenze(indirizzo_cartella):
    """Questa funzione utilizza Entrez per scaricare i file in formato genbank contenenti delle sequenze di DNA e le converte in file fasta"""
    Entrez.email = "francescoriccetelli41@gmail.com"

    id_list = ["1511219887", "1741596592", "2468368831", "2468368833", "2468368829", "2468368801", "2468368805", "2095779036", "1547152726", "2192683922", "1373738108", "985481444", "51988023","910767607", "1335031357"] 
    handle = Entrez.efetch( # Viene utilizzato il modulo Entrez di Biopython per scaricare dei file .gb dati determinati id
    db="nucleotide",
    id=id_list,
    rettype="gb",
    retmode="text"
    )
    gb_data = handle.read()
    handle.close()
    indirizzo_gb = os.path.join(indirizzo_cartella,"gyrB_sequenze.gb") # Fornisco l'indirizzo della cartella dove salvare il file .gb
    with open(indirizzo_gb, "w") as f:
        f.write(gb_data) # Scrivo il file .gb
    indirizzo_fasta = os.path.join(indirizzo_cartella, "gyrB_sequenze.fasta") #fornisco l'indirizzo della cartella ed il nome del file fasta
    SeqIO.convert(indirizzo_gb, "genbank", indirizzo_fasta, "fasta") # Utilizzo la funzione convert di SeqIO per convertire il file dal formato .gb al formato .fasta
    if indirizzo_gb: # Il file formato .gb non serve più quindi lo elimino
        os.remove(indirizzo_gb)

def analisi_fasta_DNA(input_fasta, output_txt):
    """Analizza e scrive in un file .txt le sequenze di DNA"""
    with open(output_txt, 'w') as out:
        intestazione = "ID\t" + " "*7 +" Lenght" + " ""\tGC%\tA\tT\tG\tC\tFreq_A\tFreq_T\tFreq_G\tFreq_C" # Scrive l'intestazione con una tabulazione per farlo sembrare una tabella
        out.write(intestazione  + "\n")

        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = record.seq.upper() # Verifica che le stringhe siano tutte maiuscole in caso di errore nel file in input e seleziona solo la sequenza
            length = len(seq) # Calcola la lunghezza della sequenza
            gc_percent = round(gc_fraction(seq) * 100, 2) # Calcola la percentuale di G+C con la funzione diretta di biopython

            count_A = seq.count("A") # Conta il numero di A
            count_T = seq.count("T") # Conta il numero di T
            count_G = seq.count("G") # Conta il numero di G
            count_C = seq.count("C") # Conta il numero di C

            freq_A = round(count_A / length, 2) # Calcola la frequenza relativa di A su tutta la sequenza
            freq_T = round(count_T / length, 2) # Calcola la frequenza relativa di T su tutta la sequenza
            freq_G = round(count_G / length, 2) # Calcola la frequenza relativa di G su tutta la sequenza
            freq_C = round(count_C / length, 2) # Calcola la frequenza relativa di C su tutta la sequenza

            # Crea una lista che permette di scrivere facilmente l'output
            row = [record.id, str(length), str(gc_percent)+"%", str(count_A), str(count_T), str(count_G), str(count_C), str(freq_A), str(freq_T), str(freq_G), str(freq_C)]

            out.write("\t".join(row) + "\n") # Itera la lista, tabula e unisce il contenuto di row e poi va a capo ad ogni riga

def allinea_sequenze(input_fasta, output_aln):
    """Allinea delle sequenze di DNA e le scrive in un file .aln con nome definito nella variabile output_aln fornita come parametro"""
    with open(output_aln, "w") as output:
        clustalw_cline = ClustalwCommandline("clustalw2", infile=input_fasta) # Richiama clustalW installata tramite ClustalwCommandline e gli fornisce il file di input
        stdout, stderr = clustalw_cline() # Esegue il comando Clustal e raccoglie eventuali errori
        aln_filename = input_fasta.replace(".fasta", ".aln") # Copio l'indirizzo del file di allineamento cambiando l'estensione poichè il nome è lo stesso del file di input
        alignment = AlignIO.read(aln_filename, "clustal") # Legge il file di allineamento
        output.write("\n=== Allineamento ===\n") # Scrive l'intestazione nel file .aln
        AlignIO.write(alignment, output, "clustal") # Scrive effettivamente l'allineamento nel file .aln chiamato come desiderato
        os.remove(aln_filename) # ClustalW crea in automatico un file .aln ed un file .dnd con lo stesso nome del file di input ma si è voluto salvarli con un altro nome, quindi qui rimuovo il file creato automaticamente da ClustalW
        dnd_file = input_fasta.replace(".fasta",".dnd") # Trova l'indirizzo del file .dnd
        if os.path.exists(dnd_file): # Non ho bisogno dell'albero ma solo dell'allineamento, quindi lo elimino dalla cartella
            os.remove(dnd_file)



def from_gene_to_prot_analyses(file_fasta_path, output_filename,output_clustal):
    """
    Legge un file FASTA di sequenze nucleotidiche, le traduce in proteine, 
    analizza ciascuna sequenza proteica (lunghezza, peso, pI, composizione),
    e scrive i risultati su un file di testo. Inoltre, genera un allineamento multiplo 
    delle proteine usando ClustalW e lo scrive anch'esso nel file.
    """

    # Dizionario per la conversione dei codici amminoacidici da 1 a 3 lettere
    nomi_amminoacidi = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
    }

    
    proteine = []
    records = []
        # Legge le sequenze nucleotidiche dal file FASTA e le traduce in amminoacidi
    for record in SeqIO.parse(file_fasta_path, "fasta"):
        tradotta = record.seq.translate(to_stop = True)  # traduce fino al primo codone di stop
        proteine.append(tradotta)
        records.append(SeqRecord(tradotta, id=record.id))

        
    with open(output_filename, "w") as output:
        for idx, record in enumerate(records):
            output.write(f"\n--- Sequenza {idx + 1} ---\n")
            output.write(f"{record.seq}\n")

            lunghezza = len(record.seq)
            output.write(f"Lunghezza: {lunghezza} aa\n")

            peso = molecular_weight(record.seq, seq_type="protein")
            output.write(f"Peso molecolare: {peso:.2f} Da\n")

            pi = round(ProteinAnalysis(str(record.seq)).isoelectric_point(), 2)
            output.write(f"Punto isoelettrico (pI): {pi}\n")

            # Calcola e scrive la percentuale di ciascun amminoacido nella sequenza
            for aa in nomi_amminoacidi:
                if aa in record.seq:
                    count = record.seq.count(aa)
                    perc = round((count / lunghezza) * 100, 2)
                    output.write(f"{aa} ({nomi_amminoacidi[aa]}): {perc}%\n")

    with open(output_clustal, "w") as output:
        # Scrive tutte le proteine in un file FASTA (necessario per ClustalW)
        proteine_fasta_file = file_fasta_path.replace(".fasta", "_tradotte.fasta")
        SeqIO.write(records, proteine_fasta_file, "fasta") # la funzione SeqIO.write vuole un file records ovvere un SeqRecord non un file con solo le sequenze

        # Esegue ClustalW per allineare le sequenze proteiche
        clustalw_cline = ClustalwCommandline("clustalw2", infile=proteine_fasta_file)
        stdout, stderr = clustalw_cline()  # genera .aln e .dnd

        # Legge l'allineamento risultante
        aln_filename = proteine_fasta_file.replace(".fasta", ".aln")
        alignment = AlignIO.read(aln_filename, "clustal")

        # Scrive l'allineamento nel file di output
        output.write("\n=== Allineamento ===\n")
        AlignIO.write(alignment, output, "clustal")
        os.remove(aln_filename)

def genera_albero(albero_dnd_input, nome_xml_output, colori):   

    # Se non viene fornita una lista di colori, usa una lista predefinita
    if colori is None:
        colori = ["blue", "red", "green", "purple", "fuchsia", "orange"]

    # Carica l'albero dal file in formato Newick
    albero = Phylo.read(albero_dnd_input, "newick")
    
    # Crea un ciclo infinito sugli elementi della lista dei colori
    color_cycle = itertools.cycle(colori)
    
    # Itera su tutti i cladi (nodi) dell'albero
    for clade in albero.find_clades():
        # Se il nodo è terminale (cioè una foglia, ovvero una sequenza osservata)
        if clade.is_terminal():
            # Assegna un colore dalla lista in ordine ciclico
            clade.color = next(color_cycle)
    
    # Imposta l'albero come radicato (utile per visualizzazione e salvataggio)
    albero.rooted = True

    # Riordina i rami dell'albero in modo gerarchico (ladderize)
    albero.ladderize(True)

    # Salva l'albero in formato phyloXML (include colori e altre annotazioni)
    Phylo.write(albero, nome_xml_output, "phyloxml")

    # Disegna l'albero con le lunghezze dei rami come etichette
    Phylo.draw(albero, branch_labels=lambda c: c.branch_length)


script_directory = os.path.dirname(os.path.abspath(__file__)) # Questa linea permette di risalire al percorso dove si sta eseguendo lo script, così che possa usarlo per gli indirizzi delle funzioni successive
fasta_file = os.path.join(script_directory, "gyrB_sequenze.fasta")
output_info = os.path.join(script_directory, "analisi_sequenze_DNA.txt")
output_aln = os.path.join(script_directory, "allineamento_sequenze_DNA.aln")
output_prot_analisi = os.path.join(script_directory, "analisi_sequenze_proteiche.txt")
output_ClustalW_proteine = os.path.join(script_directory, "allineamento_sequenze_proteiche.aln")
download_sequenze(script_directory)
analisi_fasta_DNA(fasta_file, output_info)
allinea_sequenze(fasta_file, output_aln)
from_gene_to_prot_analyses(fasta_file,output_prot_analisi,output_ClustalW_proteine)
sequenze_dnd = os.path.join(script_directory, "gyrB_sequenze_tradotte.dnd")
output_albero_xml = os.path.join(script_directory, "albero.xml")
genera_albero(sequenze_dnd, output_albero_xml, colori=None)