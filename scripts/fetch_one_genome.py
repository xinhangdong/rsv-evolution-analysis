"""
fetch_one_genome.py

Step 1 milestone: download one RSV reference genome from NCBI and read what is
inside it. Proves I can pull a real genome from a public database and parse its
annotations with Biopython.
"""

from Bio import Entrez, SeqIO

# NCBI asks that you identify yourself when using their API.
Entrez.email = "dgary0501@gmail.com"

# NC_001781 = the RefSeq reference genome for RSV-B (human orthopneumovirus B).
# (Try NC_038235 later for the RSV-A reference.)
accession = "NC_038235"

# efetch = "fetch a record from NCBI". Ask the 'nucleotide' database for this
# accession, in GenBank format ("gb"), returned as plain text.
print(f"Downloading {accession} from NCBI ...")
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

# SeqIO.read parses that single GenBank record into a SeqRecord object,
# which bundles the sequence together with all of its annotations.
record = SeqIO.read(handle, "genbank")
handle.close()

# --- Basic information about the genome ---
print("\n--- Basic info ---")
print("ID:          ", record.id)
print("Description: ", record.description)
print("Organism:    ", record.annotations.get("organism"))
print("Molecule:    ", record.annotations.get("molecule_type"))
print("Length (bp): ", len(record.seq))

# --- The genes ---
# 'features' are annotations along the genome (genes, coding sequences, etc.).
# We keep only the gene features and print each gene's name and location.
print("\n--- Genes ---")
gene_count = 0
for feature in record.features:
    if feature.type == "gene":
        gene_count += 1
        gene_name = feature.qualifiers.get("gene", ["?"])[0]
        print(f"{gene_count:2d}. {gene_name:6s} {feature.location}")

print(f"\nTotal genes annotated: {gene_count}")
