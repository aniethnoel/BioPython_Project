#!/usr/bin/env python
# coding: utf-8

# ### BioPython_Project
# 
# Name: Anieth Noel A

# In[ ]:


pip install biopython


# In[ ]:


from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# ### Gene Accession Number: AB000824.1

# **(i) Download the Gene Sequence & Save it in Fasta**

# In[ ]:


# Query of the Required Data using the Gene Accession Number
Entrez.email= 'user_email_id'
search=Entrez.efetch(db='nucleotide', id='AB000824.1', rettype='fasta', retmode='text')
gene= SeqIO.read(search,'fasta')
print(gene)


# In[ ]:


print(gene.id)
print(gene.description)
print(gene.seq)
print(len(gene.seq))


# In[ ]:


# Export the Query Gene sequence as fasta file for future use
gene_output=open('gene.fasta', 'w')
SeqIO.write(gene,gene_output,'fasta')
gene_output.close()


# In[ ]:


# Retrieving sequence only from the fasta file
for sequence in SeqIO.parse("gene.fasta",'fasta'):
  raw_seq=sequence.seq


# **(ii) Convert it to mRNA**

# In[ ]:


# Convert the retrieved sequence to mRNA
convert_mrna = raw_seq.transcribe()
print ("Gene:", raw_seq)
print ("mRNA:", convert_mrna)


# **(iii) Translate the gene to protein sequence – Vary the ORF**

# In[ ]:


# Translate the gene to protein Sequence
print(raw_seq.translate(to_stop=True))

#Vary the ORF
print(raw_seq[2:].translate(to_stop=True))
print(raw_seq[9:].translate(to_stop=True))


# **(iv) Design Forward & Reverse Primer of 30 Nucleotide for this gene**

# In[ ]:


# Design Forward and Backward primer of 30 Nucleotide
my_seq=gene.seq
print("Sequence:",my_seq)
print ("Forward Primer:", my_seq.complement() [:30])
print("Reverse Primer:",my_seq.complement() [-30:])


# **(v) From the Genbank file, get the protein id and download the protein sequence**

# In[ ]:


# Query of the Required Data using the Gene Accession Number as genbank format
genbank_search=Entrez.efetch(db='nucleotide', id='AB000824.1', rettype='gb', retmode='text')
genbank_format= SeqIO.read(genbank_search,'gb')
print(genbank_format)


# In[ ]:


len (genbank_format.features)


# In[ ]:


print (genbank_format.features [1])


# In[ ]:


# Retrieving the Protein Id and Protein Sequence from the Genbank File
protein_sequence = None
for details in genbank_format.features:
    if details.type == "CDS":
        protein_id = details.qualifiers["protein_id"][0]
        protein_sequence = details.qualifiers["translation"][0]
        break

if protein_sequence:
    print("Protein ID:", protein_id)
    print("Protein Sequence:", protein_sequence)
else:
    print("No protein id found.")
    print("No protein sequence found.")


# In[ ]:


# Export the Query Protein sequence as fasta file for future use
def fasta (output_file):
    with open(output_file, 'w') as f_out:
        sequence = protein_sequence
        prot_id = protein_id
        f_out.write(f'>{prot_id}\n')
        f_out.write(sequence + '\n')

output_file = ('protein_sequence_output.fasta')

fasta(output_file)


# In[ ]:


len (protein_sequence)


# **(vi) Extract the catalytic domain of this protein (use interproscan to get the domain info)**

# In[ ]:


from Bio import SearchIO


# In[ ]:


# Analyze the QueryResult from the InterProScan as xml format
input_file = "interproscan_output.xml"
interproscan_records = SearchIO.parse(input_file, "interproscan-xml")

for record in interproscan_records:
    print(record)


# In[ ]:


record [3]


# In[ ]:


c_domain = record [3]
print (c_domain)


# **(vii) Run BlastP of the domain sequence**

# In[ ]:


# Extract Domain Sequence from the Raw Protein Sequence
domain_sequence = protein_sequence[45:551]
print (domain_sequence)


# In[ ]:


# Export the domain sequence as fasta file for blastp
def dom_fasta (output_file):
    with open(output_file, 'w') as out:
        sequence = domain_sequence
        prot_name = "domain_sequence"
        out.write(f'>{prot_name}\n')
        out.write(sequence + '\n')

output_file = 'domain_sequence_output.fasta'

dom_fasta(output_file)


# In[ ]:


# Bio.Blast and NCBIWWW are used for BLAST Query
from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Load your query sequence
query_sequence = SeqIO.read("domain_sequence_output.fasta", format="fasta")

# Perform BLASTP
result_handle = NCBIWWW.qblast("blastp", "nr", query_sequence.seq)

# Write the results to a xml file
with open("domain_sequence_blastp_results.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

# Close the result handle
result_handle.close()


# In[ ]:


# Analyze the QueryResult from the BLASTP as xml format
input_file = "domain_sequence_blastp_results.xml"
blastp_records = SearchIO.parse(input_file, "blast-xml")

for blastprecord in blastp_records:
    print(blastprecord)


# In[ ]:


blastprecord [4]


# In[ ]:


# Extract the ID's from the QueryResult of BLASTP xml file
blast_xml_file = 'domain_sequence_blastp_results.xml'
def extract_id(blastp_xml):
    accession_numbers = []
    with open(blastp_xml, 'r') as result:
        blast_records = NCBIXML.parse(result)
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    accession = alignment.hit_id.split('|')[1]
                    accession_numbers.append(accession)
    return accession_numbers
accession_numbers = extract_id(blast_xml_file)
print(accession_numbers)


# **(viii) Get 5 homologs sequences of different organisms**

# In[ ]:


# Selecting the Specific ID's for Homologs Selection
selected_proteins = [0, 4, 5, 6, 8]
homologs_select = [accession_numbers[i] for i in selected_proteins]
print (homologs_select)


# In[ ]:


# Retrieve the Sequence of the Homologs by the selected ID's
Entrez.email= 'user_email_id'
homologs_search = Entrez.efetch (db ='protein', id = homologs_select , rettype='fasta', retmode='text')
homologs_seq = SeqIO.parse (homologs_search, 'fasta')
all_homologs_seq = [i for i in homologs_seq]
len (all_homologs_seq)


# In[ ]:


# Exporting the Retrieved Sequence as Fasta File
output = open('homologs.fasta', 'w')
for i in all_homologs_seq:
    SeqIO.write (i, output, 'fasta')
output.close ()


# **(ix) Change the accession ids to organism name**

# In[ ]:


#Create a dictionary with accession id as the key and organism name as the value
org_name={}
for i in open("homologs_headers.txt"): # This txt file is created by all the ID's and name in the BlastP Result XML file
  j= i.split("		")[0]
  k= i.split("[")[1].split("]")[0]
  if(j not in org_name):
    org_name[j]=k

print(org_name)


# In[ ]:


# Export the Fasta file by Changing the headers to Organism Name
out=open("homologs_org_names.fasta", 'w')
for i in SeqIO.parse("homologs.fasta" ,'fasta'):
  j=i
  i.id=org_name[j.id]
  i.description=""
  SeqIO.write(i,out,'fasta')

out.close()


# **(x) Do MSA using clustal Omega**

# In[ ]:


# Analyze the Alignment File from the Clustal Omega Website
from Bio import AlignIO
from io import StringIO

def read_aln_clustal_num(filename):
  with open(filename) as handle:
    lines = handle.readlines()
  alignment_lines = [line for line in lines if not line.startswith(';')]
  try:
    alignment = AlignIO.read(StringIO(''.join(alignment_lines)), "clustal")
    return alignment
  except ValueError:
    # Handle potential parsing errors
    print(f"Error parsing {filename}")
    return None

# Example usage
alignment = read_aln_clustal_num("clustalo-I20240504-163627-0753-24191988-p1m.aln-clustal_num")

print (alignment)

