
#LAB1-SYAZA_FARHANA.py
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


Entrez.email = "syazamohdzaki@gmail.com"

def retrieve_data(id):

    try:
        handle = Entrez.efetch(db='protein', id=protein_id, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        #handle.close()
        return record
    except Exception as e:
        st.error(f"Error retrieving data: {e}")
        return None


def get_basic_analysis(sequence):
    seq_analysis = ProteinAnalysis(str(sequence))
    analysis_results = {
        "length": len(sequence),
        "amino_acid_composition": seq_analysis.count_amino_acids(),
        "molecular_weight": seq_analysis.molecular_weight(),
        "isoelectric_point": seq_analysis.isoelectric_point()
    }
    return analysis_results


st.title('Lab 1 - Protein Sequence Analysis')
st.write("Bioinformatics II Lab1 Syaza Farhana")

protein_id = st.text_input('Enter UniProt ID', '')
retrieve = st.button('Retrieve')


if retrieve:
    if protein_id:
        record = retrieve_data(protein_id)
        if record:
            # Create two columns for layout
            col1, col2 = st.columns(2)

            # Display the protein information in the first column
            with col1:
                st.subheader("Retrieved Protein")
                st.write(f"**Name**: {record.name}")
                st.write(f"**Description**: {record.description}")
                st.write(f"**Sequence**: {record.seq}")
                
            # Perform and display basic analysis in the second column
            with col2:
                analysis = get_basic_analysis(record.seq)
                st.subheader("Basic Protein Analysis")
                st.write("**Sequence Length:**", analysis["length"])
                st.write("**Amino Acid Composition:**", analysis["amino_acid_composition"])
                st.write("**Molecular Weight:**", f"{analysis['molecular_weight']:.2f} Da")
                st.write("**Isoelectric Point:**", f"{analysis['isoelectric_point']:.2f}")
    else:
        st.warning('Please enter a UniProt ID')