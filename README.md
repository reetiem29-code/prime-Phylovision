# PhyloVision — Interactive Phylogenetic Tree Builder

## Run locally
1. Install Python 3.10 or newer.
2. Open a terminal in this folder.
3. Install dependencies:
   `pip install -r requirements.txt`
4. Start the app:
   `streamlit run app.py`

## What it does
- Accepts FASTA sequences by upload or paste.
- Validates sequence records.
- Uses MAFFT or Clustal Omega if installed; otherwise equal-length sequences are treated as already aligned.
- Builds a Neighbor Joining tree using Biopython.
- Displays a circular interactive tree with zoom/pan controls.
- Exports a Newick tree and aligned FASTA.

## Scientific limitations
- The included sample sequences are illustrative teaching data, not verified reference sequences.
- Bootstrap support is not calculated in this prototype. Do not interpret branch labels as bootstrap values unless genuine support values are present.
- The fallback for equal-length sequences does not perform alignment; it assumes the sequences are already aligned.
- Use appropriate homologous sequences, document their sources, inspect alignments, and verify results before using the output in a report.
- Protein input is supported at the interface level, but choose a distance model compatible with your data and validate the selected model for your analysis. For rigorous protein phylogeny, configure a suitable amino-acid distance model.
