# PhyloVision — Updated dashboard

A student-friendly educational bioinformatics app with sequence input, alignment explorer, Neighbor-Joining tree view, tree interpretation, pairwise sequence-identity heatmap, CSV/Newick/FASTA exports, and a Creators & About tab.

## Creator credits
- Creator & project lead: Anureet Dutta
- Co-creators: Sreejeeta Kundu, Aishwarya Banik, Arya Ghosh, Ram Argha Nayak

## Run locally
1. Install Python 3.10 or newer from https://www.python.org/downloads/ (on Windows, enable **Add Python to PATH**).
2. Extract this folder.
3. Open a terminal in the folder containing `app.py` and `requirements.txt`.
4. Run:
   ```bash
   python -m pip install -r requirements.txt
   python -m streamlit run app.py
   ```
5. Open the local URL printed in the terminal (usually http://localhost:8501).

## Important scientific notes
- The built-in sequences are illustrative teaching data and are not verified biological references.
- The current version expects equal-length, already-aligned sequences. It displays and exports that alignment; it does not run MAFFT or Clustal Omega internally.
- The tree is constructed with Biopython Neighbor Joining using an uncorrected p-distance (identity) calculation.
- Bootstrap support is not implemented. No bootstrap values are displayed or claimed.
- Pairwise identity excludes gap and ambiguous positions from each pair's denominator. It is a descriptive similarity measure, not a confidence score.
- For protein analyses, use a validated protein-specific model/tool for scientific conclusions. This app is an educational demonstrator.
