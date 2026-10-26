import io
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import streamlit as st
import plotly.graph_objects as go
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

st.set_page_config(page_title="PhyloVision", page_icon="🧬", layout="wide")

st.markdown("""
<style>
:root { color-scheme: dark; }
.stApp { background: radial-gradient(circle at 48% 20%, #0a1d39 0%, #061326 45%, #030b18 100%); color: #e6f1ff; }
[data-testid="stSidebar"] {background: #071426; border-right: 1px solid #17345b;}
.block-container {padding-top: 1.25rem; max-width: 1600px;}
h1,h2,h3 {color:#eaf3ff;}
.small-muted {color:#9db4d5;font-size:0.88rem;}
.brand {font-size:2.0rem;font-weight:800;letter-spacing:-1px;}
.brand span {color:#7c8cff;}
.panel {background:rgba(8,24,47,.78); border:1px solid #193b65; border-radius:14px; padding:16px;}
.metric {background:#0a1c34;border:1px solid #1b3c65;border-radius:12px;padding:13px;}
div.stButton > button {background:linear-gradient(90deg,#3457ff,#6c48ff);color:white;border:0;border-radius:10px;font-weight:650;}
div.stDownloadButton > button {border-radius:9px;border:1px solid #315582;}
hr {border-color:#1b3658;}
</style>
""", unsafe_allow_html=True)

def parse_fasta(text):
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if len(records) < 3:
        raise ValueError("Please provide at least 3 sequences in FASTA format.")
    names = [r.id for r in records]
    if len(set(names)) != len(names):
        raise ValueError("Each FASTA record must have a unique identifier.")
    for r in records:
        if not str(r.seq).strip():
            raise ValueError(f"Sequence '{r.id}' is empty.")
    return records

def run_alignment(records):
    """Use MAFFT/Clustal Omega if installed; otherwise accept equal-length sequences as pre-aligned."""
    for executable, args in [
        ("mafft", ["--auto", "--quiet"]),
        ("clustalo", ["-i", "-", "-o", "-", "--outfmt=fasta", "--force"])
    ]:
        exe = shutil.which(executable)
        if not exe:
            continue
        with tempfile.TemporaryDirectory() as td:
            inp = Path(td) / "input.fasta"
            inp.write_text("".join(f">{r.id}\n{str(r.seq)}\n" for r in records))
            try:
                if executable == "mafft":
                    proc = subprocess.run([exe, *args, str(inp)], capture_output=True, text=True, timeout=90)
                else:
                    fasta = inp.read_text()
                    proc = subprocess.run([exe, *args], input=fasta, capture_output=True, text=True, timeout=90)
                if proc.returncode == 0 and proc.stdout.strip():
                    aligned = list(SeqIO.parse(io.StringIO(proc.stdout), "fasta"))
                    if len(aligned) == len(records) and len({len(r.seq) for r in aligned}) == 1:
                        return aligned, executable.upper()
            except (subprocess.TimeoutExpired, OSError):
                pass
    lengths = {len(r.seq) for r in records}
    if len(lengths) == 1:
        return records, "Input sequences (equal length; assumed pre-aligned)"
    raise RuntimeError(
        "These sequences have different lengths, but MAFFT/Clustal Omega was not found. "
        "Install MAFFT or Clustal Omega, or upload an already aligned FASTA file where all sequences have equal length."
    )

def make_tree(aligned, model):
    alignment = MultipleSeqAlignment([
        SeqRecord(Seq(str(r.seq).upper()), id=r.id, description="") for r in aligned
    ])
    calculator = DistanceCalculator(model)
    matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(matrix)
    tree.rooted = False
    return tree, matrix, alignment

def circular_tree_figure(tree, show_support=True):
    terminals = tree.get_terminals()
    n = len(terminals)
    if n < 2:
        return go.Figure()
    # assign terminal angles; internal angles are circular means of descendants
    angles = {}
    for i, tip in enumerate(terminals):
        angles[tip] = 2 * math.pi * i / n
    def set_angles(clade):
        if clade in angles:
            return angles[clade]
        child_angles = [set_angles(c) for c in clade.clades]
        x = sum(math.cos(a) for a in child_angles)
        y = sum(math.sin(a) for a in child_angles)
        angles[clade] = math.atan2(y, x) % (2 * math.pi)
        return angles[clade]
    set_angles(tree.root)
    depths = tree.depths()
    max_depth = max(depths.values()) or 1.0
    radii = {clade: 0.12 + 0.82 * depths[clade] / max_depth for clade in depths}
    xs, ys = [], []
    colors = ["#55c8ff", "#8b7bff", "#ff668f", "#ffb454", "#53d6b2", "#70e08a", "#d58bff"]
    color_for = {}
    for i, tip in enumerate(terminals):
        color_for[tip] = colors[i % len(colors)]
    # give descendants a shared clade color based on the first descendant tip
    def tip_color(clade):
        if clade in color_for:
            return color_for[clade]
        c = tip_color(clade.clades[0])
        color_for[clade] = c
        return c
    edge_traces = []
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        child_angles = [angles[c] for c in clade.clades]
        # radial spokes
        for child in clade.clades:
            a = angles[child]
            r0, r1 = radii[clade], radii[child]
            edge_traces.append(go.Scatter(
                x=[r0*math.cos(a), r1*math.cos(a)],
                y=[r0*math.sin(a), r1*math.sin(a)],
                mode="lines", line=dict(color=tip_color(child), width=2.2),
                hoverinfo="skip", showlegend=False
            ))
        # arc connecting child branches at parent radius
        if len(child_angles) >= 2:
            ordered = sorted(child_angles)
            # shortest arc through sorted angles; use unwrapped range
            a0, a1 = min(child_angles), max(child_angles)
            if a1 - a0 > math.pi:
                a0, a1 = a1, a0 + 2*math.pi
            arc = [a0 + (a1-a0)*j/35 for j in range(36)]
            edge_traces.append(go.Scatter(
                x=[radii[clade]*math.cos(a) for a in arc],
                y=[radii[clade]*math.sin(a) for a in arc],
                mode="lines", line=dict(color="#355579", width=1.4),
                hoverinfo="skip", showlegend=False
            ))
    fig = go.Figure()
    for trace in edge_traces:
        fig.add_trace(trace)
    for tip in terminals:
        a, r = angles[tip], radii[tip]
        x, y = r*math.cos(a), r*math.sin(a)
        fig.add_trace(go.Scatter(
            x=[x], y=[y], mode="markers+text",
            text=[tip.name], textposition="middle right",
            textfont=dict(size=11, color="#dceaff"),
            marker=dict(size=8, color=tip_color(tip), line=dict(width=1, color="#e7f3ff")),
            hovertemplate=f"<b>{tip.name}</b><br>Distance from root: {depths[tip]:.4f}<extra></extra>",
            showlegend=False
        ))
    # branch labels only if present in tree; NJ tree usually has no bootstrap support
    if show_support:
        for clade in tree.find_clades(order="preorder"):
            if clade.is_terminal() or clade.confidence is None:
                continue
            a, r = angles[clade], radii[clade]
            fig.add_annotation(x=(r+0.035)*math.cos(a), y=(r+0.035)*math.sin(a),
                               text=str(round(clade.confidence)), showarrow=False,
                               font=dict(size=10, color="#ffcf78"))
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=10,r=10,t=10,b=10), height=540,
        xaxis=dict(visible=False, scaleanchor="y", scaleratio=1, range=[-1.25,1.4]),
        yaxis=dict(visible=False, range=[-1.15,1.15]),
        dragmode="pan", showlegend=False,
        font=dict(family="Inter, Arial", color="#dceaff")
    )
    return fig

def newick_text(tree):
    out = io.StringIO()
    Phylo.write(tree, out, "newick")
    return out.getvalue()

SAMPLE = """>Homo_sapiens
ATGCTAGCTAGCTACGATCGATCGTACGATCGATCG
>Pan_troglodytes
ATGCTAGCTAGCTACGATCGATCGTACGATCGATCA
>Gorilla_gorilla
ATGCTAGCTAGCTACGATCGATCGTATGATCGATCA
>Pongo_abelii
ATGCTAGCTAGTTACGATCGATCGTATGATCGATCA
>Mus_musculus
TTGCCAGCTTGCTACGATAGATCGTATGATCAATCA
>Rattus_norvegicus
TTGCCAGCTTGCTACGATAGATCGTATGATCAATCT
"""

if "fasta_text" not in st.session_state:
    st.session_state.fasta_text = SAMPLE
if "analysis" not in st.session_state:
    st.session_state.analysis = None

with st.sidebar:
    st.markdown('<div class="brand">🧬 Phylo<span>Vision</span></div>', unsafe_allow_html=True)
    st.caption("FROM SEQUENCES TO EVOLUTION")
    st.markdown("---")
    page = st.radio("NAVIGATION", ["Tree Builder", "Alignment", "Analysis Results", "About"], label_visibility="visible")
    st.markdown("---")
    st.markdown("### Input sequences")
    upload = st.file_uploader("Upload FASTA", type=["fasta", "fa", "fas", "txt"])
    if upload is not None:
        st.session_state.fasta_text = upload.getvalue().decode("utf-8", errors="replace")
    st.session_state.fasta_text = st.text_area(
        "Or paste FASTA sequences", value=st.session_state.fasta_text, height=210,
        help="Each sequence must start with a header such as >sequence_name."
    )
    seq_type = st.selectbox("Sequence type", ["DNA", "Protein"])
    method = st.selectbox("Tree method", ["Neighbor Joining (NJ)"])
    model = st.selectbox("Distance model", ["identity", "blastn", "trans"])
    show_support = st.checkbox("Show branch support if available", value=True)
    st.caption("Bootstrap support is not automatically calculated by this basic version.")
    run = st.button("✦  Build Phylogenetic Tree", use_container_width=True)

st.markdown('<div class="brand">Phylo<span>Vision</span> <span style="font-size:1rem;color:#91a9cc;font-weight:400;">/ Interactive Tree Lab</span></div>', unsafe_allow_html=True)
st.markdown('<div class="small-muted">Explore evolutionary patterns. Inspect branches. Keep the evidence reproducible.</div>', unsafe_allow_html=True)
st.markdown("")

if run:
    try:
        records = parse_fasta(st.session_state.fasta_text)
        aligned, alignment_method = run_alignment(records)
        tree, matrix, alignment = make_tree(aligned, model)
        st.session_state.analysis = {
            "records": records, "aligned": aligned, "alignment": alignment,
            "tree": tree, "matrix": matrix, "alignment_method": alignment_method,
            "model": model, "seq_type": seq_type, "method": method
        }
        st.success(f"Analysis completed for {len(records)} sequences.")
    except Exception as e:
        st.error(str(e))

data = st.session_state.analysis
if data is None:
    st.info("Load the sample sequences or upload your own FASTA file, then select **Build Phylogenetic Tree**. The sample data are illustrative teaching sequences.")
    st.markdown("### Your workspace")
    a,b,c = st.columns(3)
    a.metric("Sequence input", "FASTA")
    b.metric("Tree style", "Circular / radial")
    c.metric("Export", "Newick")
    st.markdown("#### What this prototype can do")
    st.write("Validate FASTA input, align equal-length sequences or use MAFFT/Clustal Omega if installed, infer a Neighbor Joining tree, inspect a circular tree, review alignment, and export Newick.")
else:
    records = data["records"]
    aligned = data["aligned"]
    tree = data["tree"]
    alignment = data["alignment"]
    if page == "Tree Builder":
        left, center, right = st.columns([1.15, 2.4, 1.2], gap="medium")
        with left:
            st.markdown("### Dataset")
            st.metric("Sequences", len(records))
            st.metric("Alignment length", f"{alignment.get_alignment_length()} bp/aa")
            st.caption(f"Alignment: {data['alignment_method']}")
            st.markdown("### Sequences")
            for rec in records:
                st.markdown(f"🟢 `{rec.id}`")
            st.markdown("### View controls")
            st.caption("Use the chart toolbar to zoom, pan, reset and save an image.")
        with center:
            st.markdown("### ✧ Circular phylogenetic tree")
            st.caption("Tip labels identify sequences. Branch lengths reflect the selected distance calculation.")
            fig = circular_tree_figure(tree, show_support)
            st.plotly_chart(fig, use_container_width=True, config={"scrollZoom": True, "displaylogo": False, "modeBarButtonsToAdd": ["drawopenpath"]})
            st.caption("Important: this basic Neighbor Joining run does not calculate bootstrap support. Support labels appear only if the tree includes confidence values.")
        with right:
            st.markdown("### Analysis summary")
            st.success("Tree generated")
            st.markdown(f"**Method**  \n{data['method']}")
            st.markdown(f"**Distance model**  \n{data['model']}")
            st.markdown(f"**Input / alignment**  \n{data['alignment_method']}")
            st.markdown("### Key notes")
            st.write("• Clusters depend on sequence selection and alignment.")
            st.write("• A tree shows inferred relationships, not direct proof of ancestry.")
            st.write("• Use bootstrap or another support analysis before making strong claims.")
            st.download_button("↓ Download Tree (Newick)", newick_text(tree), file_name="phylovision_tree.nwk", mime="text/plain", use_container_width=True)
            st.download_button("↓ Download Alignment (FASTA)", "".join(f">{r.id}\n{r.seq}\n" for r in aligned), file_name="phylovision_alignment.fasta", mime="text/plain", use_container_width=True)
    elif page == "Alignment":
        st.markdown("### Multiple sequence alignment")
        st.caption(f"Alignment method: {data['alignment_method']}")
        st.code("".join(f">{r.id}\n{r.seq}\n" for r in aligned), language="text")
        st.write("Conserved columns can indicate shared sequence features. Gaps may represent inferred insertions or deletions, but alignment quality should be checked.")
    elif page == "Analysis Results":
        st.markdown("### Analysis results")
        st.write(f"**Sequences:** {len(records)}")
        st.write(f"**Alignment length:** {alignment.get_alignment_length()}")
        st.write(f"**Tree method:** {data['method']}")
        st.write(f"**Distance model:** {data['model']}")
        st.markdown("#### Pairwise distances")
        st.dataframe(data["matrix"].format_phylip().splitlines(), use_container_width=True)
        st.download_button("Download Newick tree", newick_text(tree), file_name="phylovision_tree.nwk")
    else:
        st.markdown("### About PhyloVision")
        st.write("PhyloVision is an educational bioinformatics dashboard connecting sequence input, multiple sequence alignment and phylogenetic tree construction.")
        st.write("It uses Biopython for sequence parsing, distance calculation and Neighbor Joining tree inference, with Plotly for interactive visualisation.")
        st.warning("This prototype does not automatically calculate bootstrap confidence values. Install MAFFT or Clustal Omega for unequal-length sequences; otherwise, the fallback only accepts equal-length sequences and assumes they are already aligned.")
st.markdown("---")
st.markdown('<div class="small-muted">PhyloVision • Educational prototype • Verify sequence sources and methods before reporting scientific conclusions.</div>', unsafe_allow_html=True)
