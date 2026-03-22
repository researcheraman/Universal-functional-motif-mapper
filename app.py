import streamlit as st
import pandas as pd
import requests
from Bio import SeqIO
from io import StringIO
from collections import Counter
import numpy as np
import torch
from esm.pretrained import load_model_and_alphabet
from sklearn.cluster import HDBSCAN
from sklearn.metrics.pairwise import cosine_distances
import logomaker
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import tempfile
import os

st.set_page_config(page_title="UFMM Explorer v2", layout="wide")
st.title("🌍 UFMM Explorer v2 – Universal Functional Motif Mapper")
st.markdown("**Real-time conserved motifs across any taxa • ESM embeddings • Logos • CRISPR designs • PDF report**")

# Sidebar
st.sidebar.header("Enter Taxonomic Names")
bact_input = st.sidebar.text_input("Bacteria (e.g. Escherichia, Bacillus, Proteobacteria)", "Escherichia coli, Salmonella, Bacillus subtilis")
virus_input = st.sidebar.text_input("Viruses", "SARS-CoV-2, Influenza A, HIV-1")
mge_input = st.sidebar.text_input("Mobile Genetic Elements", "Escherichia plasmid, Staphylococcus plasmid")
host_input = st.sidebar.text_input("Optional: Host for virus-host motifs", "Homo sapiens")

params = st.sidebar.expander("Advanced Parameters")
WINDOW = params.slider("Motif window (aa)", 8, 20, 12)
MIN_OCCUR = params.slider("Min occurrences", 2, 10, 3)

run = st.sidebar.button("🚀 Run Full Advanced Analysis", type="primary")

if run:
    with st.spinner("Fetching real UniProt data + ESM embeddings + motif analysis..."):
        def fetch(group, limit=60):
            if not group.strip(): return []
            spp = [s.strip() for s in group.split(",")]
            seqs = []
            for sp in spp:
                url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=organism:{sp} AND keyword:KW-1194 AND reviewed:true&size={limit}"
                try:
                    r = requests.get(url, timeout=15)
                    recs = list(SeqIO.parse(StringIO(r.text), "fasta"))
                    seqs.extend([str(r.seq) for r in recs if 50 < len(r.seq) < 700 and "X" not in str(r.seq)])
                except: pass
            return seqs[:limit]

        bact_seqs = fetch(bact_input)
        virus_seqs = fetch(virus_input)
        mge_seqs = fetch(mge_input)
        host_seqs = fetch(host_input) if host_input else []

        all_seqs = bact_seqs + virus_seqs + mge_seqs + host_seqs
        sources = ["Bacteria"]*len(bact_seqs) + ["Virus"]*len(virus_seqs) + ["MGE"]*len(mge_seqs) + ["Host"]*len(host_seqs)

        # === ESM Functional Clustering (tiny model – CPU friendly) ===
        model, alphabet = load_model_and_alphabet("esm2_t6_8M_UR50D")
        model = model.eval()
        batch_converter = alphabet.get_batch_converter()
        data = [(i, s[:1022]) for i, s in enumerate(all_seqs)]
        _, _, tokens = batch_converter(data)
        with torch.no_grad():
            results = model(tokens, repr_layers=[6])
        embs = [results["representations"][6][i, 1:len(s)+1].mean(0).cpu().numpy() for i, s in enumerate(all_seqs)]
        emb_array = np.array(embs)

        dist = cosine_distances(emb_array)
        clusterer = HDBSCAN(min_cluster_size=4, metric='precomputed')
        labels = clusterer.fit_predict(dist)
        df_cluster = pd.DataFrame({"seq": all_seqs, "source": sources, "cluster": labels})

        # Motif discovery + logos + classification + CRISPR
        map_data = []
        for cl in set(labels):
            if cl == -1: continue
            sub = df_cluster[df_cluster.cluster == cl]
            if len(sub) < MIN_OCCUR: continue

            # Top conserved motif
            kmer_c = Counter()
            for seq in sub["seq"]:
                for i in range(len(seq) - WINDOW + 1):
                    kmer_c[seq[i:i+WINDOW]] += 1
            top_motif = kmer_c.most_common(1)[0][0] if kmer_c else ""

            # Refined classification
            if "GXXXXGK" in top_motif or "GK[T/S]" in top_motif:
                func = "Replication origins / ATP-binding site"
            elif sum(1 for a in top_motif if a in "KR") >= 4:
                func = "Protein binding motifs / Promoter-operator elements"
            elif sum(1 for a in top_motif if a in "DE") >= 3:
                func = "Enzymatic active-site signatures"
            elif sum(1 for a in top_motif if a in "LVIFW") >= 6:
                func = "Host–pathogen interaction / Membrane motifs"
            else:
                func = "Regulation / Origin-binding site"

            # Logo generation
            logo_fig = None
            if top_motif:
                # Simple frequency matrix for logo
                pos_counts = {pos: Counter() for pos in range(WINDOW)}
                for seq in sub["seq"]:
                    for i in range(len(seq) - WINDOW + 1):
                        if seq[i:i+WINDOW] == top_motif:
                            for j, aa in enumerate(top_motif):
                                pos_counts[j][aa] += 1
                pwm = pd.DataFrame(pos_counts).T.fillna(0)
                pwm = pwm.div(pwm.sum(axis=1), axis=0).fillna(0.05)
                fig, ax = plt.subplots(figsize=(6, 2))
                logo = logomaker.Logo(pwm, ax=ax)
                logo.ax.set_title(f"Cluster {cl} – {top_motif}")
                logo_fig = fig

            # CRISPR suggestions
            crispr_guides = []
            if top_motif:
                for i in range(3):
                    guide = top_motif[max(0, i):min(len(top_motif), i+20)] + "NGG"[:20]
                    crispr_guides.append(guide[:23])

            map_data.append({
                "Motif": top_motif,
                "Occurrences": kmer_c[top_motif],
                "Function": func,
                "Groups": ", ".join(sub["source"].unique()),
                "CRISPR gRNA suggestions": " | ".join(crispr_guides),
                "Logo": logo_fig
            })

        final_df = pd.DataFrame(map_data)

        st.success("✅ Full advanced analysis complete!")
        st.dataframe(final_df.drop(columns=["Logo"]), use_container_width=True)

        # Display logos
        st.subheader("Motif Logos")
        for i, row in final_df.iterrows():
            if row["Logo"] is not None:
                st.pyplot(row["Logo"])

        # PDF Report
        if st.button("📥 Download Full Professional PDF Report"):
            with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
                c = canvas.Canvas(tmp.name, pagesize=A4)
                c.setFont("Helvetica-Bold", 16)
                c.drawString(50, 800, "UFMM Explorer v2 – Motif-Function Map")
                c.setFont("Helvetica", 12)
                y = 750
                for i, row in final_df.iterrows():
                    c.drawString(50, y, f"Motif: {row['Motif']} | Function: {row['Function']}")
                    y -= 20
                    if row["Logo"] is not None:
                        # Save logo temporarily and embed
                        logo_path = f"logo_{i}.png"
                        row["Logo"].savefig(logo_path)
                        c.drawImage(logo_path, 50, y-80, width=400, height=80)
                        y -= 100
                        os.remove(logo_path)
                c.save()

                with open(tmp.name, "rb") as f:
                    pdf_bytes = f.read()
                st.download_button("Click here to download PDF", pdf_bytes, "UFMM_Full_Report.pdf", "application/pdf")

        st.caption("Real ESM embeddings • Logos • CRISPR designs • PDF report • Open for all researchers & students")
