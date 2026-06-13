import streamlit as st
import re, requests
import plotly.graph_objects as go
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# --- ENGINE ---
KD_SCALE = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
            'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
            'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

def predict_localization(seq):
    """Rule-based heuristic model for subcellular localization tracking."""
    n_term = seq[:30]
    hydrophobic_count = sum(1 for aa in n_term if KD_SCALE.get(aa, 0) > 1.5)
    nls_count = len(re.findall(r'[KR]{4,}', seq))
    
    if hydrophobic_count >= 12:
        return "Secreted / Extracellular or Membrane-Bound (Signal Peptide detected)"
    elif nls_count > 0:
        return "Nuclear (Nuclear Localization Signal motif detected)"
    else:
        return "Cytoplasmic / Soluble Intracellular Space"

def run_advanced_analysis(sequence, name="Manual Input", db_ptms=None):
    sequence = sequence.upper().strip().replace("\n", "").replace(" ", "")
    analysed = ProteinAnalysis(sequence)

    # 1. Structural Fractions
    helix, turn, sheet = analysed.secondary_structure_fraction()

    # 2. Predicted Glycosylation Sites (Regex Scan)
    predicted_glyco = [m.start() + 1 for m in re.finditer(f"(?=N[^P][ST])", sequence)]
    
    # 3. Predict Subcellular Localization
    localization = predict_localization(sequence)

    # --- RENDER METADATA ---
    st.markdown("### 📔 PROTEIN METRICS")
    st.markdown(f"🧬 **Name:** {name}")
    st.markdown(f"📍 **Predicted Cellular Compartment:** `{localization}`")
    st.markdown(f"📏 **Length:** {len(sequence)} AA | ⚖️ **MW:** {analysed.molecular_weight():.2f} Da")
    st.markdown(f"⚡ **pI:** {analysed.isoelectric_point():.2f}")

    # --- COMPOSITION BREAKDOWN ---
    aa_counts = analysed.count_amino_acids()
    total = len(sequence)
    charged = (aa_counts['R'] + aa_counts['K'] + aa_counts['D'] + aa_counts['E']) / total * 100
    hydrophobic = (aa_counts['A'] + aa_counts['I'] + aa_counts['L'] + aa_counts['M'] + aa_counts['F'] + aa_counts['V'] + aa_counts['W']) / total * 100
    polar = (aa_counts['N'] + aa_counts['C'] + aa_counts['Q'] + aa_counts['S'] + aa_counts['T'] + aa_counts['Y']) / total * 100
    
    st.markdown("### 📊 BIOCHEMICAL COMPOSITION")
    col1, col2, col3 = st.columns(3)
    col1.metric("⚡ Charged Residues", f"{charged:.1f}%")
    col2.metric("💧 Hydrophobic Core", f"{hydrophobic:.1f}%")
    col3.metric("🧲 Polar Surface", f"{polar:.1f}%")

    st.markdown("### 🏗️ SECONDARY STRUCTURE FRACTIONS")
    st.markdown(f"🌀 **Alpha Helix:** {helix*100:.1f}% | ↩️ **Turn:** {turn*100:.1f}% | 🟦 **Beta Sheet:** {sheet*100:.1f}%")

    # --- PTM HANDLING COMPONENT ---
    st.markdown("### 🏷️ POST-TRANSLATIONAL MODIFICATIONS (PTMs)")
    
    phos_sites = []
    other_ptms = {}
    
    if db_ptms:
        for feature in db_ptms:
            f_type = feature.get('type')
            desc = feature.get('description', '')
            pos_start = feature.get('begin')
            
            if pos_start and pos_start.isdigit():
                pos = int(pos_start)
                
                # Segregate Phosphorylation
                if "phospho" in desc.lower() or f_type == "MOD_RES" and "phospho" in desc.lower():
                    phos_sites.append(pos)
                # Map other modifications
                elif f_type == "MOD_RES":
                    clean_desc = desc.split(';')[0].strip()
                    other_ptms[pos] = clean_desc

        st.write(f"🔬 **Verified Phosphorylation Sites:** {sorted(list(set(phos_sites))) if phos_sites else 'None listed'}")
        if other_ptms:
            st.write("✨ **Other Identified Modifications:**")
            st.json(other_ptms)
    else:
        st.caption("💡 Raw sequence input provided. Displaying canonical predicted motifs only.")
        st.markdown(f"🔗 **Predicted Glycosylation Sites (N-X-S/T):** {predicted_glyco if predicted_glyco else 'None'}")

    # --- GRAPHING BUILDER ---
    vals = [KD_SCALE.get(aa, 0) for aa in sequence]
    residues = list(range(1, len(vals)+1))
    
    fig = go.Figure()
    
    # Base Hydrophobicity Line
    fig.add_trace(go.Scatter(x=residues, y=vals, mode='lines', line=dict(color='#2ca02c', width=1.2), name='Hydrophobicity'))
    
    # 1. Plot Phosphorylation Sites (Your Custom #FF007F)
    if phos_sites:
        valid_phos = [p for p in phos_sites if p <= len(vals)]
        fig.add_trace(go.Scatter(
            x=valid_phos, 
            y=[vals[i-1] for i in valid_phos], 
            mode='markers', 
            marker=dict(color='#FF007F', size=10, symbol='circle'), 
            name='Phosphorylation',
            text=[f"Pos {p}: Phosphorylation" for p in valid_phos],
            hoverinfo='text+y'
        ))
    
    # 2. Plot Predicted Glycosylation (Your Custom #FF0000)
    if predicted_glyco:
        fig.add_trace(go.Scatter(
            x=predicted_glyco, 
            y=[vals[i-1] for i in predicted_glyco], 
            mode='markers', 
            marker=dict(color='#FF0000', size=9, symbol='diamond'), 
            name='Glyco Site (Pred)',
            text=[f"Pos {g}: Predicted Glycosylation" for g in predicted_glyco],
            hoverinfo='text+y'
        ))
                                 
    # 3. Dynamic Plotting for other PTMs
    if other_ptms:
        grouped_mods = {}
        for pos, mod_name in other_ptms.items():
            if pos <= len(vals):
                grouped_mods.setdefault(mod_name, []).append(pos)
        
        for mod_name, positions in grouped_mods.items():
            fig.add_trace(go.Scatter(
                x=positions,
                y=[vals[p-1] for p in positions],
                mode='markers',
                marker=dict(size=10, symbol='triangle-up'),
                name=mod_name,
                text=[f"Pos {p}: {mod_name}" for p in positions],
                hoverinfo='text+y'
            ))

    fig.add_hline(y=0, line_dash="dash", line_color="black", opacity=0.4)
    fig.update_layout(title=f"Complete Biophysical Profile & Multi-PTM Map: {name}", 
                      xaxis=dict(rangeslider=dict(visible=True), range=[0, min(200, len(sequence))]), 
                      yaxis=dict(title="Kyte-Doolittle Hydrophobicity"),
                      template="plotly_white", height=600)
    
    st.plotly_chart(fig, use_container_width=True)

# --- STREAMLIT WEB PAGE SETUP ---
st.set_page_config(page_title="Protein Analyzer Pro", layout="wide")
st.title("🧬 Protein Biophysics & Advanced PTM Portal")

USER_INPUT = st.text_input("Enter UniProt ID (e.g., P04040) or Paste AA Sequence:", value="P04040").strip()

if USER_INPUT:
    if re.fullmatch(r'[ARNDCEQGHILKMFPSTWYV\s\n]+', USER_INPUT.upper()):
        run_advanced_analysis(USER_INPUT)
    else:
        with st.spinner("🛰️ Querying biological APIs..."):
            meta_res = requests.get(f"https://www.ebi.ac.uk/proteins/api/features/{USER_INPUT}")
            seq_res = requests.get(f"https://www.uniprot.org/uniprot/{USER_INPUT}.fasta")

            if seq_res.status_code == 200:
                lines = seq_res.text.split('\n')
                header = lines[0]
                seq = "".join(lines[1:])
                
                # --- NEW PARSER LOGIC ---
                # Extracts full common name from the fasta header string cleanly
                name_match = re.search(r'_HUMAN\s+(.*?)\s+OS=', header)
                if name_match:
                    p_name = name_match.group(1)
                else:
                    # Fallback fallback if protein isn't human source material
                    name_match_alt = re.search(r'>\s*\S+\s+(.*?)\s+OS=', header)
                    p_name = name_match_alt.group(1) if name_match_alt else USER_INPUT
                
                features_data = []
                if meta_res.status_code == 200:
                    data = meta_res.json()
                    features_data = data.get('features', [])

                run_advanced_analysis(seq, p_name, features_data)
            else:
                st.error("❌ Invalid identifier or connection timed out.")
