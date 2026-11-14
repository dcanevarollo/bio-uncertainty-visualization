from io import StringIO

import streamlit as st

from src import parser, entropy, visualization

st.set_page_config(page_title="Análise Visual de Incerteza", layout="wide")
st.title("🧬 Análise Visual de Incerteza e Variabilidade em Proteínas")
st.markdown("Ferramenta de exploração de alinhamentos proteicos baseada na **Entropia de Shannon**.")

st.sidebar.header("Configurações")
alignment_file = st.sidebar.file_uploader("Carregar arquivo FASTA/Seed", type=["fasta", "seed", "txt"])

if alignment_file is not None:
    # The uploader returns a binary file-like object (UploadedFile, a BytesIO subclass);
    # Biopython's Stockholm parser needs a text handle.
    alignment_text = StringIO(alignment_file.getvalue().decode("utf-8"))
    df_freq, df_probabilities = parser.parse_alignment(alignment_text)
    df_probabilities["entropy_bits"] = entropy.calc_shannon_entropy(df_probabilities)
    df_probabilities["uncertainty_norm"] = entropy.normalize_entropy(df_probabilities["entropy_bits"])

    st.subheader("1. Perfil Geral de Incerteza (Entropia Normalizada)")
    fig_1d = visualization.plot_entropy(df_probabilities)
    st.plotly_chart(fig_1d, use_container_width=True)

    st.subheader("2. Mapa de Calor 2D (Posição vs. Aminoácido)")
    st.caption(
        "💡 **Codificação Visual:** O Matiz (Cor) representa a frequência do aminoácido."
        " A Opacidade (Alpha) atenua posições incertas para reduzir a poluição visual."
    )
    fig_heatmap = visualization.plot_uncertainty_heatmap(df_probabilities, df_freq)
    st.plotly_chart(fig_heatmap, use_container_width=True)
else:
    st.info("Faça o upload do arquivo FASTA/Seed contendo os dados de alinhamento da proteína.", icon="ℹ️")
