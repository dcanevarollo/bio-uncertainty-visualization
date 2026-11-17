from io import StringIO

import streamlit as st

from src import parser, entropy, visualization

st.set_page_config(page_title="Análise Visual de Incerteza", layout="wide")
st.title("🧬 Análise Visual de Incerteza e Variabilidade em Proteínas")
st.markdown("Ferramenta de exploração de alinhamentos proteicos baseada na **Entropia de Shannon**.")

st.sidebar.header("⚙️ Configurações")
alignment_file = st.sidebar.file_uploader("Carregar arquivo FASTA/Seed", type=["fasta", "seed", "txt"])

st.sidebar.divider()

st.sidebar.header("🎨 Codificação Visual")
visualization_mode = st.sidebar.radio(
    "Metáfora de Incerteza:",
    options=[
        "Com Transparência (Canais Ortogonais)",
        "Sem Transparência (Opacidade Fixa)"
    ],
    help="Permite comparar o efeito da atenuação visual por transparência nas posições de alta entropia."
)
use_fade = visualization_mode == "Com Transparência (Canais Ortogonais)"

st.sidebar.divider()

st.sidebar.header("🔍 Filtro por Incerteza")
uncertainty_threshold = st.sidebar.slider(
    "Incerteza Mínima (Entropia Normalizada):",
    min_value=0.0,
    max_value=1.0,
    value=0.5,
    step=0.05,
    help="Exibe apenas as posições com nível de incerteza maior ou igual ao limiar escolhido."
)

if alignment_file is not None:
    # The uploader returns a binary file-like object (UploadedFile, a BytesIO subclass);
    # Biopython's Stockholm parser needs a text handle.
    alignment_text = StringIO(alignment_file.getvalue().decode("utf-8"))
    df_freq, df_probabilities = parser.parse_alignment(alignment_text)
    df_probabilities["entropy_bits"] = entropy.calc_shannon_entropy(df_probabilities)
    df_probabilities["uncertainty_norm"] = entropy.normalize_entropy(df_probabilities["entropy_bits"])

    # Apply threshold filters on DataFrame
    mask = df_probabilities["uncertainty_norm"] >= uncertainty_threshold
    df_probabilities_filtered = df_probabilities[mask]
    df_freq_filtered = df_freq[mask]

    total_pos = len(df_probabilities)
    filtered_pos = len(df_probabilities_filtered)
    st.sidebar.info(
        f"📌 **Posições exibidas:** {filtered_pos} de {total_pos}"
        f" ({(filtered_pos/total_pos):.1%} do total)"
    )

    st.subheader("1. Perfil Geral de Incerteza (Entropia Normalizada)")
    fig_1d = visualization.plot_entropy(df_probabilities)
    if uncertainty_threshold > 0.0:
        fig_1d.add_hline(
            y=uncertainty_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Limiar: {uncertainty_threshold:.2f}",
            annotation_position="top right",
        )
    st.plotly_chart(fig_1d, use_container_width=True)

    st.subheader("2. Mapa de Calor 2D (Posição vs. Aminoácido)")
    st.caption(
        "💡 **Codificação Visual:** O Matiz (Cor) representa a frequência do aminoácido."
        " A Opacidade (Alpha) atenua posições incertas para reduzir a poluição visual."
    )
    if filtered_pos > 0:
        fig_heatmap = visualization.plot_uncertainty_heatmap(df_probabilities_filtered, df_freq_filtered, use_fade)
        st.plotly_chart(fig_heatmap, use_container_width=True)
    else:
        st.warning("⚠️ Nenhuma posição atinge o limiar de incerteza selecionado. Tente reduzir o valor do slider.")
else:
    st.info("Faça o upload do arquivo FASTA/Seed contendo os dados de alinhamento da proteína.", icon="ℹ️")
