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
    threshold_mask = df_probabilities["uncertainty_norm"] >= uncertainty_threshold

    st.sidebar.header("🗺️ Navegação na Sequência")

    total_pos = len(df_probabilities)
    pos_min = int(df_probabilities.index.min())
    pos_max = int(df_probabilities.index.max())

    pos_window = st.sidebar.slider(
        "Intervalo de Posições (Zoom):",
        min_value=pos_min,
        max_value=pos_max,
        value=(pos_min, pos_max),
        step=1,
        help="Selecione uma região específica da proteína para dar zoom nas visualizações."
    )
    start_pos, end_pos = pos_window

    window_mask = (df_probabilities.index >= start_pos) & (df_probabilities.index <= end_pos)
    final_mask = window_mask & threshold_mask

    df_probabilities_filtered = df_probabilities[final_mask]
    df_freq_filtered = df_freq[final_mask]

    st.sidebar.caption(
        f"🔍 **Janela selecionada:** Posições {start_pos} a {end_pos} "
        f"({end_pos - start_pos + 1} posições | {len(df_probabilities_filtered)} exibidas no heatmap)"
    )

    filtered_pos = len(df_probabilities_filtered)
    st.sidebar.info(
        f"📌 **Posições exibidas:** {filtered_pos} de {total_pos}"
        f" ({(filtered_pos/total_pos):.1%} do total)"
    )

    st.subheader("1. Perfil Geral de Incerteza & Panorama de Navegação")
    fig_1d = visualization.plot_entropy(df_probabilities, start_pos, end_pos, uncertainty_threshold)
    st.plotly_chart(fig_1d, use_container_width=True)

    st.subheader("2. Mapa de Calor 2D (Posição vs. Aminoácido)")
    st.caption(
        "💡 **Codificação Visual:** O Matiz (Cor) representa a frequência do aminoácido."
        " A Opacidade (Alpha) atenua posições incertas para reduzir a poluição visual."
    )
    if filtered_pos > 0:
        fig_heatmap = visualization.plot_uncertainty_heatmap(df_probabilities_filtered, df_freq_filtered, use_fade)
        st.plotly_chart(fig_heatmap, use_container_width=True)

        st.markdown("---")
        st.subheader("📌 Inspeção de Detalhes do Trecho Selecionado")

        col1, col2, col3, col4 = st.columns(4)
        if len(df_probabilities_filtered) > 0:
            unc_max_pos = df_probabilities_filtered["uncertainty_norm"].idxmax()
            unc_max_val = df_probabilities_filtered["uncertainty_norm"].max()
            unc_mean = df_probabilities_filtered["uncertainty_norm"].mean()
            conserved_pos = (df_probabilities_filtered["uncertainty_norm"] < 0.1).sum()

            col1.metric(label="Pico de Incerteza", value=f"Posição {unc_max_pos}", delta=f"{unc_max_val:.2f} norm")
            col2.metric(label="Incerteza Média", value=f"{unc_mean:.2f}", delta="no trecho", delta_color="off")
            col3.metric(
                label="Posições Altamente Conservadas",
                value=f"{conserved_pos}",
                delta=f"{(conserved_pos/len(df_probabilities_filtered)):.1%}"
            )
            col4.metric(label="Total de Posições no Foco", value=f"{len(df_probabilities_filtered)}")
    else:
        st.warning("⚠️ Nenhuma posição atinge o limiar de incerteza selecionado. Tente reduzir o valor do slider.")
else:
    st.info("Faça o upload do arquivo FASTA/Seed contendo os dados de alinhamento da proteína.", icon="ℹ️")
