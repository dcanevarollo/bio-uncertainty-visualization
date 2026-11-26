import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pandas import DataFrame
from plotly.graph_objects import Figure


def plot_entropy(df_probabilities: DataFrame, start_pos: int, end_pos: int, uncertainty_threshold: float) -> Figure:
    fig = px.line(
        df_probabilities.reset_index(),
        x=df_probabilities.index,
        y="uncertainty_norm",
        labels={
            "index": "Posição na Sequência",
            "uncertainty_norm": "Incerteza Normalizada"
        },
        title="<b>Perfil Interativo de Incerteza (Entropia de Shannon)</b>"
    )

    fig.add_vrect(
        x0=start_pos,
        x1=end_pos,
        fillcolor="#2b5c8f",
        opacity=0.15,
        layer="below",
        line_width=1,
        annotation_text="Região em Foco",
        annotation_position="top left"
    )

    if uncertainty_threshold > 0.0:
        fig.add_hline(
            y=uncertainty_threshold,
            line_dash="dash",
            line_color="red",
            annotation_text=f"Limiar: {uncertainty_threshold:.2f}",
            annotation_position="top right",
        )

    fig.update_layout(
        template="plotly_white",
        height=330,
        xaxis=dict(rangeslider=dict(visible=True), type="linear")
    )

    return fig


def plot_uncertainty_heatmap(
        df_probabilities: DataFrame,
        df_freq: DataFrame,
        use_fade = True,
) -> Figure:
    columns_aa = [c for c in df_probabilities.columns if c not in ["entropy_bits", "uncertainty_norm"]]
    df_aa = df_probabilities[columns_aa]

    amino_acids: list[str] = df_aa.columns.tolist()
    positions: list[int] = df_freq.index.tolist()

    # Values array for the heatmap (amino acids x positions)
    z_values = df_aa.T.values.astype(float)
    n_aa, n_pos = z_values.shape

    # Hue channel: frequency -> viridis (RGBA floats in [0, 1])
    rgba = plt.get_cmap("viridis")(z_values)

    # Transparency channel: uncertain positions fade (broadcast over amino-acid rows)
    uncertainty_pos = df_probabilities["uncertainty_norm"].to_numpy(dtype=float)

    entropy_bits = (
        df_probabilities["entropy_bits"].to_numpy(dtype=float)
        if "entropy_bits" in df_probabilities.columns
        else uncertainty_pos
    )

    if use_fade:
        rgba[..., 3] = 1.0 - uncertainty_pos[np.newaxis, :] * 0.75
    else:
        rgba[..., 3] = 1.0

    # "rgba256" expects 0-255 for every channel, including alpha ("rgba" expects alpha in [0, 1])
    z_img = (rgba * 255).round().astype(np.uint8)

    hover_text = [
        [
            f"<b>📍 Posição:</b> {positions[j]}<br />"
            f"<b>🧬 Aminoácido / Gap:</b> {amino_acids[i]}<br />"
            f"<b>📊 Frequência (Probabilidade):</b> {z_values[i, j]:.2%}<br />"
            f"<b>⚙️ Entropia de Shannon:</b> {entropy_bits[j]:.3f} bits<br />"
            f"<b>📉 Incerteza Normalizada:</b> {uncertainty_pos[j]:.2f}"
            for j in range(n_pos)
        ]
        for i in range(n_aa)
    ]

    # go.Image places pixels at x0 + j*dx / y0 + i*dy, so columns are contiguous. After the uncertainty
    # filter, positions are NOT contiguous: the x axis uses the column index and ticks are relabeled
    # with the real sequence positions.
    tick_step = max(1, n_pos // 10)
    x_tickvals = list(range(0, n_pos, tick_step))
    x_ticktext = [str(positions[j]) for j in x_tickvals]

    fig = go.Figure(
        data=go.Image(
            z=z_img,
            colormodel="rgba256",
            x0=0,
            dx=1,
            y0=0,
            dy=1,
            hovertext=hover_text,
            hovertemplate="%{hovertext}<extra></extra>"
        )
    )

    title_mode = (
        "Canais Ortogonais (Cor + Transparência)" if use_fade
        else "Sem Transparência (Opacidade 100%)"
    )

    fig.update_layout(
        title=f"<b>Mapa de Calor 2D — {title_mode}</b>",
        xaxis_title="Posição na Sequência",
        yaxis_title="Aminoácido / Gap",
        yaxis=dict(
            tickmode="array",
            tickvals=list(range(n_aa)),
            ticktext=amino_acids,
            scaleanchor=False
        ),
        xaxis=dict(constrain="domain", tickmode="array", tickvals=x_tickvals, ticktext=x_ticktext),
        template="plotly_white",
        height=550,
        hoverlabel=dict(
            bgcolor="#ffffff",
            font_size=13,
            font_family="sans-serif",
            font_color="#575757",
            bordercolor="#2b5c8f"
        )
    )

    return fig
