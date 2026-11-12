import numpy as np
from matplotlib.figure import Figure
from matplotlib.pyplot import colormaps
from pandas import DataFrame
import plotly.express as px


def plot_entropy(df_probabilities: DataFrame) -> Figure:
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

    fig.update_traces(line_color="#2b5c8f", hovertemplate="<b>Posição:</b> %{x}<br /><b>Incerteza:</b> %{y:.2f}")
    fig.update_layout(yaxis_range=[0, 1], template="plotly_white", height=500)

    return fig


def plot_uncertainty_heatmap(df_probabilities: DataFrame, df_freq: DataFrame) -> Figure:
    columns_aa = [c for c in df_probabilities.columns if c not in ["entropy_bits", "uncertainty_norm"]]
    df_aa = df_probabilities[columns_aa]

    aminoacids = df_aa.columns.tolist()
    positions = df_freq.index.tolist()

    # Values array for the heatmap
    z_values = df_aa.T.values
    n_aa, n_pos = z_values.shape
    rgba_matrix = np.zeros((n_aa, n_pos, 4), dtype=np.uint8)

    uncertainty_pos = df_probabilities["uncertainty_norm"].values
    cmap = colormaps["viridis"]

    for j in range(n_pos):
        uncertainty = uncertainty_pos[j]
        alpha = 1.0 - (uncertainty * 0.75)

        for i in range(n_aa):
            freq = z_values[i, j]
            r, g, b, _ = cmap(freq)

            rgba_matrix[i, j, 0] = int(r * 255)
            rgba_matrix[i, j, 1] = int(g * 255)
            rgba_matrix[i, j, 2] = int(b * 255)
            rgba_matrix[i, j, 3] = int(alpha * 255)

    fig = px.imshow(
        rgba_matrix,
        x=positions,
        labels=dict(x="Posição na Sequência", y="Aminoácido / Gap", color="Frequência"),
        title="<b>Mapa de Calor 2D — Codificação por Canais Ortogonais (Cor + Transparência)</b>",
        aspect="auto"
    )

    fig.update_yaxes(
        tickmode="array",
        tickvals=list(range(len(aminoacids))),
        ticktext=aminoacids,
        showgrid=False
    )

    fig.update_layout(
        template="plotly_white",
        height=1000,
        xaxis=dict(showgrid=False),
    )

    return fig
