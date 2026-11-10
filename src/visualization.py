from matplotlib.figure import Figure
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
    fig.update_layout(yaxis_range=[0, 1], template="plotly_white", height=400)

    return fig
