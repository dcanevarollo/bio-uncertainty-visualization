from pathlib import Path

from src import parser, entropy

alignment_file = Path(__file__).parent / "data" / "raw" / "spike.alignment.seed"
df_freq, df_probabilities = parser.parse_alignment(alignment_file)

df_probabilities["entropy_bits"] = entropy.calc_shannon_entropy(df_probabilities)
df_probabilities["uncertainty_norm"] = entropy.normalize_entropy(df_probabilities["entropy_bits"])

print(df_probabilities[["entropy_bits", "uncertainty_norm"]].head())
