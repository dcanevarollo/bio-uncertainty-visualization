from collections import Counter
from pathlib import Path
from typing import TextIO, Union

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from pandas import DataFrame


def parse_alignment(alignment_file: Union[Path, TextIO], debug=False) -> tuple[DataFrame]:
    alignment = AlignIO.read(alignment_file, "stockholm")

    seq_num = len(alignment)
    pos_num = alignment.get_alignment_length()
    print(f"Number of sequences: {seq_num}")
    print(f"Alignment length: {pos_num}\n")

    if debug:
        print("Alignment:")
        record: SeqRecord
        for record in alignment:
            print(record.id) # Sequence identification (species)
            print(record.seq) # Sequence (aminoacid string + gaps '-'
        print()

    pos_data: list[dict] = []
    for i in range(pos_num):
        # Extract all characters from column i
        column = alignment[:, i]
        # Count absolute occurrences of each aminoacid/gap
        pos_data.append(Counter(column))

    # lines = sequence positions (0 to n-1); columns = aminoacids / gaps
    df_freq: DataFrame = pd.DataFrame(pos_data).fillna(0)

    # Calculate p(x) for each aminoacid present in the column (0.0 to 1.0)
    df_probabilities = df_freq.div(seq_num)

    if debug:
        print(df_freq.head())
        print(df_probabilities.head())

    return df_freq, df_probabilities

