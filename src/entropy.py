import numpy as np
import pandas as pd
from pandas import DataFrame, Series


def calc_shannon_entropy(df_probabilities: DataFrame) -> Series:
    """
    Calculate the Shannon entropy (in bits) of each line (position) of a probabilities DataFrame.

    :param df_probabilities: DataFrame where each line is a position and each column the probability of a given
    aminoacid/gap
    """
    # Convert to NumPy matrix for vectorized operation
    p = df_probabilities.values

    # Mask zeroes to avoid log2(0); np.where replace 0 by 1 temporarily
    # As log2(1) == 0, the term p * log2(0) for p=0 will result exactly in 0 * 0 = 0.
    p_safe = np.where(p > 0, p, 1.0)


    # Apply formula: -sum(p * log2(p)) for each column (axis=1)
    entropy = -np.sum(p * np.log2(p_safe), axis=1)

    return pd.Series(entropy, index=df_probabilities.index, name="Entropy")


def normalize_entropy(entropy: Series, num_chars_alphabet=21) -> Series:
    """
    Normalize the entropy values between 0 and 1. The maximum entropy occurs when all possible aminoacids have the same
    probability h_max = log2(n), where n is the size of the alphabet, that is, 20 aminoacids + 1 gap = 21.
    """
    h_max = np.log2(num_chars_alphabet)
    return entropy / h_max
