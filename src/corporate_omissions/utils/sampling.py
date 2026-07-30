
import pandas as pd
from sklearn.utils import resample


class SamplingUtils:
    """
    Utilities for sample construction and balancing.
    """

    @staticmethod
    def balance_sample(df: pd.DataFrame, target_col: str = "co2_bool") -> pd.DataFrame:
        """
        Balance a binary sample by undersampling the majority class.

        Args:
            df (pd.DataFrame): Input DataFrame.
            target_col (str): Binary target column.

        Returns:
            pd.DataFrame: Balanced DataFrame.
        """
        if target_col not in df.columns:
            raise KeyError(f"Column '{target_col}' not found in DataFrame.")

        df_majority = df[df[target_col] == 0]
        df_minority = df[df[target_col] == 1]

        if len(df_minority) == 0:
            raise ValueError("Minority class is empty; cannot balance sample.")

        df_majority_undersampled = resample(
            df_majority,
            replace=False,
            n_samples=len(df_minority),
            random_state=42,
        )

        return pd.concat([df_minority, df_majority_undersampled], axis=0).reset_index(drop=True)
