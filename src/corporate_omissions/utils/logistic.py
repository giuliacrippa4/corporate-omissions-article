
import numpy as np


class LogisticUtils:
    """
    Logistic / logit-related mathematical utilities.
    """

    @staticmethod
    def logit_a_b(x: np.ndarray, a: float, b: float) -> np.ndarray:
        """
        Exponential logit-like transform used in simulation exercises.
        """
        return a * np.exp(-b * x)

    @staticmethod
    def logistic_func(x: np.ndarray, a: float, b: float) -> np.ndarray:
        """
        Standard logistic function with slope a and offset b.
        """
        return 1.0 / (1.0 + np.exp(-a * x + b))

    @staticmethod
    def logistic_function(x: np.ndarray, params: np.ndarray) -> np.ndarray:
        """
        Logistic function evaluated at x with parameter vector params.
        params[0] is the intercept.
        """
        return 1.0 / (1.0 + np.exp(-(params[0] + np.dot(params[1:], x))))

    @staticmethod
    def negative_log_likelihood(
        params: np.ndarray,
        x: np.ndarray,
        R: np.ndarray,
    ) -> float:
        """
        Negative log-likelihood for binary logistic regression.
        """
        p = LogisticUtils.logistic_function(x, params)
        eps = 1e-12  # numerical safety
        return -np.sum(R * np.log(p + eps) + (1 - R) * np.log(1 - p + eps))
