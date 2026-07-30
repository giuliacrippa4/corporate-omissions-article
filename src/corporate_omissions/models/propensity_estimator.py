

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Literal, Any, Dict

import pandas as pd

from corporate_omissions.models.propensity_lm import PropensityModel  # your LM class file
from corporate_omissions.models.propensity_gam import PropensityGAM


Method = Literal["lm", "gam"]


@dataclass
class PropensityEstimator:
    """
    Unified interface to run LM or GAM propensity + tilt with identical arguments.

    Example:
      est = PropensityEstimator(method="lm")
      df_out = est.fit_predict(df, y_col="sc1_disclosed_complete", iv_spec="year_only")
    """
    method: Method = "lm"
    lm: Optional[PropensityModel] = None
    gam: Optional[PropensityGAM] = None
    model_: Optional[Any] = None  # the fitted underlying object (lm or gam)

    def fit_predict(self, df: pd.DataFrame, **kwargs) -> pd.DataFrame:
        m = self.method.lower()
        if m == "lm":
            self.lm = PropensityModel()
            out = self.lm.fit_predict(df, **kwargs)
            self.model_ = self.lm
            return out
        elif m == "gam":
            self.gam = PropensityGAM()
            out = self.gam.fit_predict(df, **kwargs)
            self.model_ = self.gam
            return out
        else:
            raise ValueError("method must be one of {'lm','gam'}")

    @property
    def estimated_params(self):
        # LM has params; GAM doesn't in the same way
        if isinstance(self.model_, PropensityModel):
            return self.model_.estimated_params
        return None

    @property
    def CI(self):
        if isinstance(self.model_, PropensityModel):
            return self.model_.CI
        return None

    @property
    def spec_(self) -> Optional[Dict]:
        if self.model_ is None:
            return None
        return getattr(self.model_, "spec_", None)
