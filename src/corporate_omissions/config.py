from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Columns:
    firm_id: str = "gvkey"
    year: str = "year"
    scope1: str = "scope1_emissions"
    disclosed: str = "scope1_disclosed"
    sector: str = "sector"
    market_cap: str = "market_cap"
    assets: str = "at"

COLS = Columns()

RANDOM_SEED: int = 123
WINSOR_PCTS = (0.01, 0.99)
