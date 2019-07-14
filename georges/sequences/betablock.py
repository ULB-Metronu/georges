"""TODO"""
from __future__ import annotations
from typing import Optional
from dataclasses import dataclass


class BetaBlockType(type):
    """TODO"""
    pass


@dataclass
class BetaBlock(metaclass=BetaBlockType):
    """TODO"""
    beta11: float = 1.0
    alpha11: float = 0.0
    gamma11: Optional[float] = None
    beta22: float = 1.0
    alpha22: float = 0.0
    gamma22: Optional[float] = None
    disp1: float = 0.0
    disp2: float = 0.0
    disp3: float = 0.0
    disp4: float = 0.0
    emit1: float = 1E-9
    emit2: float = 1E-9
    emit3: float = 1E-9

    def __post_init__(self):
        if self.gamma11 is None:
            self.gamma11 = (1 + self.alpha11**2) / self.beta11
        if self.gamma22 is None:
            self.gamma22 = (1 + self.alpha22**2) / self.beta22

    def __getitem__(self, item):
        return getattr(self, item.lower())
