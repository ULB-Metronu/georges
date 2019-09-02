"""
TODO
"""
from typing import Optional, List


class Input:
    def __init__(self, sequence: Optional[List] = None):
        self._sequence = sequence

    @property
    def sequence(self):
        return self._sequence
