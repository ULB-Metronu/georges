"""
TODO
"""
from typing import Optional, List


class Input:
    def __init__(self,
                 name: str = 'LINE',
                 sequence: Optional[List] = None):
        self._sequence = sequence
