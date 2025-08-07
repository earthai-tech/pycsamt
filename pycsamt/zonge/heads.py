# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#       Author: LKouadio <etanoyau@gmail.com>
#       License: LGPL-3.0-or-later
"""
Group all header-level metadata containers and expose uniform
read/write hooks so higher-level code treats the header just like any
other AVG component.
"""
from __future__ import annotations

from typing import List, Dict, Any, Sequence

from .properties import (
    Hardware, SurveyAnnotation, SurveyConfiguration,
    Receiver, Transmitter, SkipFlag,
)
from .base   import AVGComponentBase
from .utils  import _block_to_dict, _dict_to_lines


__all__ =['Head']

class Head(AVGComponentBase):
    """Container for Hardware, Annotation, Config, Tx, Rx, Skip."""

    def __init__(self) -> None:
        super().__init__(data=None, meta={})
        self.hardware   = Hardware()
        self.annotation = SurveyAnnotation()
        self.config     = SurveyConfiguration()
        self.tx         = Transmitter(station=0, length_m=0.0)
        self.rx         = Receiver  (station=0, length_m=0.0)
        self.skip       = SkipFlag()

    def read(self, header_lines: Sequence[str]) -> None:
        """Populate sub-objects from raw AVG header text."""
        hw_block, ann_block, cfg_block, tx_block, rx_block = (
            [], [], [], [], []
        )
        cur = None
        for ln in header_lines:
            low = ln.lower()
            if low.startswith('\\ amtavg') or low.startswith('\\ astatic'):
                cur = hw_block
            elif '$survey.type' in low or '$survey.array' in low:
                cur = cfg_block
            elif '$rx.' in low:
                cur = rx_block
            elif '$tx.' in low or '$ xmtr' in low:
                cur = tx_block
            elif 'job.' in low:
                cur = ann_block
            if cur is not None:
                cur.append(ln)

        if hw_block:
            self.hardware.set(raw=hw_block)
        if ann_block:
            self.annotation.set(raw=ann_block)
        if cfg_block:
            self.config.set(raw=cfg_block)
        if tx_block:
            self.tx.set(**_block_to_dict(tx_block))
        if rx_block:
            self.rx.set(**_block_to_dict(rx_block))

    def write(self) -> List[str]:
        """Return ``key=value`` strings ready to prepend to an AVG file."""
        lines: List[str] = []
        lines += _dict_to_lines(self.hardware.to_json())
        lines += _dict_to_lines(self.annotation.to_json())
        lines += _dict_to_lines(self.config.to_json())
        lines += _dict_to_lines(self.tx.__dict__)
        lines += _dict_to_lines(self.rx.__dict__)
        return lines

    def asdict(self) -> Dict[str, Any]:
        return {
            "hardware"   : self.hardware.to_json(),
            "annotation" : self.annotation.to_json(),
            "config"     : self.config.to_json(),
            "transmitter": self.tx.__dict__,
            "receiver"   : self.rx.__dict__,
            "skip_flag"  : self.skip.code,
        }

    def __str__(self) -> str:
        return ("Head(hardware, annotation, config, "
                "tx, rx, skip_flag)")

