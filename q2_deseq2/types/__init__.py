# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._formats import (
    DESeq2StatsFormat,
    DESeq2RunDirectoryFormat,
    DESeq2StatsDirectoryFormat,
    DESeq2RunMetadataFormat,
)
from ._types import DESeq2Run, DESeq2Stats

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = [
    "DESeq2Run",
    "DESeq2Stats",
    "DESeq2StatsFormat",
    "DESeq2RunDirectoryFormat",
    "DESeq2RunMetadataFormat",
    "DESeq2StatsDirectoryFormat",
]
