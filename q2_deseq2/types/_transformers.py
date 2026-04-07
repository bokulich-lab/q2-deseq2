# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Lesser General Public License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_deseq2.plugin_setup import plugin
from . import DESeq2StatsFormat


@plugin.register_transformer
def _deseq2_stats_to_dataframe(ff: DESeq2StatsFormat) -> pd.DataFrame:
    return pd.read_csv(str(ff), sep="\t")


@plugin.register_transformer
def _dataframe_to_deseq2_stats(data: pd.DataFrame) -> DESeq2StatsFormat:
    ff = DESeq2StatsFormat()
    data.to_csv(str(ff), sep="\t", index=False)
    return ff
