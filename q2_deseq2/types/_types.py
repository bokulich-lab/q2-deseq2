# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data import FeatureData
from qiime2.plugin import SemanticType

DESeq2Stats = SemanticType("DESeq2Stats", variant_of=FeatureData.field["type"])
DESeq2Run = SemanticType("DESeq2Run")
