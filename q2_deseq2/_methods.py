# ----------------------------------------------------------------------------
# Copyright (c) 2024, Michal Ziemski.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd

from q2_deseq2._deseq2 import run_deseq2


def differential_expression_table(
    table: biom.Table,
    condition,
    test_level: str = '',
    reference_level: str = '',
    min_total_count: int = 10,
    fit_type: str = 'parametric',
    alpha: float = 0.05,
    cooks_cutoff: bool = True,
    independent_filtering: bool = True
) -> pd.DataFrame:
    run_result = run_deseq2(
        table=table,
        condition=condition,
        test_level=test_level,
        reference_level=reference_level,
        min_total_count=min_total_count,
        fit_type=fit_type,
        alpha=alpha,
        cooks_cutoff=cooks_cutoff,
        independent_filtering=independent_filtering
    )
    return run_result.results
