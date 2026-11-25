import unittest

import pandas as pd

from analysis import prepare, evaluate
from common import DATA, RESULTS


class BindingBenchmarkTests(unittest.TestCase):
    def test_frozen_assays_and_predictions_agree(self):
        reference = "".join((DATA / "ova.fasta").read_text().splitlines()[1:])
        eligible, conflicting = prepare(pd.read_csv(DATA / "binding_assays.csv"),
                                        pd.read_csv(DATA / "external_scores.csv"), reference)
        self.assertGreater(len(eligible), 30)
        self.assertGreaterEqual(conflicting, 0)
        self.assertFalse(eligible.duplicated(["allele", "peptide"]).any())
        current = evaluate(eligible)
        saved = pd.read_csv(RESULTS / "metrics.csv")
        pd.testing.assert_frame_equal(current, saved, atol=1e-6, check_dtype=False)
        self.assertGreaterEqual(current.roc_auc.notna().sum(), 6)
