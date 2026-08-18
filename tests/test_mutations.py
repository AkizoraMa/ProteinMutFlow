from __future__ import annotations

import unittest

from mutflow.mutations import canonical_variant, parse_mutation, parse_variant


class MutationParsingTests(unittest.TestCase):
    def test_chain_qualified_multi_site_variant(self) -> None:
        mutations = parse_variant("A:S10T_A:D11E")
        self.assertEqual(canonical_variant(mutations), "A:S10T_A:D11E")
        self.assertEqual(len(mutations), 2)

    def test_legacy_chainless_token_remains_parseable(self) -> None:
        mutation = parse_mutation("S10A")
        self.assertIsNone(mutation.chain)
        self.assertEqual(mutation.residue_number, 10)

    def test_rejects_identity_mutation(self) -> None:
        with self.assertRaisesRegex(ValueError, "target equals WT"):
            parse_mutation("A:S10S")

    def test_rejects_duplicate_site(self) -> None:
        with self.assertRaisesRegex(ValueError, "same residue"):
            parse_variant("A:S10T_A:S10A")


if __name__ == "__main__":
    unittest.main()
