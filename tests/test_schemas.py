from __future__ import annotations

import json
import unittest
from importlib import resources

from jsonschema import Draft202012Validator


class ContractSchemaTests(unittest.TestCase):
    def test_all_packaged_schemas_are_valid_draft_2020_12(self) -> None:
        schema_root = resources.files("mutflow.schemas")
        for name in ("workflow.schema.json", "results.schema.json", "run.schema.json"):
            with self.subTest(schema=name):
                schema = json.loads(schema_root.joinpath(name).read_text(encoding="utf-8"))
                Draft202012Validator.check_schema(schema)


if __name__ == "__main__":
    unittest.main()
