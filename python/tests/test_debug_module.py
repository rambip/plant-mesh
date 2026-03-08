import unittest

import tubulin.debug as debug_module


class DebugModuleTest(unittest.TestCase):
    def test_debug_module_prints_and_labels(self):
        print(debug_module)
        self.assertEqual(debug_module.label(), "tubulin.debug")
