"""
Example unit tests for sims_TruthCatalog package
"""
import unittest
import desc.sims_truthcatalog

class sims_TruthCatalogTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.sims_truthcatalog.sims_TruthCatalog(self.message)
        self.assertEqual(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.sims_truthcatalog.sims_TruthCatalog)
        foo = desc.sims_truthcatalog.sims_TruthCatalog(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
