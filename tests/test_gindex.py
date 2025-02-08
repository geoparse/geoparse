import unittest

from geoparse.gindex import cellpoint, pointcell


class TestPointCell(unittest.TestCase):
    def setUp(self):
        self.lats = [37.7749, 34.0522]
        self.lons = [-122.4194, -118.2437]

    def test_geohash(self):
        result = pointcell(self.lats, self.lons, "geohash", 6)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2(self):
        result = pointcell(self.lats, self.lons, "s2", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2_int(self):
        result = pointcell(self.lats, self.lons, "s2_int", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, int) for x in result))

    def test_h3(self):
        result = pointcell(self.lats, self.lons, "h3", 8)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            pointcell(self.lats, self.lons, "invalid", 6)


class TestCellPoint(unittest.TestCase):
    def test_geohash(self):
        result = cellpoint(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(isinstance(x, tuple) and len(x) == 2 for x in result))

    def test_h3(self):
        result = cellpoint(["8928308280fffff"], cell_type="h3")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2_int(self):
        result = cellpoint([9744573459660040192], cell_type="s2_int")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2(self):
        result = cellpoint(["89c25c"], cell_type="s2")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            cellpoint(["invalid"], cell_type="unknown")


if __name__ == "__main__":
    unittest.main()
