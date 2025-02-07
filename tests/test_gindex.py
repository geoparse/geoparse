import unittest

from geoparse.gindex import pointcell


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


if __name__ == "__main__":
    unittest.main()
