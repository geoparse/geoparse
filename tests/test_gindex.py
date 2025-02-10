import unittest

from shapely.geometry import MultiPolygon, Polygon

from geoparse.gindex import SpatialIndexer  # Import the SpatialIndexer class from geoparse

indexer = SpatialIndexer()  # Create an instance of SpatialIndexer


class TestPointCell(unittest.TestCase):
    def setUp(self):
        self.lats = [37.7749, 34.0522]
        self.lons = [-122.4194, -118.2437]

    def test_geohash(self):
        result = indexer.pointcell(self.lats, self.lons, "geohash", 6)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2(self):
        result = indexer.pointcell(self.lats, self.lons, "s2", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2_int(self):
        result = indexer.pointcell(self.lats, self.lons, "s2_int", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, int) for x in result))

    def test_h3(self):
        result = indexer.pointcell(self.lats, self.lons, "h3", 8)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            indexer.pointcell(self.lats, self.lons, "invalid", 6)


class TestCellPoint(unittest.TestCase):
    def test_geohash(self):
        result = indexer.cellpoint(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(isinstance(x, tuple) and len(x) == 2 for x in result))

    def test_h3(self):
        result = indexer.cellpoint(["8928308280fffff"], cell_type="h3")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2_int(self):
        result = indexer.cellpoint([9744573459660040192], cell_type="s2_int")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2(self):
        result = indexer.cellpoint(["89c25c"], cell_type="s2")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            indexer.cellpoint(["invalid"], cell_type="unknown")


class TestPolyCell(unittest.TestCase):
    def test_polycell_geohash(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = indexer.polycell(geometries, cell_type="geohash", res=6)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_s2(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = indexer.polycell(geometries, cell_type="s2", res=10)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_h3(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = indexer.polycell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_multipolygon(self):
        geometries = [MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])])]
        result = indexer.polycell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_polycell_invalid_cell_type(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with self.assertRaises(ValueError) as context:
            indexer.polycell(geometries, cell_type="invalid", res=6)
        self.assertIn("Unsupported cell type", str(context.exception))

    def test_polycell_dump(self):
        import os
        import tempfile

        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with tempfile.TemporaryDirectory() as tmpdirname:
            result = indexer.polycell(geometries, cell_type="h3", res=9, dump=tmpdirname)
            self.assertIsNone(result)
            self.assertTrue(os.path.exists(os.path.join(tmpdirname, "h3", "9")))


if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main()
