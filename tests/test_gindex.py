import unittest

from shapely.geometry import MultiPolygon, Polygon

from geoparse.gindex import CellGeom, GeomCell

# Create an instance of SpatialIndexer
geomcell = GeomCell()
cellgeom = CellGeom()


class TestPointCell(unittest.TestCase):
    def setUp(self):
        self.lats = [37.7749, 34.0522]
        self.lons = [-122.4194, -118.2437]

    def test_geohash(self):
        result = geomcell.pointcell(self.lats, self.lons, "geohash", 6)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2(self):
        result = geomcell.pointcell(self.lats, self.lons, "s2", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2_int(self):
        result = geomcell.pointcell(self.lats, self.lons, "s2_int", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, int) for x in result))

    def test_h3(self):
        result = geomcell.pointcell(self.lats, self.lons, "h3", 8)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            geomcell.pointcell(self.lats, self.lons, "invalid", 6)


class TestCellPoint(unittest.TestCase):
    def test_geohash(self):
        result = cellgeom.cellpoint(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(isinstance(x, tuple) and len(x) == 2 for x in result))

    def test_h3(self):
        result = cellgeom.cellpoint(["8928308280fffff"], cell_type="h3")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2_int(self):
        result = cellgeom.cellpoint([9744573459660040192], cell_type="s2_int")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2(self):
        result = cellgeom.cellpoint(["89c25c"], cell_type="s2")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            cellgeom.cellpoint(["invalid"], cell_type="unknown")


class TestPolyCell(unittest.TestCase):
    def test_polycell_geohash(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = geomcell.polycell(geometries, cell_type="geohash", res=6)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_s2(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = geomcell.polycell(geometries, cell_type="s2", res=10)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_h3(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = geomcell.polycell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_polycell_multipolygon(self):
        geometries = [MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])])]
        result = geomcell.polycell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_polycell_invalid_cell_type(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with self.assertRaises(ValueError) as context:
            geomcell.polycell(geometries, cell_type="invalid", res=6)
        self.assertIn("Unsupported cell type", str(context.exception))

    def test_polycell_dump(self):
        import os
        import tempfile

        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with tempfile.TemporaryDirectory() as tmpdirname:
            result = geomcell.polycell(geometries, cell_type="h3", res=9, dump=tmpdirname)
            self.assertIsNone(result)
            self.assertTrue(os.path.exists(os.path.join(tmpdirname, "h3", "9")))


class TestCellPoly(unittest.TestCase):
    def test_cellpoly_geohash(self):
        cells = ["ezs42", "ezs43"]
        cell_type = "geohash"
        res, geoms = cellgeom.cellpoly(cells, cell_type)
        self.assertEqual(res, [5, 5])
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cellpoly_h3(self):
        cells = ["8928308280fffff", "8928308283bffff"]
        cell_type = "h3"
        res, geoms = cellgeom.cellpoly(cells, cell_type)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cellpoly_s2(self):
        cells = ["89c25c", "89c25d"]
        cell_type = "s2"
        res, geoms = cellgeom.cellpoly(cells, cell_type)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cellpoly_invalid_type(self):
        cells = ["ezs42"]
        with self.assertRaises(ValueError):
            cellgeom.cellpoly(cells, "invalid")


if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main()
