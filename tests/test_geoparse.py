import unittest

import folium
from shapely.geometry import MultiPolygon, Polygon

from geoparse.geoparse import Karta, SpatialIndex


class TestKarta(unittest.TestCase):
    def setUp(self):
        # Example coordinates for setup (not used in this test)
        self.lats = [37.7749, 34.0522]
        self.lons = [-122.4194, -118.2437]

    def test_base_map(self):
        # Define southwest and northeast coordinates for the bounding box
        sw = [51.2652, -0.5426]  # Southwest coordinate (London, UK)
        ne = [51.7225, 0.2824]  # Northeast coordinate (London, UK)

        # Call the _base_map method from the Karta class
        karta = Karta._base_map(sw, ne)

        # Verify that the returned object is a Folium Map
        self.assertIsInstance(karta, folium.Map)

        # Verify that the map has the expected tile layers
        tile_layers = [
            layer.tile_name for layer in karta._children.values() if isinstance(layer, folium.raster_layers.TileLayer)
        ]
        expected_tile_layers = ["Light", "Dark", "Outdoors", "Satellite", "OSM"]
        for expected_layer in expected_tile_layers:
            self.assertIn(expected_layer, tile_layers)


class TestPointCell(unittest.TestCase):
    def setUp(self):
        self.lats = [37.7749, 34.0522]
        self.lons = [-122.4194, -118.2437]

    def test_geohash(self):
        result = SpatialIndex.point_cell(self.lats, self.lons, "geohash", 6)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2(self):
        result = SpatialIndex.point_cell(self.lats, self.lons, "s2", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_s2_int(self):
        result = SpatialIndex.point_cell(self.lats, self.lons, "s2_int", 10)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, int) for x in result))

    def test_h3(self):
        result = SpatialIndex.point_cell(self.lats, self.lons, "h3", 8)
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), len(self.lats))
        self.assertTrue(all(isinstance(x, str) for x in result))

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            SpatialIndex.point_cell(self.lats, self.lons, "invalid", 6)


class TestCellPoint(unittest.TestCase):
    def test_geohash(self):
        result = SpatialIndex.cell_point(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)
        self.assertTrue(all(isinstance(x, tuple) and len(x) == 2 for x in result))

    def test_h3(self):
        result = SpatialIndex.cell_point(["8928308280fffff"], cell_type="h3")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2_int(self):
        result = SpatialIndex.cell_point([9744573459660040192], cell_type="s2_int")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_s2(self):
        result = SpatialIndex.cell_point(["89c25c"], cell_type="s2")
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 1)
        self.assertTrue(isinstance(result[0], tuple) and len(result[0]) == 2)

    def test_invalid_cell_type(self):
        with self.assertRaises(ValueError):
            SpatialIndex.cell_point(["invalid"], cell_type="unknown")


class TestPolyCell(unittest.TestCase):
    def test_poly_cell_geohash(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = SpatialIndex.poly_cell(geometries, cell_type="geohash", res=6)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_poly_cell_s2(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = SpatialIndex.poly_cell(geometries, cell_type="s2", res=10)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_poly_cell_h3(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        result = SpatialIndex.poly_cell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertTrue(all(isinstance(cell, str) for cell in result))
        self.assertGreater(len(result), 0)

    def test_poly_cell_multipolygon(self):
        geometries = [MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])])]
        result = SpatialIndex.poly_cell(geometries, cell_type="h3", res=9)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)

    def test_poly_cell_invalid_cell_type(self):
        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with self.assertRaises(ValueError) as context:
            SpatialIndex.poly_cell(geometries, cell_type="invalid", res=6)
        self.assertIn("Unsupported cell type", str(context.exception))

    def test_poly_cell_dump(self):
        import os
        import tempfile

        geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]
        with tempfile.TemporaryDirectory() as tmpdirname:
            result = SpatialIndex.poly_cell(geometries, cell_type="h3", res=9, dump=tmpdirname)
            self.assertIsNone(result)
            self.assertTrue(os.path.exists(os.path.join(tmpdirname, "h3", "9")))


class TestCellPoly(unittest.TestCase):
    def test_cell_poly_geohash(self):
        cells = ["ezs42", "ezs43"]
        cell_type = "geohash"
        res, geoms = SpatialIndex.cell_poly(cells, cell_type)
        self.assertEqual(res, [5, 5])
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cell_poly_h3(self):
        cells = ["8928308280fffff", "8928308283bffff"]
        cell_type = "h3"
        res, geoms = SpatialIndex.cell_poly(cells, cell_type)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cell_poly_s2(self):
        cells = ["89c25c", "89c25d"]
        cell_type = "s2"
        res, geoms = SpatialIndex.cell_poly(cells, cell_type)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(geoms), 2)
        self.assertTrue(all(isinstance(g, Polygon) for g in geoms))

    def test_cell_poly_invalid_type(self):
        cells = ["ezs42"]
        with self.assertRaises(ValueError):
            SpatialIndex.cell_poly(cells, "invalid")


if __name__ == "__main__":
    unittest.main()


if __name__ == "__main__":
    unittest.main()
