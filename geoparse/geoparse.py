import collections
import json
import math
import os
from datetime import datetime
from math import atan2, cos, radians, sin, sqrt
from multiprocessing import Pool, cpu_count
from time import time
from typing import List, Optional, Tuple, Union

import geopandas as gpd
import h3
import numpy as np
import pandas as pd
import pygeohash
import pyproj
import requests
from polygon_geohasher.polygon_geohasher import geohash_to_polygon, polygon_to_geohashes
from s2 import s2
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform


class GeomUtils:
    """
    A utility class for performing various geometric operations on Shapely geometry objects.

    This class provides methods for determining UTM projections, transforming geometries between
    coordinate reference systems (CRS), calculating geometric statistics, and computing bearings
    for LineString geometries.

    Methods
    -------
    find_proj(geom: Union[Point, LineString, Polygon, MultiPolygon]) -> str
        Determines the appropriate UTM zone projection for a given geometry.

    trans_proj(geom: BaseGeometry, proj1: str, proj2: str) -> BaseGeometry
        Transforms a Shapely geometry object from one CRS to another.

    geom_stats(geom: Optional[Union[Polygon, MultiPolygon]] = None, projection=None, unit: str = "m") -> Optional[List[Union[int, float]]]
        Computes geometric statistics for a Polygon or MultiPolygon geometry.

    bearing(geom: LineString) -> tuple[int, int, str, str]
        Calculates the bearing and cardinal directions of a LineString.

    Notes
    -----
    - This class relies on the `shapely`, `pyproj`, and `math` libraries for geometric operations.
    - Ensure that input geometries are valid Shapely objects and that CRS definitions are valid
      and supported by `pyproj`.

    Examples
    --------
    >>> from shapely.geometry import Polygon, LineString
    >>> geom_utils = GeomUtils()

    >>> # Example for find_proj
    >>> polygon = Polygon([(-120, 35), (-121, 35), (-121, 36), (-120, 36), (-120, 35)])
    >>> utm_proj = geom_utils.find_proj(polygon)
    >>> print(utm_proj)
    'EPSG:32610'

    >>> # Example for trans_proj
    >>> point = Point(10, 50)
    >>> transformed_point = geom_utils.trans_proj(point, "EPSG:4326", "EPSG:32632")
    >>> print(transformed_point)
    <Point object at 0x...>

    >>> # Example for geom_stats
    >>> stats = geom_utils.geom_stats(polygon, unit="km")
    >>> print(stats)
    [1, 0, 4, 12322.539175581376, 444.0301771896464, 'EPSG:32631']

    >>> # Example for bearing
    >>> line = LineString([(0, 0), (1, 1)])
    >>> bearing_info = geom_utils.bearing(line)
    >>> print(bearing_info)
    (45, 45, 'NE', 'NE-SW')
    """

    @staticmethod
    def find_proj(geom: Union[Point, LineString, Polygon, MultiPolygon]) -> str:
        """
        Determines the appropriate UTM zone projection for a given geometry.

        Calculates the Universal Transverse Mercator (UTM) zone projection based on the centroid
        coordinates of the input geometry. The function returns the corresponding EPSG code for
        the UTM zone in which the geometry is located.

        Parameters
        ----------
        geom : Point, LineString, Polygon, or MultiPolygon
            A Shapely geometry object, which can be a Point, Polygon, or MultiPolygon.

        Returns
        -------
        str
            The EPSG code representing the UTM projection for the geometry's location. For the
            northern hemisphere, the function returns codes in the format 'EPSG:326XX'. For the
            southern hemisphere, it returns 'EPSG:327XX', where 'XX' is the UTM zone number.

        Notes
        -----
        The UTM (Universal Transverse Mercator) system divides the Earth into 60 longitudinal zones,
        each 6 degrees wide. This function uses the centroid of the input geometry to determine the
        appropriate zone and EPSG code.

        Examples
        --------
        >>> from shapely.geometry import Polygon
        >>> geom = Polygon([(-120, 35), (-121, 35), (-121, 36), (-120, 36), (-120, 35)])
        >>> find_proj(geom)
        'EPSG:32610'
        """
        if geom.geom_type != "Point":
            # If the geometry is not a Point, use its centroid
            geom = geom.centroid

        # Extract latitude and longitude from the geometry
        lat = geom.y
        lon = geom.x

        # Determine the base EPSG code depending on the hemisphere
        if lat >= 0:
            proj = "EPSG:326"  # Northern Hemisphere
        else:
            proj = "EPSG:327"  # Southern Hemisphere

        # Calculate the UTM zone number based on longitude
        utm = math.ceil(30 + lon / 6)

        # Return the complete EPSG code for the UTM projection
        return proj + str(utm)

    @staticmethod
    def trans_proj(geom: BaseGeometry, proj1: str, proj2: str) -> BaseGeometry:
        """
        Transforms a Shapely geometry object from one CRS to another.

        Uses `pyproj` to create a transformation pipeline that converts the input geometry
        from the source CRS (`proj1`) to the target CRS (`proj2`). The resulting geometry
        is returned in the new coordinate reference system.

        Parameters
        ----------
        geom : BaseGeometry
            A Shapely geometry object to be transformed. This can include Point, Polygon,
            MultiPolygon, LineString, or any other Shapely geometry type.
        proj1 : str
            The EPSG code or PROJ string representing the source CRS of the input geometry.
        proj2 : str
            The EPSG code or PROJ string representing the target CRS for the transformed geometry.

        Returns
        -------
        BaseGeometry
            The transformed Shapely geometry object in the target projection.

        Notes
        -----
        - The function requires `pyproj` and `shapely` libraries.
        - Ensure that the input and output CRS definitions are valid and supported by `pyproj`.

        Examples
        --------
        >>> from shapely.geometry import Point
        >>> geom = Point(10, 50)
        >>> trans_proj(geom, "EPSG:4326", "EPSG:32632")
        <Point object at 0x...>

        """
        # Create a transformation function using pyproj's Transformer
        project = pyproj.Transformer.from_crs(pyproj.CRS(proj1), pyproj.CRS(proj2), always_xy=True).transform

        # Apply the transformation to the geometry and return the transformed geometry
        return transform(project, geom)

    @staticmethod
    def geom_stats(
        geom: Optional[Union[Polygon, MultiPolygon]] = None, projection=None, unit: str = "m"
    ) -> Optional[List[Union[int, float]]]:
        """
        Computes geometric statistics for a Polygon or MultiPolygon geometry.

        Calculates various statistics for a given Shapely geometry, such as the number of shells (outer boundaries),
        number of holes, number of shell points, total area, and total perimeter length. If no geometry is provided,
        the function will print a usage example.

        Parameters
        ----------
        geom : Polygon or MultiPolygon, optional
            A Shapely geometry object (Polygon or MultiPolygon) for which to compute the statistics. If not provided,
            the function will print a usage example and not perform any computations. Default is None.
        projection: str, optional
            The EPSG code used for calculating perimeter and area of the geom in meters or kilometers, and square meters or square kilometers.
            If None, the UTM zone will be calculated.
        unit : str, optional
            The unit for area and length calculations. Accepts "m" for meters and "km" for kilometers. Default is "m".

        Returns
        -------
        list of int or float, optional
            A list containing the following statistics in order:
                - Number of shells (int)
                - Number of holes (int)
                - Number of shell points (int)
                - Total area (float)
                - Total perimeter length (float)
                - The projection used for calculating the area and perimeter
        If no geometry is provided, the function returns None.

        Examples
        --------
        >>> from shapely.geometry import Polygon
        >>> geom = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        >>> geom_stats(geom, unit="km")
        [1, 0, 4, 12322.539175581376, 444.0301771896464, 'EPSG:32631']
        """
        if not geom:  # Print usage help if geom is None
            print(
                "mdf[['nshells', 'nholes', 'nshell_points', 'area', 'perimeter', 'projection']] = [gutils.geom_stats(geom, unit='km') for geom in mdf.geometry]"
            )
            return

        # Identify the appropriate UTM zone if the projection is not provided.
        if not projection:
            projection = GeomUtils.find_proj(geom)
        projection = projection.upper()

        # Handle different geometry types
        if geom.geom_type == "Polygon":
            polylist = [geom]
        elif geom.geom_type == "MultiPolygon":
            polylist = list(geom.geoms)
        else:
            raise ValueError("The input geometry must be a Polygon or MultiPolygon.")

        # Initialize variables for calculating statistics
        n_shells = len(polylist)
        n_holes = n_shell_points = perimeter = area = 0

        # Iterate through each Polygon in the list to calculate statistics
        for poly in polylist:
            n_holes += len(poly.interiors)  # Count the number of holes
            n_shell_points += len(poly.exterior.coords) - 1  # Count the number of shell points
            # Transform geometry to the appropriate projection and calculate length/area
            perimeter += GeomUtils.trans_proj(poly, "EPSG:4326", projection).exterior.length
            area += GeomUtils.trans_proj(poly, "EPSG:4326", projection).area

        # Return statistics based on the specified unit
        if unit == "m":  # If unit is meters
            return [n_shells, n_holes, n_shell_points, area, perimeter, projection]
        else:  # If unit is kilometers
            return [n_shells, n_holes, n_shell_points, area / 1_000_000, perimeter / 1000, projection]

    @staticmethod
    def flatten_3d(geom: gpd.GeoSeries) -> List[Union[Polygon, MultiPolygon]]:
        """
        Flattens a GeoSeries of 3D Polygons or MultiPolygons into 2D geometries.

        This function removes the z-coordinate from each 3D geometry in the input GeoSeries,
        converting it into a 2D Polygon or MultiPolygon. The result is a list of 2D geometries.

        Parameters
        ----------
        geom : gpd.GeoSeries
            A GeoSeries containing 3D Polygons or MultiPolygons (geometries with z-coordinates).

        Returns
        -------
        List[Union[Polygon, MultiPolygon]]
            A list of 2D Polygons or MultiPolygons with the z-coordinates removed.

        Examples
        --------
        >>> gdf.geometry = flatten_3d(gdf.geometry)
            Converts all 3D geometries in the GeoSeries `gdf.geometry` to 2D geometries.

        Notes
        -----
        The function is useful when working with datasets that contain 3D geometries but
        only 2D geometries are needed for further spatial analysis or visualization.

        """
        new_geom = []
        for p in geom:
            if p.has_z:
                if p.geom_type == "Polygon":
                    lines = [xy[:2] for xy in list(p.exterior.coords)]
                    new_p = Polygon(lines)
                    new_geom.append(new_p)
                elif p.geom_type == "MultiPolygon":
                    new_multi_p = []
                    for ap in p:
                        lines = [xy[:2] for xy in list(ap.exterior.coords)]
                        new_p = Polygon(lines)
                        new_multi_p.append(new_p)
                    new_geom.append(MultiPolygon(new_multi_p))
        return new_geom

    @staticmethod
    def bearing(geom: LineString) -> tuple[int, int, str, str]:
        """
        Calculate the bearing and cardinal directions of a LineString.

        Parameters
        ----------
        geom : shapely.geometry.LineString
            The input geometry for which the bearing is calculated.

        Returns
        -------
        tuple[int, int, str, str]
            A tuple containing:
            - `bearing` (int): The bearing angle in degrees (0 to 359).
            - `axis_bearing` (int): The axis bearing in degrees (0 to 179).
            - `cardinal_direction` (str): The cardinal direction (e.g., "N", "NE").
            - `axis_direction` (str): The cardinal axis direction (e.g., "N-S", "E-W").
            Returns (-1, -1, "LOOP", "LOOP") if the LineString forms a loop.

        Notes
        -----
        The bearing is calculated based on the angle between the starting and
        ending points of the LineString, relative to North.

        Examples
        --------
        >>> from shapely.geometry import LineString
        >>> geom = LineString([(0, 0), (1, 1)])
        >>> bearing(geom)
        (45, 45, 'NE', 'NE-SW')

        >>> loop_geom = LineString([(0, 0), (1, 1), (0, 0)])
        >>> bearing(loop_geom)
        (-1, -1, 'LOOP', 'LOOP')
        """
        cardinal_dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        axis_dirs = ["N-S", "NE-SW", "E-W", "NW-SE"]

        # Check if the LineString is a loop
        if geom.coords[0] == geom.coords[-1]:
            return -1, -1, "LOOP", "LOOP"
        else:
            x0, y0 = geom.coords[0]
            x1, y1 = geom.coords[-1]

            angle = math.atan2(y1 - y0, x1 - x0) * 180 / math.pi
            bearing = (450 - angle) % 360  # Angle between North and the road segment

            ix = round(bearing / 45)
            return round(bearing), round(bearing) % 180, cardinal_dirs[ix % 8], axis_dirs[ix % 4]

    @staticmethod
    def line_to_points(row: gpd.GeoSeries) -> gpd.GeoDataFrame:
        """
        Splits a LineString geometry into individual Point geometries while preserving original attributes.

        This function takes a GeoSeries representing a single row of a GeoDataFrame, extracts the coordinates
        from a LineString geometry, and creates a new GeoDataFrame with each Point as a separate row. All original
        attributes from the input row are preserved in the new GeoDataFrame.

        Parameters
        ----------
        row : gpd.GeoSeries
            A GeoSeries representing a single row of a GeoDataFrame. It must include a 'geometry' column
            containing a LineString geometry.

        Returns
        -------
        gpd.GeoDataFrame
            A new GeoDataFrame where each row corresponds to a Point geometry derived from the coordinates of the LineString.
            All other columns from the original row are preserved.

        Examples
        --------
        >>> line_gdf = gpd.GeoDataFrame({"geometry": [LineString([(0, 0), (1, 1), (2, 2)])]})
        >>> point_gdf = line_to_points(line_gdf.iloc[0])
        >>> print(point_gdf)
           geometry
        0  POINT (0 0)
        1  POINT (1 1)
        2  POINT (2 2)
        """
        points = [Point(x) for x in list(row["geometry"].coords)]  # create list of Point objects
        gdf = gpd.GeoDataFrame(
            index=range(len(points)), columns=row.index
        )  # create new GeoDataFrame with all columns and Point geometry
        gdf.loc[:, "geometry"] = points
        gdf.loc[:, row.index.drop("geometry")] = row[row.index.drop("geometry")].values
        return gdf


class CellUtils:
    """
    A utility class for performing operations on spatial cells, such as compaction, uncompaction, and statistical analysis.
    It supports various spatial indexing systems, including Geohash, H3, and S2, and provides methods to manipulate and
    analyze spatial cells efficiently.

    The class provides methods to:
    1. Compact spatial cells by merging adjacent cells into parent cells, reducing the total number of cells.
    2. Uncompact S2 cells by expanding them into their child cells at a specified resolution.
    3. Compute statistics for H3 cells covering a given geometry, including the number of cells and their area.

    Supported cell types:
    - Geohash: A hierarchical spatial indexing system using base-32 encoding.
    - H3: A hexagonal hierarchical spatial indexing system.
    - S2: A spherical geometry library for spatial indexing on a sphere.

    Key Features:
    - Compaction of spatial cells to reduce redundancy and improve efficiency.
    - Uncompaction of S2 cells for detailed spatial analysis.
    - Statistical analysis of H3 cells, including cell count and area calculations.

    Examples
    --------
    >>> cellop = CellOperation()

    >>> # Compact H3 cells
    >>> h3_cells = ["8928308280fffff", "8928308280bffff"]
    >>> compacted_cells = cellop.compact_cells(h3_cells, "h3")
    >>> print(compacted_cells)

    >>> # Uncompact S2 cells
    >>> s2_tokens = ["89c2847c", "89c2847d"]
    >>> uncompacted_cells = cellop.uncompact_s2(s2_tokens, level=10)
    >>> print(uncompacted_cells)

    >>> # Compute H3 cell statistics for a geometry
    >>> from shapely.geometry import Polygon
    >>> geom = Polygon([(-122.0, 37.0), (-122.0, 38.0), (-121.0, 38.0), (-121.0, 37.0), (-122.0, 37.0)])
    >>> cell_count, cell_area = cellop.h3_stats(geom, h3_res=9, compact=True)
    >>> print(f"Cell count: {cell_count}, Cell area: {cell_area} km^2")

    Notes
    -----
    - Compaction is useful for reducing the number of cells while maintaining spatial coverage.
    - Uncompaction allows for detailed analysis by expanding cells into finer resolutions.
    - H3 cell statistics are useful for understanding the spatial distribution and coverage of a geometry.
    """

    @staticmethod
    def compact_cells(cells: list, cell_type: str) -> list:
        """
        Compacts a list of spatial cells (e.g., Geohash, S2, or H3) by merging adjacent cells into parent cells.

        The function takes a list of spatial cells and compacts them into larger cells if possible, reducing the total number
        of cells by merging adjacent cells into their parent cell at a coarser resolution. The compaction process differs based
        on the specified `cell_type` and its respective hierarchy.

        Parameters
        ----------
        cells : list
            A list of spatial cells represented as strings. Each cell corresponds to a spatial area in a specific grid system
            (e.g., Geohash, H3, or S2).

        cell_type : str
            The type of spatial cell system used. Accepted values are:
            - "geohash" : Geohash spatial indexing system.
            - "h3"      : H3 hexagonal spatial indexing system.
            - "s2"      : S2 spherical spatial indexing system.

        Returns
        -------
        list
            A list of compacted spatial cells. Each cell is represented as a string and is at the coarsest resolution possible
            based on the input cells.

        Raises
        ------
        ValueError
            If `cell_type` is not one of "geohash", "h3", or "s2".

        Notes
        -----
        - For `h3`, the function uses the built-in `h3.compact()` method.
        - For `s2`, the compaction merges cells up to their parent cells by considering the S2 hierarchy.
        - For `geohash`, cells are merged based on shared prefixes.
        """
        if cell_type == "h3":
            return h3.compact_cells(cells)
        elif cell_type == "s2":
            # Convert S2 cell IDs from tokens
            cells = [s2.CellId.from_token(item) for item in cells]
            res = cells[0].level()  # Assuming all S2 cells have the same resolution
            num_children = 4
        elif cell_type == "geohash":
            res = len(cells[0])  # Resolution is based on the length of Geohash strings
            num_children = 32
        else:
            raise ValueError(f"Invalid cell_type '{cell_type}'. Accepted values are: 'geohash', 'h3', 's2'.")

        # Initialize list to store compacted cells
        compact_cells = []
        for i in range(res, 0, -1):
            # Get parent cell IDs based on the type
            parent_ids = [cell.parent() if cell_type == "s2" else cell[: i - 1] for cell in cells]
            count_dict = collections.Counter(parent_ids)  # Count occurrences of each parent cell

            # Get indices of parent cells with the required number of children
            idx = [i for i, item in enumerate(parent_ids) if count_dict.get(item, 0) == num_children]

            # Create a mask to exclude compacted cells
            mask = [True] * len(cells)
            for ix in idx:
                mask[ix] = False
            cells = [item for i, item in enumerate(cells) if mask[i]]

            # Append compacted cells to the result
            compact_cells += cells
            cells = list(set([item for item in parent_ids if count_dict.get(item, 0) == num_children]))

        # Include any remaining cells in the compacted list
        compact_cells += cells

        if cell_type == "geohash":
            return compact_cells
        else:  # Convert S2 cells back to tokens
            return [item.to_token() for item in compact_cells]

    @staticmethod
    def uncompact_s2(compact_tokens: list, level: int) -> list:
        """
        Expands a list of compacted S2 cell tokens to a specified resolution level.

        This function takes a list of compact S2 cell tokens and generates their child cells up to the desired
        resolution level. It is used to "uncompact" S2 cells that have been previously compacted, producing a
        more detailed representation.

        Parameters
        ----------
        compact_tokens : list
            A list of S2 cell tokens represented as strings. These tokens are at a coarser resolution level and
            will be expanded into their child cells.

        level : int
            The target S2 cell resolution level to which the input tokens should be expanded. The resolution level
            determines the size of the child cells. A higher level corresponds to finer granularity (smaller cells).

        Returns
        -------
        list
            A list of S2 cell tokens represented as strings. Each token corresponds to a child cell of the input
            compact tokens, expanded to the specified resolution level.

        Raises
        ------
        ValueError
            If the provided `level` is less than or equal to the resolution level of the input `compact_tokens`.

        Example
        -------
        >>> compact_tokens = ["89c2847c", "89c2847d"]
        >>> uncompact_s2(compact_tokens, level=10)
        ["89c2847c1", "89c2847c2", "89c2847c3", ..., "89c2847d1", "89c2847d2", ...]
        """
        uncompact_tokens = []
        for token in compact_tokens:
            cell_id = s2.CellId.from_token(token)  # Convert each token to an S2 CellId object
            uncompact_tokens += list(cell_id.children(level))  # Generate child cells at the specified level
        # Convert each CellId object back to a token and remove duplicates
        uncompact_tokens = [item.to_token() for item in uncompact_tokens]
        return list(set(uncompact_tokens))

    @staticmethod
    def h3_stats(geom: BaseGeometry, h3_res: int, compact: bool = False) -> Tuple[int, float]:
        """
        Computes H3 cell statistics for a given geometry at a specified resolution.

        This function takes a Shapely geometry object and computes the number of H3 cells covering the geometry at a
        specified resolution. It also calculates the area of each H3 cell at the given resolution. Optionally, the function
        can return the compacted set of H3 cells, reducing the number of cells required to represent the geometry.

        Parameters
        ----------
        geom : shapely.geometry.base.BaseGeometry
            A Shapely geometry object (e.g., Polygon or MultiPolygon) representing the area of interest.
        h3_res : int
            The H3 resolution level for generating spatial cells. The resolution level controls the granularity of the cells.
        compact : bool, optional
            If True, the function returns a compacted set of H3 cells, reducing the number of cells needed to represent the geometry.
            Default is False.

        Returns
        -------
        tuple
            A tuple containing:
            - int: Number of H3 cells covering the given geometry.
            - float: Area of each H3 cell at the specified resolution, in square kilometers.

        Examples
        --------
        >>> from shapely.geometry import Polygon
        >>> geom = Polygon([(-122.0, 37.0), (-122.0, 38.0), (-121.0, 38.0), (-121.0, 37.0), (-122.0, 37.0)])
        >>> h3_stats(geom, h3_res=9, compact=True)
        (512, 0.001)

        Notes
        -----
        The function utilizes the H3 library for generating and compacting H3 cells and for calculating cell area. The area
        is always returned in square kilometers ("km^2").
        """
        cells = SpatialIndex.polycell([geom], cell_type="h3", res=h3_res)
        area = h3.hex_area(h3_res, unit="km^2")
        if compact:
            cells = h3.compact(cells)
        return len(cells), area


class OSMUtils:
    """
    A utility class for converting OpenStreetMap (OSM) way IDs into Shapely geometries.

    This class provides methods to retrieve and convert OSM way IDs into Shapely `Polygon`
    or `LineString` objects using the Overpass API. It supports both single and multiple
    way ID conversions.

    Methods
    -------
    way_to_geom(way_id: int) -> Optional[LineString or Polygon]
        Converts a single OSM way ID into a Shapely `Polygon` or `LineString` object.

    ways_to_geom(ids: List[int]) -> List[LineString or Polygon]
        Converts a list of OSM way IDs into a list of Shapely `Polygon` or `LineString` objects.

    Notes
    -----
    - The class uses the Overpass API to fetch geometry data for OSM ways.
    - Ensure that the provided way IDs are valid and that the Overpass API is accessible.
    - The returned geometries are in WGS84 (latitude/longitude) coordinates.

    Examples
    --------
    >>> converter = WayConverter()
    >>> # Convert a single way ID
    >>> geometry = converter.way_to_geom(123456)
    >>> print(geometry)
    POLYGON ((13.3888 52.5170, 13.3976 52.5291, 13.4286 52.5232, 13.3888 52.5170))

    >>> # Convert multiple way IDs
    >>> geometries = converter.ways_to_geom([123456, 234567])
    >>> print(geometries)
    [<shapely.geometry.polygon.Polygon object at 0x...>,
     <shapely.geometry.linestring.LineString object at 0x...>]
    """

    @staticmethod
    def way_to_geom(way_id: int, url: str = "https://overpass-api.de/api/interpreter") -> Optional[LineString or Polygon]:
        """
        Converts an OSM way ID into a Shapely Polygon or LineString object.

        This function retrieves the geometry corresponding to the given OSM way ID and
        returns it as a Shapely `Polygon` or `LineString` object based on whether the way
        forms a closed loop or not.

        Parameters
        ----------
        way_id : int
            The OpenStreetMap (OSM) way ID to be retrieved.
        url : str, optional
            The URL endpoint for the Overpass API. Defaults to "https://overpass-api.de/api/interpreter".

        Returns
        -------
        shapely.geometry.Polygon or shapely.geometry.LineString
            A Shapely `Polygon` object if the way forms a closed loop, or a `LineString`
            object otherwise.

        Notes
        -----
        - The function constructs an Overpass API query using the given way ID,
          requests the geometry, and then converts it into a Shapely geometry.
        - Assumes that the Overpass API returns data in JSON format with a "geometry" attribute.

        Examples
        --------
        >>> way_id = 123456
        >>> url = "https://overpass-api.de/api/interpreter"
        >>> geometry = way_to_geom(way_id, url)
        >>> print(geometry)
        POLYGON ((13.3888 52.5170, 13.3976 52.5291, 13.4286 52.5232, 13.3888 52.5170))
        """
        query = f"[out:json][timeout:600][maxsize:4073741824];way({way_id});out geom;"
        response = requests.get(url, params={"data": query}).json()
        response = response["elements"][0]
        geom = response["geometry"]
        coords = [(node["lon"], node["lat"]) for node in geom]
        if geom[0] == geom[-1]:  # Check if the way forms a closed loop
            return Polygon(coords)
        else:
            return LineString(coords)

    @staticmethod
    def ways_to_geom(ids: List[int], url: str = "https://overpass-api.de/api/interpreter") -> List[LineString or Polygon]:
        """
        Converts an array of OpenStreetMap (OSM) way IDs into Shapely geometries.

        This function retrieves the geometries corresponding to the given OSM way IDs and
        returns a list of Shapely `LineString` or `Polygon` objects based on the geometries
        fetched from the OSM API.

        Parameters
        ----------
        ids : list of int
            A list of OSM way IDs to be retrieved.
        url : str, optional
            The URL endpoint for the Overpass API. Defaults to "https://overpass-api.de/api/interpreter".

        Returns
        -------
        list of shapely.geometry.LineString or shapely.geometry.Polygon
            A list of Shapely `LineString` or `Polygon` objects representing the geometries
            of the OSM ways. If the way forms a closed loop, it is returned as a `Polygon`;
            otherwise, it is returned as a `LineString`.

        Notes
        -----
        - The function constructs an Overpass API query using the given IDs, requests the
          geometries, and then converts them into Shapely geometries.
        - The function assumes that the Overpass API returns data in JSON format and expects
          the "geometry" attribute to contain the coordinates.

        Examples
        --------
        >>> way_ids = [123456, 234567, 345678]
        >>> url = "https://overpass-api.de/api/interpreter"
        >>> geometries = ways_to_geom(way_ids, url)
        >>> print(geometries)
        [<shapely.geometry.polygon.Polygon object at 0x...>,
         <shapely.geometry.linestring.LineString object at 0x...>]
        """
        query = "[out:json][timeout:600][maxsize:4073741824];"
        for item in ids:
            query += f"way({item});out geom;"

        response = requests.get(url, params={"data": query}).json()
        response = response["elements"]
        nodes = response[0]["geometry"]  # used later to determine if the way is a Polygon or a LineString
        ways = [item["geometry"] for item in response]

        geoms = []
        for way in ways:
            coords = [(node["lon"], node["lat"]) for node in way]
            if nodes[0] == nodes[-1]:  # in polygons the first and last items are the same
                geoms.append(Polygon(coords))
            else:
                geoms.append(LineString(coords))
        return geoms

    @staticmethod
    def decode(encoded: str) -> list:
        """
        Decodes an encoded polyline string from Valhalla into a list of coordinates.

        Valhalla routing, map-matching, and elevation services use an encoded polyline format
        to store a series of latitude and longitude coordinates as a single string. This function
        decodes the polyline into a list of coordinates with six decimal precision.

        Parameters
        ----------
        encoded : str
            An encoded polyline string as per the Valhalla encoding format.

        Returns
        -------
        list of list of float
            A list of [longitude, latitude] pairs decoded from the input polyline string.

        Notes
        -----
        - The function uses six decimal degrees of precision for decoding Valhalla's encoded polylines.
        - The decoded coordinates are returned in [longitude, latitude] format.

        References
        ----------
        - https://github.com/valhalla/valhalla-docs/blob/master/decoding.md#decode-a-route-shape

        Examples
        --------
        >>> encoded_polyline = "_p~iF~ps|U_ulLnnqC_mqNvxq`@"
        >>> decoded_coords = decode(encoded_polyline)
        >>> print(decoded_coords)
        [[-120.2, 38.5], [-120.95, 40.7], [-126.453, 43.252]]
        """
        inv = 1.0 / 1e6  # Six decimal places of precision in Valhalla
        decoded = []
        previous = [0, 0]
        i = 0
        while i < len(encoded):  # For each byte in the encoded string
            ll = [0, 0]  # To store latitude and longitude
            for j in [0, 1]:
                shift = 0
                byte = 0x20
                while byte >= 0x20:  # Keep decoding bytes until the complete coordinate is read
                    byte = ord(encoded[i]) - 63
                    i += 1
                    ll[j] |= (byte & 0x1F) << shift
                    shift += 5
                ll[j] = previous[j] + (~(ll[j] >> 1) if ll[j] & 1 else (ll[j] >> 1))
                previous[j] = ll[j]
            # Convert to float and format the result
            decoded.append([float("%.6f" % (ll[1] * inv)), float("%.6f" % (ll[0] * inv))])
        return decoded

    @staticmethod
    def map_matching(df: pd.DataFrame, cost: str, url: str, format: str = "osrm") -> Optional[dict]:
        """
        Performs map matching using Valhalla's Meili service.

        Map matching aligns a series of GPS points onto a road network. This function takes a DataFrame
        of coordinates, sends a request to the Meili map-matching service, and returns the matched
        coordinates along with other route information.

        Parameters
        ----------
        df : pd.DataFrame
            A pandas DataFrame containing the GPS coordinates to be map-matched. It should be in the
            format of [{"lon": float, "lat": float}, ...].
        cost : str
            The routing profile to use for map matching. Common values include "auto", "bicycle",
            or "pedestrian".
        url : str
            The URL endpoint for the Meili map-matching service.
        format : str, optional
            The response format for the request, either "osrm" or "geojson". Defaults to "osrm".

        Returns
        -------
        Optional[dict]
            A dictionary representing the JSON response from the map-matching service if the request
            is successful, otherwise None.

        Examples
        --------
        >>> coordinates = [{"lon": -73.9857, "lat": 40.7484}, {"lon": -73.9851, "lat": 40.7478}]
        >>> df = pd.DataFrame(coordinates)
        >>> url = "https://valhalla.mapzen.com/trace_attributes"
        >>> matched_route = map_matching(df, "auto", url)
        >>> print(matched_route)
        {'shape': '_p~iF~ps|U_ulLnnqC_mqNvxq`@', 'confidence_score': 1.0}
        """
        meili_head = '{"shape":'  # Initial portion of the request body
        meili_coordinates = df.to_json(orient="records")  # Convert DataFrame to JSON format

        meili_tail = f', "search_radius":150, "shape_match":"map_snap", "costing":"{cost}", "format":"{format}"}}'

        # Combine the header, coordinates, and tail into a single request body
        meili_request_body = meili_head + meili_coordinates + meili_tail

        # Send the request to the Meili service
        response = requests.post(url, data=meili_request_body, headers={"Content-type": "application/json"})
        if response.status_code == 200:
            return response.json()  # Convert the JSON response to a dictionary
        else:
            return None


class SpatialIndex:
    """
    A class for performing spatial indexing operations on geographic data.

    This class provides methods to convert geographic coordinates (latitude, longitude) and geometries
    (Polygon, MultiPolygon) into spatial cell representations using various encoding systems such as
    Geohash, S2, and H3. It also supports parallel processing for efficient handling of large datasets.

    Methods:
    --------
    pointcell(lats, lons, cell_type, res):
        Converts latitude and longitude coordinates into spatial index representations.

    polycell(geoms, cell_type, res, dump):
        Converts a list of geometries into a set of unique spatial cells.

    ppointcell(lats, lons, cell_type, res):
        Converts latitude and longitude coordinates into spatial index representations in parallel.

    ppolycell(mdf, cell_type, res, compact, dump, verbose):
        Performs parallelized conversion of geometries in a GeoDataFrame to cell identifiers.

    cellpoint(cells, cell_type):
        Converts a list of cell IDs into their corresponding centroids.

    cellpoly(cells, cell_type):
        Converts a list of spatial cells to their corresponding geometries and resolution levels.

    pcellpoint(cells, cell_type):
        Converts a list of cell IDs into their corresponding latitude and longitude points in parallel.

    pcellpoly(cells, cell_type):
        Parallelized version of `cellpoly`, converting a list of spatial cells to geometries and resolution levels.
    """

    @staticmethod
    def pointcell(lats: list[float], lons: list[float], cell_type: str, res: int) -> list:
        """
        Convert latitude and longitude coordinates into spatial index representations using various cell encoding types.

        Parameters:
        -----------
        lats : list[float]
            A list of latitude values (in decimal degrees) for each point to encode.
        lons : list[float]
            A list of longitude values (in decimal degrees) for each point to encode.
        cell_type : str
            The type of spatial encoding to use. Options are:
                - 'geohash': Encodes the coordinates using the geohash format.
                - 's2': Encodes the coordinates using the S2 library, outputting a string representation.
                - 's2_int': Encodes the coordinates using the S2 library, outputting an integer representation.
                - 'h3': Encodes the coordinates using the H3 library, outputting a hex string.
        res : int
            The resolution or precision level for the encoding, specific to each encoding type.

        Returns:
        --------
        list
            A list of encoded cell identifiers. The data type of each identifier depends on `cell_type`:
                - 'geohash' and 's2': list of str
                - 's2_int': list of int
                - 'h3': list of str

        Raises:
        -------
        ValueError
            If `cell_type` is not one of 'geohash', 's2', 's2_int', or 'h3'.

        Example:
        --------
        >>> lats = [37.7749, 34.0522]
        >>> lons = [-122.4194, -118.2437]
        >>> pointcell(lats, lons, "geohash", 6)
        ['9q8yy', '9qh0b']
        """
        if cell_type == "geohash":
            return [pygeohash.encode(lat, lon, res) for lat, lon in zip(lats, lons)]
        elif cell_type == "s2":
            return [s2.geo_to_s2(lat, lon, res) for lat, lon in zip(lats, lons)]  # string
        elif cell_type == "s2_int":
            return [int(s2.geo_to_s2(lat, lon, res), 16) for lat, lon in zip(lats, lons)]  # int data type requires less memory
        elif cell_type == "h3":
            return [h3.latlng_to_cell(lat, lon, res) for lat, lon in zip(lats, lons)]
        else:
            raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', 's2_int', or 'h3'.")

    @staticmethod
    def polycell(
        geoms: List[Union[Polygon, MultiPolygon]], cell_type: str, res: int, dump: str = None
    ) -> Union[List[str], None]:
        """
        Converts a list of geometries into a set of unique spatial cells based on the specified cell type and resolution.

        This function takes a list of Shapely geometries (e.g., Polygon, MultiPolygon) and converts them into spatial cells
        using one of the supported cell systems: Geohash, S2, or H3. The resulting cells are returned as a list of unique
        cell IDs. If `dump` is set to a valid directory path, the cells are saved to a file in that directory, instead of being returned.

        Parameters
        ----------
        geoms : list of shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            A list of Shapely geometry objects (Polygon or MultiPolygon).
        cell_type : str
            The type of spatial cell system to use. Supported values are "geohash", "s2", or "h3".
        res : int
            The resolution level for the spatial cells. The resolution parameter determines the granularity of the cells.
        dump : str, optional
            If set to a valid directory path (string), the cells are saved to a file in the specified folder.
            The file will be saved in a subdirectory structure following the pattern: `/path/to/dir/cell_type/res/`.
            If `dump` is None, the function returns the list of cell IDs. Default is None.

        Returns
        -------
        list of str or None
            If `dump` is None, a list of unique cell IDs is returned.
            If `dump` is provided, None is returned after saving the cells to a file.

        Raises
        ------
        ValueError
            If `cell_type` is not one of the supported values ("geohash", "s2", "h3").

        Examples
        --------
        >>> from shapely.geometry import Polygon, MultiPolygon
        >>> geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), MultiPolygon([...])]
        >>> # Convert geometries to H3 cells at resolution 9
        >>> h3_cells = polycell(geometries, cell_type="h3", res=9)

        >>> # Convert geometries to S2 cells and save to a directory
        >>> polycell(geometries, cell_type="s2", res=10, dump="~/Desktop/spatial_cells")
        """

        polys = []
        for geom in geoms:
            if geom.geom_type == "Polygon":
                polys += [geom.__geo_interface__]
            elif geom.geom_type == "MultiPolygon":  # If MultiPolygon, extract each Polygon separately
                polys += [g.__geo_interface__ for g in geom.geoms]

        if cell_type == "geohash":
            cells = list({geohash for geom in geoms for geohash in polygon_to_geohashes(geom, precision=res, inner=False)})
        elif cell_type == "s2":
            cells = list(
                {item["id"] for poly in polys for item in s2.polyfill(poly, res, geo_json_conformant=True, with_id=True)}
            )
        elif cell_type == "h3":
            cells = [cell for poly in polys for cell in h3.geo_to_cells(poly, res)]
        else:
            raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', or 'h3'.")

        if not dump:
            return cells
        else:
            # Create the directories if they don't exist
            cells_path = os.path.expanduser(f"{dump}/{cell_type}/{res}")
            os.makedirs(cells_path, exist_ok=True)
            with open(f"{cells_path}/{datetime.now()}.txt", "w") as json_file:
                json.dump(cells, json_file)
            return None

    @staticmethod
    def ppointcell(lats: list[float], lons: list[float], cell_type: str, res: int) -> list:
        """
        Converts lists of latitude and longitude points to cell identifiers in parallel.

        This function takes lists of latitude and longitude points and converts each pair
        into a cell identifier based on the specified `cell_type` and resolution `res`.
        It leverages parallel processing to speed up the conversion, dividing the data
        into chunks and using `Pool.starmap` for concurrent execution.

        Parameters
        ----------
        lats : list of float
            List of latitude values.
        lons : list of float
            List of longitude values, corresponding element-wise to the latitude list.
        cell_type : str
            The type of spatial encoding to use. Options are:
                - 'geohash': Encodes the coordinates using the geohash format.
                - 's2': Encodes the coordinates using the S2 library, outputting a string representation.
                - 's2_int': Encodes the coordinates using the S2 library, outputting an integer representation.
                - 'h3': Encodes the coordinates using the H3 library, outputting a hex string.
        res : int
            Resolution or precision level for the cell identifiers. Higher values indicate finer precision.

        Returns
        -------
        list
            A list of cell identifiers corresponding to the input latitude and longitude points.

        Raises
        ------
        ValueError
            If `cell_type` is not one of 'geohash', 's2', 's2_int', or 'h3'.

        Notes
        -----
        This function splits the input latitude and longitude lists into chunks and performs the cell
        conversion in parallel, with each chunk processed by a separate CPU core. This can significantly
        reduce processing time for large datasets.

        Examples
        --------
        >>> lats = [37.7749, 40.7128]
        >>> lons = [-122.4194, -74.0060]
        >>> cell_type = "h3"
        >>> res = 9
        >>> ppointcell(lats, lons, cell_type, res)
        ['8928308280fffff', '8a28308280fffff']
        """
        n_cores = cpu_count()

        # Prepare arguments for parallel processing
        lat_chunks = np.array_split(lats, 4 * n_cores)
        lon_chunks = np.array_split(lons, 4 * n_cores)
        args = zip(lat_chunks, lon_chunks, [cell_type] * 4 * n_cores, [res] * 4 * n_cores)

        # Parallelize the conversion using Pool.starmap
        with Pool(n_cores) as pool:
            cells = pool.starmap(SpatialIndex.pointcell, args)
        cells = [item for sublist in cells for item in sublist]  # Flatten the list of cells

        return cells

    @staticmethod
    def ppolycell(
        mdf: gpd.GeoDataFrame, cell_type: str, res: int, compact: bool = False, dump: str = None, verbose: bool = False
    ) -> Tuple[List[str], int]:
        """
        Performs a parallelised conversion of geometries in a GeoDataFrame to cell identifiers of a specified type
        (e.g., Geohash, S2, or H3), optionally compacting the result to reduce the number of cells.

        This function first divides the bounding box of the input GeoDataFrame into smaller grid cells, then calculates
        the intersection between these grid cells and the input geometries. The resulting geometries are processed in
        parallel to generate cell identifiers according to the specified `cell_type` and `res` (resolution). The result
        can be compacted to reduce the number of cells. Optionally, if `dump` is provided, the results are saved in multiple
        files, where the number of files is 4 times the number of CPU cores available in the system.

        Parameters
        ----------
        mdf : gpd.GeoDataFrame
            A GeoDataFrame containing geometries that need to be converted to cell identifiers.

        cell_type : str
            The type of cell identifier to use. Options are:
            - "geohash": Converts geometries to Geohash identifiers.
            - "s2": Converts geometries to S2 cell tokens.
            - "h3": Converts geometries to H3 cell tokens.

        res : int
            The resolution or precision level of the cell identifiers. Higher values indicate finer precision.

        compact : bool, optional, default=False
            If True, compact the resulting cells to reduce their number. This is typically applicable for S2 and H3 cells.

        dump : str, optional
            A string representing a valid directory path. If provided, the cells are saved in multiple files
            within the directory `/path/to/dir/cell_type/res/`. The number of output files will be 4 times the number
            of CPU cores available in the system. If not provided, the function returns the list of cell identifiers
            instead of saving them to files. Default is None.

        verbose : bool, optional, default=False
            If True, print timing and progress information to the console.

        Returns
        -------
        Tuple[List[str], int]
            - A list of cell identifiers as strings, corresponding to the geometries in the input GeoDataFrame.
            - The total number of unique cell identifiers.

        Raises
        ------
        ValueError
            If an invalid `cell_type` is provided. Supported types are "geohash", "s2", and "h3".

        Example
        -------
        >>> # Assuming `mdf` is a GeoDataFrame with geometries:
        >>> cells, count = ppolycell(mdf, cell_type="s2", res=10, compact=True, dump="~/Desktop/cells", verbose=True)
        >>> print(f"Generated {count} cells: {cells}")
        """
        if verbose:
            print(datetime.now())
            print("\nSlicing the bounding box of polygons ... ", end="")
            start_time = time()

        # Determine the number of slices and grid cells based on CPU cores
        n_cores = cpu_count()
        slices = 128 * n_cores

        # Calculate the bounding box dimensions
        minlon, minlat, maxlon, maxlat = mdf.total_bounds
        dlon = maxlon - minlon
        dlat = maxlat - minlat
        ratio = dlon / dlat

        # Calculate the number of grid cells in x and y directions
        x_cells = round(sqrt(slices) * ratio)
        y_cells = round(sqrt(slices) / ratio)

        # Calculate step size for grid cells
        steplon = dlon / x_cells
        steplat = dlat / y_cells

        # Create grid polygons based on bounding box slices
        grid_polygons = []
        for lat in np.arange(minlat, maxlat, steplat):
            for lon in np.arange(minlon, maxlon, steplon):
                llon, llat, ulon, ulat = (lon, lat, lon + steplon, lat + steplat)  # lower lat, upper lat
                polygon = Polygon([(llon, llat), (ulon, llat), (ulon, ulat), (llon, ulat)])
                grid_polygons.append(polygon)

        gmdf = gpd.GeoDataFrame(geometry=grid_polygons, crs=mdf.crs)  # Create a GeoDataFrame with grid polygons

        if verbose:
            elapsed_time = round(time() - start_time)
            print(f"{elapsed_time} seconds.   {slices} slices created.")

            print("Performing intersection between grid and polygons ... ", end="")
            start_time = time()

        # Perform intersection between input geometries and grid cells
        gmdf = gpd.overlay(mdf, gmdf, how="intersection")  # grid mdf

        if verbose:
            elapsed_time = round(time() - start_time)
            print(f"{elapsed_time} seconds.   {len(gmdf)} intersected slices.")

            print("Calculating cell IDs in parallel ... ", end="")
            start_time = time()

        # Shuffle geometries for even load distribution across chunks
        gmdf = gmdf.sample(frac=1)
        geom_chunks = np.array_split(list(gmdf.geometry), 4 * n_cores)
        args = zip(geom_chunks, [cell_type] * 4 * n_cores, [res] * 4 * n_cores, [dump] * 4 * n_cores)

        # Parallel processing to generate cells
        if dump:
            with Pool(n_cores) as pool:
                pool.starmap(SpatialIndex.polycell, args)
            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"{elapsed_time} seconds.")
            return
        else:
            with Pool(n_cores) as pool:
                cells = pool.starmap(SpatialIndex.polycell, args)
            cells = [item for sublist in cells for item in sublist]  # Flatten the list of cells

            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"{elapsed_time} seconds.")

            # Remove duplicates based on cell type
            if cell_type in {"geohash", "s2"}:
                if verbose:
                    print("Removing duplicate cells ... ", end="")
                    start_time = time()
                cells = list(set(cells))  # Remove duplicate cells
                if verbose:
                    elapsed_time = round(time() - start_time)
                    print(f"{elapsed_time} seconds.")

            cell_counts = len(cells)  # Total unique cell count

            # Compact the cells if needed
            if compact:
                if verbose:
                    print("Compacting cells ... ", end="")
                    start_time = time()
                cells = CellUtils.compact_cells(cells, cell_type)
                if verbose:
                    elapsed_time = round(time() - start_time)
                    print(f"{elapsed_time} seconds.")

            return cells, cell_counts

    @staticmethod
    def cellpoint(cells: List[Union[str, int]], cell_type: str) -> List[Tuple[float, float]]:
        """
        Converts a list of cell IDs into their corresponding centroids.

        This function supports various cell ID types: 'geohash', 'h3', 's2_int' (integer-based S2 cells),
        and 's2' (token-based S2 cells). For each cell in the list, it returns the latitude and longitude
        of the cell's center.

        Parameters
        ----------
        cells : list of str or int
            List of cell identifiers. The format of each cell ID depends on the specified `cell_type`:
            - For 'geohash': `cells` should be a list of geohash strings.
            - For 'h3': `cells` should be a list of H3 cell ID strings.
            - For 's2_int': `cells` should be a list of integer-based S2 cell IDs.
            - For 's2': `cells` should be a list of S2 token strings.
        cell_type : str
            Type of the cell ID format. Should be one of the following:
            - 'geohash': Geohash encoding.
            - 'h3': H3 hexagonal grid encoding.
            - 's2_int': Integer-based S2 cell ID.
            - 's2': S2 token (string-based cell ID).

        Returns
        -------
        list of tuple of float
            A list of tuples where each tuple contains the latitude and longitude (in degrees) of the
            center point for each cell ID.

        Raises
        ------
        ValueError
            If the `cell_type` is not one of 'geohash', 'h3', 's2', or 's2_int'.

        Examples
        --------
        >>> cellpoint(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        [(42.6, -5.6), (57.64911, 10.40744)]

        >>> cellpoint(["8928308280fffff"], cell_type="h3")
        [(37.775938728915946, -122.41795063018799)]

        >>> cellpoint([9744573459660040192], cell_type="s2_int")
        [(37.7749, -122.4194)]

        >>> cellpoint(["89c25c"], cell_type="s2")
        [(37.7749, -122.4194)]
        """
        if cell_type == "geohash":
            return [pygeohash.decode(cell) for cell in cells]
        elif cell_type == "h3":
            return [h3.cell_to_latlng(cell) for cell in cells]
        elif cell_type == "s2_int":
            return [(s2.CellId(cell).to_lat_lng().lat().degrees, s2.CellId(cell).to_lat_lng().lng().degrees) for cell in cells]
        elif cell_type == "s2":
            return [
                (s2.CellId.from_token(cell).to_lat_lng().lat().degrees, s2.CellId.from_token(cell).to_lat_lng().lng().degrees)
                for cell in cells
            ]
        else:
            raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', 's2_int', or 'h3'.")

    @staticmethod
    def cellpoly(cells: list, cell_type: str) -> tuple:
        """
        Converts a list of spatial cells to their corresponding geometries and resolution levels.

        The function takes a list of spatial cells (e.g., Geohash, H3, or S2) and converts each cell
        into a geometry object (Polygon) based on the specified cell type. It also calculates the resolution
        level for each cell.

        Parameters
        ----------
        cells : list
            A list of spatial cells represented as strings. Each cell corresponds to a spatial area
            in a specific grid system (e.g., Geohash, H3, or S2).

        cell_type : str
            The type of spatial cell system used. Accepted values are:
            - "geohash" : Geohash spatial indexing system.
            - "h3"      : H3 hexagonal spatial indexing system.
            - "s2"      : S2 spherical spatial indexing system.

        Returns
        -------
        tuple
            A tuple containing:
            - `res` : list of int
                A list of resolution levels corresponding to each cell in the input.
            - `geoms` : list of shapely.geometry.Polygon
                A list of Polygon geometries representing the spatial boundaries of the input cells.

        Raises
        ------
        ValueError
            If `cell_type` is not one of "geohash", "h3", or "s2".

        Example
        -------
        >>> from shapely.geometry import Polygon
        >>> cells = ["ezs42", "ezs43"]  # Geohash cells
        >>> cell_type = "geohash"
        >>> res, geoms = cellpoly(cells, cell_type)
        >>> print(res)
        [5, 5]  # Resolution levels of the input cells
        >>> print(geoms)
        [<shapely.geometry.polygon.Polygon object at 0x...>, <shapely.geometry.polygon.Polygon object at 0x...>]
        # Polygon geometries representing the spatial boundaries of the cells

        Notes
        -----
        The function supports three spatial indexing systems:
        - Geohash: Uses rectangular bounding boxes to represent cells.
        - H3: Uses hexagonal grid cells.
        - S2: Uses spherical grid cells.
        """
        # Check for valid cell_type
        if cell_type not in {"geohash", "h3", "s2"}:
            raise ValueError(f"Invalid cell_type '{cell_type}'. Accepted values are: 'geohash', 'h3', 's2'.")

        # Determine resolution level based on cell type
        res = [
            len(cell)
            if cell_type == "geohash"
            else int(cell[1], 16)
            if cell_type == "h3"
            else s2.CellId.from_token(cell).level()  # cell = token
            for cell in cells
        ]

        # Create geometry objects based on cell type
        geoms = [
            geohash_to_polygon(cell)
            if cell_type == "geohash"
            else Polygon(s2.s2_to_geo_boundary(cell, geo_json_conformant=True))
            if cell_type == "s2"
            # Shapely expects (lng, lat) format, so we reverse the coordinates returned by cell_to_boundary
            else Polygon([(lng, lat) for lat, lng in h3.cell_to_boundary(cell)])
            for cell in cells
        ]

        return res, geoms

    @staticmethod
    def pcellpoint(cells: List[Union[str, int]], cell_type: str) -> List[Tuple[float, float]]:
        """
        Converts a list of cell IDs into their corresponding latitude and longitude points in parallel.

        Parameters
        ----------
        cells : list of str or int
            List of cell identifiers.
        cell_type : str
            Type of the cell ID format.

        Returns
        -------
        list of tuple of float
            List of tuples containing the latitude and longitude (in degrees) of each cell ID.
        """
        n_cores = cpu_count()

        # Prepare arguments for parallel processing
        cell_chunks = np.array_split(cells, 4 * n_cores)
        args = zip(cell_chunks, [cell_type] * 4 * n_cores)

        # Parallelize the conversion using Pool.starmap
        with Pool(n_cores) as pool:
            points = pool.starmap(SpatialIndex.cellpoint, args)
        points = [item for sublist in points for item in sublist]  # Flatten the list of cells

        return points

    @staticmethod
    def pcellpoly(cells: List[Union[str, int]], cell_type: str) -> tuple:
        """
        Parallelized version of `cellpoly`, converting a list of spatial cells to geometries and resolution levels.

        Parameters
        ----------
        cells : list of str or int
            List of spatial cells in a specific grid system.

        cell_type : str
            Type of spatial cell system ("geohash", "h3", or "s2").

        Returns
        -------
        tuple
            A tuple containing:
            - `res` : list of int
                Resolution levels for each cell in the input.
            - `geoms` : list of shapely.geometry.Polygon
                Polygon geometries representing the boundaries of input cells.
        """
        n_cores = cpu_count()

        cell_chunks = np.array_split(cells, 4 * n_cores)
        # Convert each numpy array to a list which converts numpy.str_ to str
        cell_chunks = [arr.tolist() for arr in cell_chunks]
        args = zip(cell_chunks, [cell_type] * 4 * n_cores)

        with Pool(n_cores) as pool:
            results = pool.starmap(SpatialIndex.cellpoly, args)

        # Unpack `res` and `geoms` from the result tuples
        res = [r for result in results for r in result[0]]
        geoms = [g for result in results for g in result[1]]

        return res, geoms


class SpatialOps:
    """
    A utility class for performing advanced spatial operations on GeoDataFrames and geometries.

    This class provides methods for handling 3D geometries, converting LineStrings to Points,
    performing spatial intersections, and executing parallelized spatial overlay operations.
    It is designed to work with GeoPandas GeoDataFrames and Shapely geometries.

    Methods
    -------
    flatten_3d(geom: gpd.GeoSeries) -> List[Union[Polygon, MultiPolygon]]
        Converts 3D geometries in a GeoSeries to 2D geometries by removing the z-coordinate.

    line_to_points(row: gpd.GeoSeries) -> gpd.GeoDataFrame
        Splits a LineString geometry into individual Point geometries while preserving attributes.

    intersection(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, poly_id: Optional[str] = None) -> gpd.GeoDataFrame
        Performs a spatial intersection between two GeoDataFrames and returns the intersecting subset.

    quick_intersection(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, poly_id: Optional[str] = None) -> gpd.GeoDataFrame
        Performs an optimized spatial intersection using bounding box filtering and spatial indexing.

    poverlay(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, how: str = "intersection", keep_geom_type: bool = False) -> gpd.GeoDataFrame
        Executes a parallelized spatial overlay operation between two GeoDataFrames.

    Notes
    -----
    - The class relies on GeoPandas and Shapely for spatial operations and multiprocessing for parallelization.
    - Ensure that input GeoDataFrames have the same coordinate reference system (CRS) for accurate results.

    Examples
    --------
    >>> spatial_ops = SpatialOps()

    >>> # Example for flatten_3d
    >>> gdf_2d = spatial_ops.flatten_3d(gdf.geometry)
    >>> print(gdf_2d)

    >>> # Example for line_to_points
    >>> point_gdf = spatial_ops.line_to_points(line_gdf.iloc[0])
    >>> print(point_gdf)

    >>> # Example for intersection
    >>> result_gdf = spatial_ops.intersection(gdf1, gdf2, poly_id="region_id")
    >>> print(result_gdf)

    >>> # Example for quick_intersection
    >>> result_gdf = spatial_ops.quick_intersection(gdf1, gdf2, poly_id="region_id")
    >>> print(result_gdf)

    >>> # Example for poverlay
    >>> result_gdf = spatial_ops.poverlay(gdf1, gdf2, how="intersection")
    >>> print(result_gdf)
    """

    @staticmethod
    def intersection(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, poly_id: Optional[str] = None) -> gpd.GeoDataFrame:
        """
        Performs a spatial intersection between two GeoDataFrames and return the intersecting subset of the first GeoDataFrame.

        This function identifies geometries in `gdf1` that intersect with any geometries in `gdf2`. It adds a new column, `counts`,
        to `gdf2` representing the number of intersecting geometries for each feature in `gdf2`. If a `poly_id` column is specified,
        it also adds the geometry ID from `gdf2` to the intersected subset of `gdf1`.

        Parameters
        ----------
        gdf1 : geopandas.GeoDataFrame
            The first GeoDataFrame whose geometries are tested for intersection with `gdf2`.
        gdf2 : geopandas.GeoDataFrame
            The second GeoDataFrame containing geometries to intersect with `gdf1`.
        poly_id : str, optional
            The column name in `gdf2` containing unique geometry identifiers. If provided, the intersected subset of `gdf1`
            will include a new column `geom_id` indicating the geometry ID from `gdf2` that each feature intersects with.

        Returns
        -------
        geopandas.GeoDataFrame
            A new GeoDataFrame containing only the intersecting geometries from `gdf1` with respect to `gdf2`.
            If `poly_id` is provided, the intersected GeoDataFrame will also include a `geom_id` column.

        Examples
        --------
        >>> gdf1 = geopandas.read_file("data1.shp")
        >>> gdf2 = geopandas.read_file("data2.shp")
        >>> result_gdf = intersection(gdf1, gdf2, poly_id="region_id")

        Notes
        -----
        The function modifies `gdf2` in place by adding a `counts` column, which reflects the number of geometries
        in `gdf1` that intersect with each geometry in `gdf2`.

        """
        int_gdf = pd.DataFrame()  # Initialize an empty DataFrame to store intersecting geometries from gdf1
        counts = []  # List to store counts of intersecting geometries for each feature in gdf2

        for geom in gdf2.geometry:
            # Filter `gdf1` to retain only geometries that intersect with the current geometry in `gdf2`
            gdf = gdf1[gdf1.intersects(geom)]

            if poly_id is not None and len(gdf) > 0:
                # If `poly_id` is provided, retrieve the geometry ID from `gdf2` and assign it to `geom_id` column in `gdf`
                gid = gdf2[gdf2.geometry == geom][poly_id].iloc[0]
                gdf["geom_id"] = gid

            # Concatenate the intersecting geometries to the final DataFrame
            int_gdf = pd.concat([int_gdf, gdf])
            counts.append(len(gdf))  # Store the number of intersecting geometries for the current feature in gdf2

        gdf2["counts"] = counts  # Add the counts of intersecting geometries as a new column in gdf2
        return int_gdf  # Return the intersected subset of gdf1

    @staticmethod
    def quick_intersection(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, poly_id: Optional[str] = None) -> gpd.GeoDataFrame:
        """
        Performs a quick spatial intersection between two GeoDataFrames using bounding box optimization.

        This function identifies geometries in `gdf1` that intersect with any geometries in `gdf2`. It uses
        a spatial index to quickly filter `gdf1` geometries that are likely to intersect with the bounding
        box of each geometry in `gdf2`. It then performs a precise intersection check on this subset, improving
        the performance of the intersection operation.

        If a `poly_id` column is provided, the function adds a new `geom_id` column to the resulting intersected
        GeoDataFrame, storing the geometry ID from `gdf2` that each feature in `gdf1` intersects with. It also
        modifies `gdf2` by adding a `counts` column to indicate the number of intersecting geometries.

        Parameters
        ----------
        gdf1 : geopandas.GeoDataFrame
            The first GeoDataFrame whose geometries are tested for intersection with `gdf2`.
        gdf2 : geopandas.GeoDataFrame
            The second GeoDataFrame containing geometries to intersect with `gdf1`.
        poly_id : str, optional
            The column name in `gdf2` containing unique geometry identifiers. If provided, the intersected subset of `gdf1`
            will include a new column `geom_id` indicating the geometry ID from `gdf2` that each feature intersects with.

        Returns
        -------
        geopandas.GeoDataFrame
            A new GeoDataFrame containing only the intersecting geometries from `gdf1` with respect to `gdf2`.
            If `poly_id` is provided, the intersected GeoDataFrame will also include a `geom_id` column.

        Examples
        --------
        >>> gdf1 = geopandas.read_file("data1.shp")
        >>> gdf2 = geopandas.read_file("data2.shp")
        >>> result_gdf = quick_intersection(gdf1, gdf2, poly_id="region_id")

        Notes
        -----
        - This function modifies `gdf2` in place by adding a `counts` column, which reflects the number of geometries
          in `gdf1` that intersect with each geometry in `gdf2`.
        - It leverages spatial indexing using the `sindex` attribute of `gdf1` to quickly identify candidates for
          intersection, which significantly improves performance for large datasets.

        """
        int_gdf = pd.DataFrame()  # Initialize an empty DataFrame to store intersecting geometries from gdf1
        counts = []  # List to store counts of intersecting geometries for each feature in gdf2

        for geom in gdf2.geometry:
            # Get the indices of geometries in `gdf1` that are likely to intersect the bounding box of `geom` in `gdf2`
            pos_idx = list(gdf1.sindex.intersection(geom.bounds))

            # Select the subset of `gdf1` based on these indices
            pos_gdf = gdf1.iloc[pos_idx]

            # Filter the subset to retain only geometries that precisely intersect with `geom`
            pre_gdf = pos_gdf[pos_gdf.intersects(geom)]

            if poly_id is not None and len(pre_gdf) > 0:
                # If `poly_id` is provided, assign the geometry ID from `gdf2` to the `geom_id` column in `pre_gdf`
                gid = gdf2[gdf2.geometry == geom][poly_id].iloc[0]
                pre_gdf["geom_id"] = gid

            # Concatenate the precise intersecting geometries to the final intersected DataFrame
            int_gdf = pd.concat([int_gdf, pre_gdf])
            counts.append(len(pre_gdf))  # Store the number of intersecting geometries for the current feature in gdf2

        gdf2["counts"] = counts  # Add the counts of intersecting geometries as a new column in gdf2
        return int_gdf  # Return the intersected subset of gdf1

    @staticmethod
    def poverlay(
        gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, how: str = "intersection", keep_geom_type: bool = False
    ) -> gpd.GeoDataFrame:
        """
        Performs a spatial overlay operation between two GeoDataFrames in parallel using multiple CPU cores.

        This function divides the first GeoDataFrame into chunks according to the number of available CPU cores
        and applies the specified overlay operation (e.g., intersection, union, difference) in parallel on each chunk
        with respect to the second GeoDataFrame. The results are then concatenated and returned as a single GeoDataFrame.

        Parameters
        ----------
        gdf1 : gpd.GeoDataFrame
            The first GeoDataFrame to be used in the spatial overlay operation.
        gdf2 : gpd.GeoDataFrame
            The second GeoDataFrame to be used in the spatial overlay operation.
        how : str, optional
            The type of overlay operation to perform. Options include "intersection", "union", "difference",
            "symmetric_difference", and "identity". Defaults to "intersection".
        keep_geom_type : bool, optional
            Whether to retain the original geometry type (e.g., Polygon, LineString) in the resulting overlay.
            If set to True, only features of the same geometry type are retained. Defaults to False.

        Returns
        -------
        gpd.GeoDataFrame
            A new GeoDataFrame resulting from the spatial overlay operation, with the same coordinate reference system
            (CRS) as the first input GeoDataFrame (`gdf1`).

        Examples
        --------
        >>> gdf1 = gpd.GeoDataFrame({"geometry": [Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])]})
        >>> gdf2 = gpd.GeoDataFrame({"geometry": [Polygon([(1, 1), (3, 1), (3, 3), (1, 3)])]})
        >>> result_gdf = poverlay(gdf1, gdf2, how="intersection")
        >>> print(result_gdf)
                                                     geometry
        0  POLYGON ((2.00000 1.00000, 2.00000 2.00000, 1....

        Notes
        -----
        - The spatial overlay operation is performed using the `geopandas.overlay` function. The parallelization is achieved
          using the `multiprocessing` library to divide and distribute the overlay operations across multiple CPU cores.
        - Ensure that both GeoDataFrames (`gdf1` and `gdf2`) have the same coordinate reference system (CRS) before applying
          the overlay operation to avoid unexpected results.

        Raises
        ------
        ValueError
            If the `how` parameter is not one of the supported overlay operation types: "intersection", "union",
            "difference", "symmetric_difference", or "identity".
        """
        # Determine the number of CPU cores available for parallel processing
        n_cores = cpu_count()

        # Split the first GeoDataFrame into chunks for parallel processing
        gdf1_chunks = np.array_split(gdf1, n_cores)

        # Create a list of the second GeoDataFrame repeated for each chunk
        gdf2_chunks = [gdf2] * n_cores

        # Prepare inputs for the parallel processing pool
        inputs = zip(gdf1_chunks, gdf2_chunks, [how] * n_cores, [keep_geom_type] * n_cores)

        # Create a multiprocessing pool and apply the overlay function in parallel on each chunk
        with Pool(n_cores) as pool:
            gdf = pd.concat(pool.starmap(gpd.overlay, inputs), ignore_index=True)

        return gdf

    @staticmethod
    def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        """
        Calculates the great-circle distance between two points on the Earth's surface.

        The haversine formula determines the shortest distance over the Earth's surface
        between two points given their latitudes and longitudes. The result is the
        distance in meters, based on a mean Earth radius.

        Parameters
        ----------
        lat1 : float
            Latitude of the first point in decimal degrees.
        lon1 : float
            Longitude of the first point in decimal degrees.
        lat2 : float
            Latitude of the second point in decimal degrees.
        lon2 : float
            Longitude of the second point in decimal degrees.

        Returns
        -------
        float
            The great-circle distance between the two points in meters.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Haversine_formula
        .. [2] https://en.wikipedia.org/wiki/Longitude

        Examples
        --------
        >>> haversine(52.2296756, 21.0122287, 41.8919300, 12.5113300)
        1319743.483

        Notes
        -----
        The mean Earth radius is taken as 6,371,008.8 meters.
        a = 6378137.0        # Equatorial radius
        b = 6356752.3142     # Polar radius
        R = (2*a + b)/3      # Mean radius = 6371008.7714
        """
        r = 6371008.8  # Mean Earth radius in meters
        lat1, lat2, dlon = radians(lat1), radians(lat2), radians(lon2 - lon1)
        dlat = lat2 - lat1

        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))  # Angular distance in radians
        return r * c  # Distance in meters

    @staticmethod
    def vincenty(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
        """
        Calculates the geodesic distance between two points on the Earth's surface
        using the Vincenty formula, which accounts for the Earth's ellipsoidal shape.

        Parameters:
        - lat1, lon1: Latitude and longitude of the first point (in degrees).
        - lat2, lon2: Latitude and longitude of the second point (in degrees).

        Returns:
        - Distance between the two points in meters.

        Notes:
        - This implementation may encounter numerical issues, such as divide-by-zero errors,
          in edge cases where the points are on opposite sides of the Earth or on the same meridian
          e.g., from (0,0) to (0,90).However, for points (0,0) to (0.001,90), the distance calculation
          is accurate within a small error margin (about 9.3e-06 meters).

        - The error in the above approach can be significant for very small distances,
          such as between (0,0) and (0,0.001).
        """

        # Constants for WGS-84 ellipsoid
        a = 6378137.0  # Equatorial radius in meters
        f = 1 / 298.257223563  # Flattening
        b = a * (1 - f)  # Polar radius

        # Convert degrees to radians
        lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

        # Differences in longitude
        ll = lon2 - lon1

        # Iterative Vincenty formula
        u1 = math.atan((1 - f) * math.tan(lat1))
        u2 = math.atan((1 - f) * math.tan(lat2))
        sin_u1 = math.sin(u1)
        cos_u1 = math.cos(u1)
        sin_u2 = math.sin(u2)
        cos_u2 = math.cos(u2)

        lambda_ = ll
        lambda_prev = 0
        max_iterations = 1000
        tolerance = 1e-12

        for _ in range(max_iterations):
            sin_lambda = math.sin(lambda_)
            cos_lambda = math.cos(lambda_)
            sin_sigma = math.sqrt((cos_u2 * sin_lambda) ** 2 + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) ** 2)
            cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
            sigma = math.atan2(sin_sigma, cos_sigma)
            sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma
            cos2_alpha = 1 - sin_alpha**2
            cos2_sigma_m = cos_sigma - 2 * sin_u1 * sin_u2 / cos2_alpha
            cc = f / 16 * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))
            lambda_prev = lambda_
            lambda_ = ll + (1 - cc) * f * sin_alpha * (
                sigma + cc * sin_sigma * (cos2_sigma_m + cc * cos_sigma * (-1 + 2 * cos2_sigma_m**2))
            )

            if abs(lambda_ - lambda_prev) < tolerance:
                break
        else:
            raise ValueError("Vincenty formula did not converge")

        u2 = cos2_alpha * (a**2 - b**2) / (b**2)
        aa = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        bb = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        delta_sigma = (
            bb
            * sin_sigma
            * (
                cos2_sigma_m
                + bb
                / 4
                * (
                    cos_sigma * (-1 + 2 * cos2_sigma_m**2)
                    - bb / 6 * cos2_sigma_m * (-3 + 4 * sin_sigma**2) * (-3 + 4 * cos2_sigma_m**2)
                )
            )
        )
        s = b * aa * (sigma - delta_sigma)

        return s

    @staticmethod
    def geocoding_google(address_or_zipcode: str, api_key: str) -> pd.Series:
        """
        Returns geographic coordinates (latitude and longitude) for a given address or zip code using the Google Geocoding API.

        This function utilizes the Google Geocoding API to convert a given address or zip code into geographic coordinates.
        The function returns the latitude and longitude as a pandas Series. If the request is unsuccessful or the address
        is not found, the function returns a Series with `(None, None)`.

        Parameters
        ----------
        address_or_zipcode : str
            A text-based address or zip code that needs to be geocoded.
        api_key : str
            A valid Google Maps API key required to access the Google Geocoding service.

        Returns
        -------
        pd.Series
            A pandas Series containing the latitude and longitude as floats. If the request fails or the address is not found,
            returns a Series with `(None, None)`.

        Examples
        --------
        >>> df[["lat", "lon"]] = df.apply(lambda row: geocoding_google(row.address, "your_api_key"), axis=1)
        >>> result = geocoding_google("1600 Amphitheatre Parkway, Mountain View, CA", "your_api_key")
        >>> print(result)
        lat    37.4224764
        lon   -122.0842499
        dtype: float64

        Notes
        -----
        - Make sure to enable the Google Geocoding API in your Google Cloud Console and provide a valid API key.
        - The API might return ambiguous results if the input address is incomplete or vague.
        - Consider handling `None` values in the returned Series if the API fails to find the address or the request limit is exceeded.

        Raises
        ------
        Exception
            If there is an error in the API request or response parsing, an exception is raised with an error message.
        """
        lat, lon = None, None
        base_url = "https://maps.googleapis.com/maps/api/geocode/json"
        endpoint = f"{base_url}?address={address_or_zipcode}&key={api_key}"
        r = requests.get(endpoint)
        if r.status_code not in range(200, 299):
            return None, None
        try:
            """
            This try block incase any of our inputs are invalid. This is done instead
            of actually writing out handlers for all kinds of responses.
            """
            results = r.json()["results"][0]
            lat = results["geometry"]["location"]["lat"]
            lon = results["geometry"]["location"]["lng"]
        except Exception:
            pass  # Handle any errors that may occur
        return pd.Series([lat, lon])

    @staticmethod
    def reverse_geocoding_google(lat: float, lon: float, api_key: str) -> str:
        """
        Returns the postal code for a given geographic coordinate (latitude, longitude) using the Google Geocoding API.

        This function makes a reverse geocoding request to the Google Geocoding API to obtain the postal code associated
        with the provided latitude and longitude. If the postal code is found, it is returned as a string. If not,
        `None` is returned.

        Parameters
        ----------
        lat : float
            The latitude of the location to reverse geocode.
        lon : float
            The longitude of the location to reverse geocode.
        api_key : str
            A valid Google Maps API key for accessing the geocoding service.

        Returns
        -------
        str
            The postal code corresponding to the input geographic coordinates, if found. Returns `None` if no postal code
            is found or if the request fails.

        Examples
        --------
        >>> reverse_geocoding_google(37.4224764, -122.0842499, "your_api_key")
        '94043'

        >>> df["postcode"] = df.apply(lambda row: reverse_geocoding_google(row.lat, row.lon, "your_api_key"), axis=1)
        """
        lat = 0 if abs(lat) < 0.0001 else lat  # Prevent invalid 'latlng' error for very small values.
        lon = 0 if abs(lon) < 0.0001 else lon

        # Make the reverse geocoding request
        url = f"https://maps.googleapis.com/maps/api/geocode/json?latlng={lat},{lon}&key={api_key}"
        response = requests.get(url)
        data = response.json()

        # Parse the response to extract the postal code
        if "results" in data and len(data["results"]) > 0:
            for item in data["results"]:
                for component in item.get("address_components", []):
                    if "postal_code" in component.get("types", []) and len(component["types"]) == 1:
                        return component.get("long_name", None)
        return None
