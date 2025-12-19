import collections
import json
import math
import os
import re
from collections import deque
from datetime import datetime
from math import atan2, cos, radians, sin, sqrt
from multiprocessing import Pool, cpu_count
from time import time

import folium  # Folium is a Python library used for visualizing geospatial data. Actually, it's a Python wrapper for Leaflet which is a leading open-source JavaScript library for plotting interactive maps.
import geopandas as gpd
import h3
import lonboard as lb
import matplotlib
import matplotlib.colors
import numpy as np
import pandas as pd
import pygeohash
import pyproj
import requests
import shapely
from branca.element import MacroElement, Template
from folium import plugins
from lonboard.basemap import CartoStyle
from s2 import s2
from scipy.spatial import KDTree
from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform, unary_union
from shapely.prepared import prep


class Karta:
    @staticmethod
    def _base_map(sw: list, ne: list) -> folium.Map:
        """
        Creates a base map with multiple tile layers and fits the map to the specified bounding box.

        This function initializes a Folium map object with multiple tile layers, including:
        - `Bright Mode` (CartoDB Positron)
        - `Dark Mode` (CartoDB Dark Matter)
        - `Satellite` (Esri World Imagery)
        - `OpenStreetMap` (OSM)

        It then fits the map's view to the bounding box defined by the southwest (`sw`) and northeast (`ne`) coordinates.

        Parameters
        ----------
        sw : list
            The southwest coordinate [latitude, longitude] of the bounding box to fit the map view.

        ne : list
            The northeast coordinate [latitude, longitude] of the bounding box to fit the map view.

        Returns
        -------
        folium.Map
            A Folium map object with multiple tile layers and the view fitted to the provided bounding box.

        Examples
        --------
        >>> sw = [51.2652, -0.5426]  # Southwest coordinate (London, UK)
        >>> ne = [51.7225, 0.2824]  # Northeast coordinate (London, UK)
        >>> karta = Karta._base_map(sw, ne)
        >>> karta.save("map.html")  # Save the map to an HTML file
        """
        # Initialize the base map without any default tiles
        karta = folium.Map(tiles=None)

        # Add OpenStreetMap (OSM) tile layer
        folium.TileLayer("openstreetmap", name="OSM", max_zoom=19).add_to(karta)

        # Add a satellite tile layer (Esri World Imagery)
        folium.TileLayer(
            name="Satellite",
            attr='© <a href="https://www.esri.com/en-us/legal/overview">Esri</a>',
            tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
            overlay=False,
            control=True,
            max_zoom=19,
        ).add_to(karta)

        # Add OpenTopoMap as a tile layer
        folium.TileLayer(
            tiles="https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
            attr='© <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors © <a href="https://opentopomap.org/">OpenTopoMap</a>',
            name="Outdoors",
        ).add_to(karta)

        # Dictionary of tile layers to be added
        tiles = {
            "cartodbdark_matter": "Dark",
            "cartodbpositron": "Light",
        }

        # Add each tile layer to the map
        for item in tiles:
            folium.TileLayer(item, name=tiles[item], max_zoom=21).add_to(karta)

        # Fit the map's view to the bounding box defined by the southwest and northeast coordinates
        karta.fit_bounds([sw, ne])

        return karta

    @staticmethod
    def _select_color(col: int | float | str, head: int | None = None, tail: int | None = None) -> str:
        """
        Generates a consistent color based on the input column value by mapping it to a predefined color palette.

        This function uses a set color palette and maps the given column value to a color. If the column value is a string,
        a substring can be selected using `head` and `tail` indices, and it will be converted to a numerical index. If the
        column value is an integer, it will directly be mapped to a color using modulo arithmetic.

        Parameters
        ----------
        col : int, float or str
            The column value to be mapped to a color. It can be either a number or a string.
            - If a number, it is directly used for color mapping.
            - If a string, it will be cleaned of non-alphanumeric characters, and a substring defined by `head` and `tail`
              can be selected for mapping.

        head : int, optional
            The starting index of the substring to be used for color mapping if `col` is a string. Default is None.

        tail : int, optional
            The ending index of the substring to be used for color mapping if `col` is a string. Default is None.

        Returns
        -------
        str
            A hexadecimal color code selected from the predefined palette corresponding to the input column value.


        Examples
        --------
        >>> Karta._select_color("Category1")
        '#e6194b'  # Red color from the palette

        >>> Karta._select_color(5)
        '#3cb44b'  # Green color from the palette

        >>> Karta._select_color("Example", head=0, tail=3)
        '#e12348'  # Bright Red from the palette
        """
        # Predefined color palette
        palette = [
            "#e6194b",  # red
            "#4363d8",  # blue
            "#3cb44b",  # green
            "#800000",  # maroon (dark red)
            "#008080",  # teal (dark green)
            "#000080",  # navy (dark blue)
            "#f58231",  # orange
            "#911eb4",  # purple
            "#808000",  # olive
            "#9a6324",  # brown
            "#f032e6",  # magenta
            "#dfb119",  # dark yellow
            "#42d4f4",  # cyan
            "#808080",  # grey
            "#e12348",  # Bright Red
            "#dc2c46",  # Strong Red
            "#d73644",  # Vivid Red
            "#cd4a40",  # Deep Red
            "#c8543e",  # Intense Red
            "#c25e3c",  # Fire Red
            "#bd683a",  # Scarlet
            "#b77238",  # Fiery Orange
            "#b27c36",  # Tangerine
            "#ad8634",  # Burnt Orange
        ]

        # Handle NaN (return black)
        if col is None or (isinstance(col, float) and math.isnan(col)):
            return "#000000"

        if isinstance(col, (int, float)):  # Check for both int and float
            idx = int(col) % len(palette)  # Get color index using modulo arithmetic
        else:
            col = str(col)  # Convert to string
            col = re.sub(r"[\W_]+", "", col)  # Remove non-alphanumeric characters
            idx = int(col[head:tail], 36) % len(palette)  # Convert substring to a number base 36 (36 = 10 digits + 26 letters)

        return palette[idx]

    @staticmethod
    def _add_point(
        row: pd.Series,
        karta: folium.Map,
        color: str = "black",
        speed_field: str = "speed",
        speed_limit_field: str = "speedlimit",
        opacity: float = 0.5,
        radius: int = 3,
        weight: int = 6,
        popup_dict: dict = None,
        x: str = None,
        y: str = None,
    ) -> None:
        """
        Adds a point (marker) to a Folium map based on the specified parameters and data in the provided row.

        The function attempts to extract coordinates from a geometry column if available, or directly from `x` and `y` columns
        (longitude and latitude). It then adds a circle marker to the Folium map (`karta`) using the specified color, radius,
        and other style parameters.

        Parameters
        ----------
        row : pd.Series
            A row of data containing either a 'geometry' attribute or x/y columns for coordinates.

        karta : folium.Map
            A Folium map object to which the marker will be added.

        color : str
            Specifies the color of the marker. If "speed" is passed, the marker color is determined by comparing
            the 'speed' and 'speedlimit' values in the row (e.g., blue for under the speed limit, black for very high speeds).
            Otherwise, it can be a column name in the `row` to create a unique color from that column's value.

        opacity : float, optional
            Opacity of the marker (default is 0.5).

        radius : int, optional
            Radius of the circle marker (default is 3).

        weight : int, optional
            Weight (thickness) of the circle marker's border (default is 6).

        popup_dict : dict, optional
            A dictionary where keys are labels and values are column names in the row. This dictionary is used to create
            an HTML popup with the specified labels and values (default is None).

        x : str, optional
            Column name for longitude, if 'geometry' attribute is not present (default is None).

        y : str, optional
            Column name for latitude, if 'geometry' attribute is not present (default is None).

        Returns
        -------
        None
            The function modifies the Folium map in place and does not return anything.

        Examples
        --------
        >>> row = pd.Series({"geometry": Point(40.748817, -73.985428), "color": "red"})
        >>> karta = folium.Map(location=[40.748817, -73.985428], zoom_start=12)
        >>> Karta._add_point(row, karta, "color")
        """
        try:
            # Attempt to extract coordinates from the geometry column if present
            location = [row.geometry.y, row.geometry.x]
        except Exception:
            # If geometry is not present, use x and y columns for location
            location = [row[y], row[x]]  # x, y: lon, lat column names in DataFrame

        if any(math.isnan(item) for item in location):
            return None

        if color == speed_field:
            if pd.isna(row[speed_limit_field]) or row[speed_limit_field] <= 0:
                color = "purple"
            elif row[speed_field] <= row[speed_limit_field]:
                color = "blue"
            elif row[speed_field] < 1.1 * row[speed_limit_field]:
                color = "green"
            elif row[speed_field] < 1.2 * row[speed_limit_field]:
                color = "yellow"
            elif row[speed_field] < 1.3 * row[speed_limit_field]:
                color = "orange"
            elif row[speed_field] < 1.4 * row[speed_limit_field]:
                color = "red"
            else:
                color = "black"
        # Determine color if column is specified
        elif color in row.index:  # color in DataFrame columns
            color = Karta._select_color(row[color])

        # Create a popup HTML if popup_dict is provided
        popup = "".join(f"{item}: <b>{row[popup_dict[item]]}</b><br>" for item in popup_dict) if popup_dict else None

        # Add a CircleMarker to the map with the specified parameters
        folium.CircleMarker(
            location=location, radius=radius, color=color, opacity=opacity, weight=weight, tooltip=popup
        ).add_to(karta)

    @staticmethod
    def _add_line(
        row: pd.Series,
        karta: folium.Map,
        color: str = "blue",
        opacity: float = 0.5,
        weight: int = 6,
        popup_dict: dict = None,
    ) -> None:
        """
        Adds a polyline (line) to a Folium map based on geometry data from a Pandas Series row.

        Parameters
        ----------
        row : pd.Series
            A Pandas Series containing a 'geometry' column with a LineString object.
        karta : folium.Map
            The Folium map to which the polyline will be added.
        color : str, optional
            The color of the polyline. Defaults to "blue".
        opacity : float, optional
            The opacity of the polyline (range: 0.0 to 1.0). Defaults to 0.5.
        weight : int, optional
            The thickness of the polyline in pixels. Defaults to 6.
        popup_dict : dict, optional
            A dictionary mapping column names to popup labels.

        Returns
        -------
        None
            The function modifies the Folium map in place and does not return anything.

        Example
        -------
        >>> row = pd.Series({"geometry": LineString([(-74.006, 40.7128), (-73.9352, 40.7306)]), "road_name": "Broadway"})
        >>> karta = folium.Map(location=[40.72, -74.00], zoom_start=12)
        >>> Karta._add_line(row, karta, color="red", popup_dict={"Name": "road_name"})
        """
        coordinates = [(coord[1], coord[0]) for coord in row.geometry.coords]

        # Handle color selection
        color = Karta._select_color(row[color]) if color in row.index else color

        # Create popup if popup_dict is provided
        popup = "".join(f"{item}: <b>{row[popup_dict[item]]}</b><br>" for item in popup_dict) if popup_dict else None

        # Add polyline to the map
        folium.PolyLine(coordinates, color=color, weight=weight, opacity=opacity, tooltip=popup).add_to(karta)

    @staticmethod
    def _add_poly(
        row: pd.Series,
        karta: folium.Map,
        fill_color: str = "red",
        highlight_color: str = "green",
        fill_opacity: float = 0.25,
        highlight_opacity: float = 0.5,
        line_width: float = 0.3,
        popup_dict: dict = None,
    ) -> None:
        """
        Adds a polygon to a Folium map based on the specified parameters and data in the provided row.

        This function creates a polygon (GeoJson) object for the specified row's geometry and adds it to the Folium map (`karta`).
        It allows customization of fill color, line width, and popups. The function also defines style and highlight properties
        for the polygon.

        Parameters
        ----------
        row : pd.Series
            A row of data containing a 'geometry' attribute that defines the polygon shape.

        karta : folium.Map
            A Folium map object to which the polygon will be added.

        fill_color : str
            Column name to determine the fill color of the polygon. If the column is present in the row, the color is extracted
            using the `_select_color` function.

        line_width : int
            The width of the border (outline) of the polygon.

        popup_dict : dict, optional
            A dictionary where keys are labels and values are column names in the row. This dictionary is used to create an
            HTML popup with the specified labels and values for the polygon (default is None).

        Returns
        -------
        None
            The function modifies the Folium map in place and does not return anything.

        Examples
        --------
        >>> row = pd.Series({"geometry": Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), "fill_color": "blue"})
        >>> karta = folium.Map(location=[0.5, 0.5], zoom_start=10)
        >>> Karta._add_poly(row, karta, "fill_color", line_width=2)
        """
        # Determine fill color if specified column is present
        fill_color = Karta._select_color(row[fill_color]) if fill_color in row.index else fill_color

        # Style function to apply to the polygon
        def style_function(x):
            return {
                "fillColor": fill_color,
                "color": "black",  # Border color
                "fillOpacity": fill_opacity,
                "weight": line_width,
            }

        # Highlight style function when the polygon is hovered over
        def highlight_function(x):
            return {
                "fillColor": highlight_color,  # fill_color,
                "color": "black",  # Border color
                "fillOpacity": highlight_opacity,
                "weight": line_width,
            }

        # Create a popup if a popup dictionary is provided
        popup = "".join(f"{item}: <b>{row[popup_dict[item]]}</b><br>" for item in popup_dict) if popup_dict else None

        # Create a GeoJson object from the row's geometry and add it to the map
        gjson = row.geometry.__geo_interface__
        gjson = folium.GeoJson(data=gjson, style_function=style_function, highlight_function=highlight_function, tooltip=popup)
        gjson.add_to(karta)

    @staticmethod
    def plp(
        gdf_list: pd.DataFrame | gpd.GeoDataFrame | list[pd.DataFrame | gpd.GeoDataFrame] | None = None,
        # Point
        cluster: bool = False,
        heatmap: bool = False,
        heatmap_radius: int = 12,
        line: bool = False,
        antpath: bool = False,
        point_color: str = "blue",
        speed_field: str = "speed",
        speed_limit_field: str = "speedlimit",
        point_opacity: float = 0.5,
        point_radius: int = 3,
        point_weight: int = 6,
        point_popup: dict | None = None,
        buffer_radius: int = 0,
        ring_inner_radius: int = 0,
        ring_outer_radius: int = 0,
        x: str | None = None,
        y: str | None = None,
        # LineString
        line_color: str = "blue",
        line_opacity: float = 0.5,
        line_weight: int = 6,
        line_popup: dict | None = None,
        # Polygon
        centroid: bool = False,  # if True it shows centroids of polygons on the map.
        fill_color: str = "red",
        highlight_color: str = "green",
        fill_opacity: float = 0.25,
        highlight_opacity: float = 0.5,
        line_width: float = 0.3,
        poly_popup: dict | None = None,
        geohash_res: int = 0,
        s2_res: int = -1,
        h3_res: int = -1,
        force_full_cover: bool = True,
        geohash_inner: bool = False,
        compact: bool = False,
        # Cells and OSM objects
        cells: list[str] | None = None,
        cell_type: str | None = None,  # list of geohash, S2 or H3 cell IDs
        osm_ways: list[int] | None = None,  # list of OSM way IDs (lines or polygons) and Overpass API URL to query from
        url: str | None = "https://overpass-api.de/api/interpreter",  # OpenStreetMap server URL
    ) -> folium.Map:
        """
        plp (points, lines, polygons) creates a Folium map with points, lines, or polygons based on the input geospatial data.
        The function `plp` allows users to add different geometrical elements (points, lines, polygons) to a Folium map.
        It supports various visual styles and configurations, such as clustering, heatmaps, and geohash or cell-based layers.

        Parameters
        ----------
        gdf_list : list of gpd.GeoDataFrame or pd.DataFrame, optional
            List of GeoDataFrames or DataFrames containing geometrical data to be plotted. If a single DataFrame is provided,
            it will be wrapped in a list internally.

        cluster : bool, default False
            If True, clusters points together based on their proximity using Folium's `MarkerCluster`.

        heatmap : bool, default False
            If True, creates a heatmap layer using Folium's `HeatMap` for points.

        heatmap_radius : int, default 12
            Radius of each point of the heatmap

        line : bool, default False
            If True, connects points using Folium's `PolyLine` to form lines.

        antpath : bool, default False
            If True, creates animated ant paths for the line geometries using Folium's `AntPath`.

        point_color : str, default "blue"
            Color of the points when displayed on the map.

        speed_field : str, default "speed"
            Name of the speed field in DataFrame or GeoDataFrame.

        speed_limit_field : str, default "speedlimit"
            Name of the speed limit field in DataFrame or GeoDataFrame.

        point_opacity : float, default 0.5
            Opacity of the points. Value should be between 0 and 1.

        point_radius : int, default 3
            Radius of the points in pixels.

        point_weight : int, default 6
            Weight (thickness) of the point outline. Typically set to twice the `point_radius`.

        point_popup : dict, optional
            Dictionary where keys are labels and values are column names in the DataFrame. Used to create HTML popups with
            attributes of each point.

        buffer_radius : float, default 0
            Buffer radius (in meters) to create a buffer around each point. Set to 0 to disable buffering.

        ring_inner_radius : float, default 0
            Inner radius of ring buffers around points. Only used if `ring_outer_radius` is set.

        ring_outer_radius : float, default 0
            Outer radius of ring buffers around points. If set, creates a ring around each point.

        x : str, optional
            Column name for the x-coordinate (longitude). Specify it to use the column other than that containing 'lon'.

        y : str, optional
            Column name for the y-coordinate (latitude). Specify it to use the column other than that containing 'lat'.

        line_color : str, default "blue"
            Color of the lines connecting points or LineString geometries.

        line_opacity : float, default 0.5
            Opacity of the lines. Value should be between 0 and 1.

        line_weight : int, default 6
            Thickness of the lines.

        line_popup : dict, optional
            Dictionary where keys are labels and values are column names in the DataFrame. Used to create HTML popups with
            attributes of each line.

        centroid : bool, default False
            If True, displays the centroids of polygon geometries on the map.

        fill_color : str, default "red"
            Fill color for polygon geometries.

        highlight_color : str, default "green"
            Color used to highlight polygons when hovered.

        line_width : float, default 0.3
            Thickness of polygon outlines.

        poly_popup : dict, optional
            Dictionary where keys are labels and values are column names in the DataFrame. Used to create HTML popups with
            attributes of each polygon.

        geohash_res : int, default 0
            Resolution for creating geohash-based polygonal layers. Set to 0 to disable.

        s2_res : int, default -1
            Resolution for creating S2-based polygonal layers. Set to -1 to disable.

        h3_res : int, default -1
            Resolution for creating H3-based polygonal layers. Set to -1 to disable.

        geohash_inner : bool, default False
            If True, shows only inner geohash cells. Does not work if `compact` is set to True.

        compact : bool, default False
            If True, creates compact representation of geohash, S2, or H3 cells.

        cells : list, optional
            List of geohash, S2, or H3 cell IDs to visualize.

        cell_type : str, optional
            Type of cells used in `cells` parameter. Can be 'geohash', 's2', or 'h3'.

        osm_ways : list of int, optional
            List of OSM way IDs to visualize as lines or polygons.

        url : str, optional
            Overpass API URL to query OSM geometries by `osm_ways` parameter.

        Returns
        -------
        folium.Map
            A Folium map object with the added geometrical features based on input parameters.

        Examples
        --------
        >>> # Example usage
        >>> gdf = gpd.read_file("path/to/shapefile.shp")
        >>> plp(gdf)
        """

        # Handle `cells` input by converting cell IDs to geometries
        if cells:
            # cell_poly is faster than the parallelized pcell_poly because, for small datasets, the overhead of parallelization increases the total runtime.
            geoms, res = SpatialIndex.cell_poly(cells, cell_type=cell_type)
            gdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")
            karta = Karta.plp(gdf, poly_popup={"ID": "id", "Resolution": "res"})
            return karta

        # Handle `osm_ways` input by converting OSM way IDs to geometries
        if osm_ways:
            geoms = OSMUtils.ways_to_geom(osm_ways, url)
            gdf = gpd.GeoDataFrame({"way_id": osm_ways, "geometry": geoms}, crs="EPSG:4326")
            if isinstance(gdf.geometry[0], LineString):
                karta = Karta.plp(gdf, line_popup={"way_id": "way_id"}, line_color="red")
            else:
                karta = Karta.plp(
                    gdf,
                    geohash_res=geohash_res,
                    s2_res=s2_res,
                    h3_res=h3_res,
                    geohash_inner=geohash_inner,
                    compact=compact,
                    poly_popup={"way_id": "way_id"},
                    fill_color="red",
                )
            return karta

        # Ensure `gdf_list` is always a list of GeoDataFrames or DataFrames
        if isinstance(gdf_list, pd.DataFrame):
            gdf_list = [gdf_list]

        # Iterate through the list of GeoDataFrames to update bounding box
        for gdf in gdf_list:
            if not isinstance(gdf, gpd.GeoDataFrame):  # if pd.DataFrame
                if not x:  # if x is not specified, determine longitude and latitude columns
                    x = [col for col in gdf.columns if "lon" in col.lower() or "lng" in col.lower()][0]
                    y = [col for col in gdf.columns if "lat" in col.lower()][0]
                lons = gdf[x]
                lats = gdf[y]
                minlat, minlon, maxlat, maxlon = min(lats), min(lons), max(lats), max(lons)
            else:  # If input is a GeoDataFrame, use total_bounds to get the bounding box
                minlon, minlat, maxlon, maxlat = gdf.total_bounds

        # Create a base map using the bounding box
        sw = [minlat, minlon]  # South West (bottom left corner)
        ne = [maxlat, maxlon]  # North East (top right corner)
        karta = Karta._base_map(sw, ne)  # Initialize folium map with the bounding box

        # Iterate through each DataFrame or GeoDataFrame in the list to add layers to the map
        for i, gdf in enumerate(gdf_list, start=1):
            geom = gdf.geometry.values[0] if isinstance(gdf, gpd.GeoDataFrame) else None
            # Handle Polygon geometries
            if isinstance(geom, Polygon) or isinstance(geom, MultiPolygon):
                group_polygon = folium.FeatureGroup(name=f"{i}- Polygon")
                gdf.apply(
                    Karta._add_poly,
                    karta=group_polygon,
                    fill_color=fill_color,
                    highlight_color=highlight_color,
                    fill_opacity=fill_opacity,
                    highlight_opacity=highlight_opacity,
                    line_width=line_width,
                    popup_dict=poly_popup,
                    axis=1,
                )
                group_polygon.add_to(karta)

                if centroid:  # Show centroids of polygons if `centroid=True`
                    group_centroid = folium.FeatureGroup(name=f"{i}- Centroid")
                    cdf = gpd.GeoDataFrame({"geometry": gdf.centroid}, crs="EPSG:4326")  # centroid df
                    cdf.apply(Karta._add_point, karta=group_centroid, axis=1)
                    group_centroid.add_to(karta)
            # Handle LineString geometries
            elif isinstance(geom, LineString):
                group_line = folium.FeatureGroup(name=f"{i}- Line")
                gdf.apply(
                    Karta._add_line,
                    karta=group_line,
                    color=line_color,
                    opacity=line_opacity,
                    weight=line_weight,
                    popup_dict=line_popup,
                    axis=1,
                )
                group_line.add_to(karta)

            # Handle DataFrame or Point geometry
            else:  # if not isinstance(gdf, gpd.GeoDataFrame) or isinstance(geom, Point):
                if not (heatmap or cluster):
                    group_point = folium.FeatureGroup(name=f"{i}- Point")
                    gdf.apply(
                        Karta._add_point,
                        karta=group_point,
                        color=point_color,
                        speed_field=speed_field,
                        speed_limit_field=speed_limit_field,
                        opacity=point_opacity,
                        radius=point_radius,
                        weight=point_weight,
                        popup_dict=point_popup,
                        x=x,
                        y=y,
                        axis=1,
                    )
                    if point_color == speed_field:
                        template = """
                        {% macro html(this, kwargs) %}

                        <!doctype html>
                        <html lang="en">
                        <head>
                          <meta charset="utf-8">
                          <meta name="viewport" content="width=device-width, initial-scale=1">
                          <title>jQuery UI Draggable - Default functionality</title>
                          <link rel="stylesheet" href="//code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">

                          <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
                          <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>

                          <script>
                          $( function() {
                            $( "#maplegend" ).draggable({
                                            start: function (event, ui) {
                                                $(this).css({
                                                    right: "auto",
                                                    top: "auto",
                                                    bottom: "auto"
                                                });
                                            }
                                        });
                        });

                          </script>
                        </head>
                        <body>


                        <div id='maplegend' class='maplegend'
                            style='position: absolute; z-index:9999; border:2px solid grey; background-color:rgba(255, 255, 255, 1);
                             border-radius:6px; padding: 10px; font-size:14px; right: 95px; top: 10px;'>

                        <!-- The next line is a comment -->
                        <!-- <div class='legend-title'>Legend (draggable!)</div> -->
                        <div class='legend-scale'>
                          <ul class='legend-labels'>
                            <li><span style='background:black;'></span>Speeding ≥ 40% </li>
                            <li><span style='background:red;'></span>30% ≤ Speeding < 40%</li>
                            <li><span style='background:orange;'></span>20% ≤ Speeding < 30%</li>
                            <li><span style='background:yellow;'></span>10% ≤ Speeding < 20%</li>
                            <li><span style='background:green;'></span>0 < Speeding < 10%</li>
                            <li><span style='background:blue;'></span>No speeding</li>
                            <li><span style='background:purple;'  ></span>Speed limit unavailable</li>
                          </ul>
                        </div>
                        </div>

                        </body>
                        </html>

                        <style type='text/css'>
                          .maplegend .legend-title {
                            text-align: left;
                            margin-bottom: 0;
                            font-weight: bold;
                            font-size: 90%;
                            }
                          .maplegend .legend-scale ul {
                            margin: 0;
                            margin-bottom: 0;
                            padding: 0;
                            float: left;
                            list-style: none;
                            }
                          .maplegend .legend-scale ul li {
                            font-size: 80%;
                            list-style: none;
                            margin-left: 0;
                            line-height: 18px;
                            margin-bottom: 1px;
                            }
                          .maplegend ul.legend-labels li span {
                            display: block;
                            float: left;
                            height: 16px;
                            width: 30px;
                            margin-right: 5px;
                            margin-left: 0;
                            border: 1px solid #999;
                            }
                          .maplegend .legend-source {
                            font-size: 80%;
                            color: #777;
                            clear: both;
                            }
                          .maplegend a {
                            color: #777;
                            }
                        </style>
                        {% endmacro %}"""

                        macro = MacroElement()
                        macro._template = Template(template)
                        karta.get_root().add_child(macro)

                    group_point.add_to(karta)

                # Add clustering, heatmap, line connections, and buffer/ring visualizations as specified

                # Create a clustering layer if `cluster=True`
                if cluster:
                    group_cluster = folium.FeatureGroup(name=f"{i}- Cluster")
                    # If the input is a regular DataFrame, use the latitude and longitude columns
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        group_cluster.add_child(plugins.MarkerCluster(locations=list(zip(lats, lons))))
                    # If it's a GeoDataFrame, use geometry coordinates
                    else:
                        group_cluster.add_child(plugins.MarkerCluster(locations=list(zip(gdf.geometry.y, gdf.geometry.x))))
                    # Add the clustering layer to the map
                    group_cluster.add_to(karta)

                # Create a heatmap layer if `heatmap=True`
                if heatmap:
                    group_heatmap = folium.FeatureGroup(name=f"{i}- Heatmap")
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        group_heatmap.add_child(plugins.HeatMap(list(zip(lats, lons)), radius=heatmap_radius))
                    else:
                        group_heatmap.add_child(
                            plugins.HeatMap(list(zip(gdf.geometry.y, gdf.geometry.x)), radius=heatmap_radius)
                        )
                    group_heatmap.add_to(karta)

                # Create a line connection layer if `line=True`
                if line:
                    group_line = folium.FeatureGroup(name=f"{i}- Line")
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        group_line.add_child(
                            folium.PolyLine(list(zip(lats, lons)), color=line_color, weight=line_weight, opacity=line_opacity)
                        )
                    else:
                        group_line.add_child(
                            folium.PolyLine(
                                list(zip(gdf.geometry.y, gdf.geometry.x)),
                                color=line_color,
                                weight=line_weight,
                                opacity=line_opacity,
                            )
                        )
                    # Add the line layer to the map
                    group_line.add_to(karta)

                # Create an animated path layer using AntPath if `antpath=True`
                if antpath:
                    group_antpath = folium.FeatureGroup(name=f"{i}- Ant Path")
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        group_antpath.add_child(plugins.AntPath(list(zip(lats, lons))))
                    else:
                        group_antpath.add_child(plugins.AntPath(list(zip(gdf.geometry.y, gdf.geometry.x))))
                    group_antpath.add_to(karta)

                # Create a buffer visualization if `buffer_radius > 0`
                if buffer_radius > 0:
                    group_buffer = folium.FeatureGroup(name=f"{i}- Buffer")
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        bgdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(lons, lats), crs="EPSG:4326")
                    else:
                        bgdf = gdf.copy()  # buffered gdf: Create a copy of the GeoDataFrame to modify geometries
                    # Apply buffer to geometries using the specified radius in meters
                    bgdf["geometry"] = (
                        bgdf.to_crs(GeomUtils.find_proj(bgdf.geometry.values[0])).buffer(buffer_radius).to_crs("EPSG:4326")
                    )
                    # Add the buffered geometries to the map as polygons
                    bgdf.apply(
                        Karta._add_poly,
                        karta=group_buffer,
                        fill_color=fill_color,
                        highlight_color=fill_color,
                        popup_dict=None,
                        axis=1,
                    )
                    # Add the buffer layer to the map
                    group_buffer.add_to(karta)

                # Create ring visualization if `ring_outer_radius > 0`
                if ring_outer_radius > 0:
                    group_ring = folium.FeatureGroup(name=f"{i}- Ring")
                    if not isinstance(gdf, gpd.GeoDataFrame):
                        bgdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(lons, lats), crs="EPSG:4326")
                    else:
                        bgdf = gdf.copy()  # buffered gdf: Create a copy of the GeoDataFrame to modify geometries
                    # Create ring shapes by applying an outer and inner buffer, subtracting the inner from the outer
                    bgdf["geometry"] = (
                        bgdf.to_crs(GeomUtils.find_proj(bgdf.geometry.values[0]))
                        .buffer(ring_outer_radius)
                        .difference(bgdf.to_crs(GeomUtils.find_proj(bgdf.geometry.values[0])).buffer(ring_inner_radius))
                        .to_crs("EPSG:4326")
                    )  # radius in meters
                    # Add the ring-shaped geometries to the map as polygons
                    bgdf.apply(
                        Karta._add_poly,
                        karta=group_ring,
                        fill_color=fill_color,
                        highlight_color=fill_color,
                        popup_dict=None,
                        axis=1,
                    )
                    # Add the ring layer to the map
                    group_ring.add_to(karta)

        # Geohash visualization if `geohash_res > 0`
        if geohash_res > 0:  # inner=False doesn't work if compact=True
            # Create a polygon for bounding box if input is not a polygon
            if isinstance(geom, Polygon) or isinstance(geom, MultiPolygon):
                cdf = gdf.copy()
            else:
                bb = Polygon([[minlon, minlat], [maxlon, minlat], [maxlon, maxlat], [minlon, maxlat], [minlon, minlat]])
                cdf = gpd.GeoDataFrame({"geometry": [bb]}, crs="EPSG:4326")  # Create a bounding box GeoDataFrame

            # Convert geometries to geohash cells and their geometries
            cells, _ = SpatialIndex.ppoly_cell(cdf, cell_type="geohash", res=geohash_res, compact=compact)
            geoms, res = SpatialIndex.cell_poly(cells, cell_type="geohash")
            cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

            # Add geohash cells to the map as a polygon layer
            group_geohash = folium.FeatureGroup(name=f"{i} - Geohash")

            cdf.apply(
                Karta._add_poly,
                karta=group_geohash,
                fill_color=fill_color,
                highlight_color=highlight_color,
                popup_dict={"ID": "id", "Resolution": "res"},
                axis=1,
            )
            group_geohash.add_to(karta)

        # S2 cell visualization if `s2_res > -1`
        if s2_res > -1:
            # Create a polygon for bounding box if input is not a polygon
            if isinstance(geom, Polygon) or isinstance(geom, MultiPolygon):
                cdf = gdf.copy()
            else:
                bb = Polygon([[minlon, minlat], [maxlon, minlat], [maxlon, maxlat], [minlon, maxlat], [minlon, minlat]])
                cdf = gpd.GeoDataFrame({"geometry": [bb]}, crs="EPSG:4326")  # cell df

            # Convert geometries to S2 cells and their geometries
            cells, _ = SpatialIndex.ppoly_cell(cdf, cell_type="s2", res=s2_res, compact=compact)
            geoms, res = SpatialIndex.cell_poly(cells, cell_type="s2")
            cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

            # Add S2 cells to the map as a polygon layer
            group_s2 = folium.FeatureGroup(name=f"{i} - S2")
            cdf.apply(
                Karta._add_poly,
                karta=group_s2,
                fill_color=fill_color,
                highlight_color=highlight_color,
                popup_dict={"ID": "id", "Resolution": "res"},
                axis=1,
            )
            group_s2.add_to(karta)

        # H3 cell visualization if `h3_res > -1`
        if h3_res > -1:
            if isinstance(geom, Polygon) or isinstance(geom, MultiPolygon):
                cdf = gdf.copy()
            # Create a bounding box GeoDataFrame
            else:
                bb = Polygon([[minlon, minlat], [maxlon, minlat], [maxlon, maxlat], [minlon, maxlat], [minlon, minlat]])
                cdf = gpd.GeoDataFrame({"geometry": [bb]}, crs="EPSG:4326")  # cell df

            # Convert geometries to H3 cells and their Shapely hexagons
            cells, _ = SpatialIndex.ppoly_cell(
                cdf, cell_type="h3", res=h3_res, force_full_cover=force_full_cover, compact=compact
            )
            geoms, res = SpatialIndex.cell_poly(cells, cell_type="h3")
            cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

            # Add H3 cells to the map as a polygon layer
            group_h3 = folium.FeatureGroup(name=f"{i} - H3")
            cdf.apply(
                Karta._add_poly,
                karta=group_h3,
                fill_color=fill_color,
                highlight_color=highlight_color,
                popup_dict={"ID": "id", "Resolution": "res"},
                axis=1,
            )
            group_h3.add_to(karta)
        folium.LayerControl(collapsed=False).add_to(karta)
        return karta

    @staticmethod
    def choropleth(mdf: gpd.GeoDataFrame, columns: list, legend: str, bins: list = None, palette: str = "YlOrRd") -> folium.Map:
        """
        Creates a choropleth map using the given GeoDataFrame and specified parameters.

        This function generates a Folium choropleth map layer by visualizing the data from a GeoDataFrame using color gradients
        to represent different data values across geographic areas.

        Parameters
        ----------
        mdf : geopandas.GeoDataFrame
            The GeoDataFrame containing multipolygon geometries and data attributes to be visualized.
        columns : list of str
            A list of two elements:
                - columns[0] : str
                    The column name in `mdf` that contains unique identifiers for each region.
                - columns[1] : str
                    The column name in `mdf` containing the data values to be visualized.
        legend : str
            The title for the legend, which describes what is represented on the map.

        bins : list of float, optional
            A list of numerical values that define the value intervals for the choropleth color categories.
            By default, the bins correspond to quantiles at [0, 0.25, 0.5, 0.75, 0.98, 1.0].
        palette : str, optional
            The color palette to be used for the choropleth (default is "YlOrRd").

        Returns
        -------
        folium.Map
            The Folium map object containing the choropleth layer.

        Examples
        --------
        >>> choropleth(
                mdf,
                ['region_id', 'population'],
                bins=[0, 100, 500, 1000, 5000],
                legend="Population by Region",
                palette="YlOrRd",
            )
        """
        # Extract the bounding coordinates of the GeoDataFrame
        minlon, minlat, maxlon, maxlat = mdf.total_bounds  # Get the total bounds of the GeoDataFrame
        sw = [minlat, minlon]  # South-west corner
        ne = [maxlat, maxlon]  # North-east corner
        karta = Karta._base_map(sw, ne)  # Create a base map using the bounding coordinates

        if bins is None:
            bins = np.quantile(mdf[columns[1]].dropna(), [0, 0.25, 0.5, 0.75, 0.98, 1])
        # Create a choropleth layer based on the GeoDataFrame
        choropleth = folium.Choropleth(
            geo_data=mdf,  # The GeoDataFrame containing geographic data
            name="Choropleth",  # Name of the layer for display in layer control
            data=mdf,  # The data source for values to be represented
            columns=columns,  # [unique_identifier_column, data_value_column] for matching regions with data
            key_on="feature.properties." + columns[0],  # Key to match GeoDataFrame regions with the data
            legend_name=legend,  # Description of the data being visualized
            bins=bins,  # Value ranges for choropleth colors
            fill_color=palette,  # Color scheme for the choropleth
            fill_opacity=0.5,  # Transparency level of filled regions
            line_opacity=0.25,  # Transparency level of borders between regions
            smooth_factor=0,  # Level of smoothing applied to the edges of regions
            highlight=True,  # Enable or disable highlighting of regions on hover
        ).add_to(karta)  # Add the choropleth layer to the map

        # Add a tooltip to display the attribute values for each region when hovered over
        folium.features.GeoJsonTooltip(fields=columns).add_to(choropleth.geojson)

        # Add layer control to the map
        folium.LayerControl(collapsed=False).add_to(karta)

        # Return the Folium map object containing the choropleth layer
        return karta

    @staticmethod
    def joint_choropleth(
        mdf: gpd.GeoDataFrame,
        gdf: gpd.GeoDataFrame,
        poly_id: str,
        legend: str,
        bins: list = None,
        palette: str = "YlOrRd",
        cell_type: str = None,
        res: int = None,
    ) -> folium.Map:
        """
        Creates a joint choropleth map by counting point features (gdf) within polygons (mdf).

        This function generates a Folium choropleth map showing the count of point features
        from gdf that fall within each polygon of mdf. Optionally converts mdf to spatial
        cells (geohash, S2, or H3) before counting.

        Parameters
        ----------
        mdf : gpd.GeoDataFrame
            GeoDataFrame containing polygon geometries to aggregate points into.
        gdf : gpd.GeoDataFrame
            GeoDataFrame containing geometry features to count within polygons.
        poly_id : str
            Column name in mdf containing unique polygon identifiers.
        legend : str
            Title for the choropleth legend.
        bins : list, optional
            Value intervals for choropleth color categories.
            By default, the bins correspond to quantiles at [0, 0.25, 0.5, 0.75, 0.98, 1.0].
        palette : str, optional
            Color palette name (default: "YlOrRd").
        cell_type : str, optional
            Spatial index type ("geohash", "s2", or "h3") to convert mdf polygons into.
        res : int, optional
            Resolution level for spatial indexing if cell_type is specified.

        Returns
        -------
        folium.Map
            Folium map object with choropleth layer.

        Examples
        --------
        >>> joint_choropleth(
               neighborhoods,  # Polygon GeoDataFrame
               restaurants,    # Point GeoDataFrame
               poly_id="neighborhood_id",
               legend="Restaurant Count",
               cell_type="h3",
               res=8
           )
        """
        if cell_type is not None:
            cell_list, _ = SpatialIndex.ppoly_cell(mdf, cell_type=cell_type, res=res)
            cells, _ = SpatialIndex.pcell_poly(cell_list, cell_type=cell_type)
            mdf = gpd.GeoDataFrame({poly_id: cell_list}, geometry=cells, crs="EPSG:4326")

        gdf = gpd.sjoin(gdf[["geometry"]], mdf[[poly_id, "geometry"]], predicate="within")
        del gdf["index_right"]

        mdf = pd.merge(mdf, gdf[poly_id].value_counts().reset_index(), on=poly_id, how="left")

        if bins is None:
            bins = np.quantile(mdf["count"].dropna(), [0, 0.25, 0.5, 0.75, 0.98, 1])
        return Karta.choropleth(
            mdf,
            columns=[poly_id, "count"],
            bins=bins,
            legend=legend,
            palette=palette,
        )


class SnabbKarta:
    @staticmethod
    def _get_color(col: int | float | str) -> list:
        """
        Generates a consistent color based on input value.
        Returns as [R, G, B] list for deck.gl compatibility.
        """
        palette = [
            [255, 0, 0],  # red
            [0, 128, 0],  # green
            [0, 0, 255],  # blue
            [255, 255, 0],  # yellow
            [128, 0, 128],  # purple
            [255, 165, 0],  # orange
            [255, 192, 203],  # pink
            [165, 42, 42],  # brown
            [128, 128, 128],  # gray
            [255, 215, 0],  # gold
            [0, 255, 255],  # cyan
            [255, 0, 255],  # magenta
            [0, 255, 0],  # lime
            [0, 0, 128],  # navy
            [0, 128, 128],  # teal
        ]

        if col is None or (isinstance(col, float) and math.isnan(col)):
            return [0, 0, 0]

        if isinstance(col, (int, float)):
            idx = int(col) % len(palette)
        else:
            col = str(col)
            col = re.sub(r"[\W_]+", "", col)
            idx = int(col, 36) % len(palette)

        return palette[idx]

    @staticmethod
    def _get_speed_colors(speeds, speed_limits, opacity):
        # Create result array filled with default values (black)
        result = np.full((len(speeds), 4), [0, 0, 0, opacity], dtype=np.uint8)

        # Create masks for each condition
        invalid_mask = pd.isna(speed_limits) | (speed_limits <= 0)
        within_limit_mask = ~invalid_mask & (speeds <= speed_limits)
        green_mask = ~invalid_mask & (speeds < 1.1 * speed_limits) & (speeds > speed_limits)
        yellow_mask = ~invalid_mask & (speeds < 1.2 * speed_limits) & (speeds >= 1.1 * speed_limits)
        orange_mask = ~invalid_mask & (speeds < 1.3 * speed_limits) & (speeds >= 1.2 * speed_limits)
        red_mask = ~invalid_mask & (speeds < 1.4 * speed_limits) & (speeds >= 1.3 * speed_limits)

        # Apply colors based on masks (all values as uint8)
        result[invalid_mask] = [128, 0, 128, opacity]  # purple
        result[within_limit_mask] = [0, 0, 255, opacity]  # blue
        result[green_mask] = [0, 255, 0, opacity]  # green
        result[yellow_mask] = [255, 255, 0, opacity]  # yellow
        result[orange_mask] = [255, 165, 0, opacity]  # orange
        result[red_mask] = [255, 0, 0, opacity]  # red

        return result

    @staticmethod
    def _create_point_layer(
        gdf: gpd.GeoDataFrame,
        color: str = "blue",
        opacity: float = 1,
        get_radius: str | int = 1,
        radius_min_pixels: int = 1,
        radius_max_pixels: int = 10,
        speed_field: str = "speed",
        speed_limit_field: str = "speedlimit",
        pickable: bool = True,
    ) -> lb.ScatterplotLayer:
        """Creates a Lonboard ScatterplotLayer from a GeoDataFrame."""

        # Convert opacity to uint8 (0-255 range)
        opacity = np.uint8(opacity * 255)

        # Handle color assignment
        if color == speed_field:
            fill_color = SnabbKarta._get_speed_colors(gdf[speed_field].values, gdf[speed_limit_field].values, opacity)

        elif color in gdf.columns:
            fill_color = np.array([SnabbKarta._get_color(item) + [opacity] for item in gdf[color]], dtype=np.uint8)

        else:
            rgb_color = [int(c * 255) for c in matplotlib.colors.to_rgb(color)]
            fill_color = np.array([rgb_color + [opacity]] * len(gdf), dtype=np.uint8)

        # Handle radius
        radius = gdf[get_radius].values.astype(float) if get_radius in gdf.columns else np.array([get_radius] * len(gdf))

        return lb.ScatterplotLayer.from_geopandas(
            gdf,
            get_fill_color=fill_color,
            get_radius=radius,
            get_line_color=fill_color,
            filled=True,
            stroked=True,
            line_width_min_pixels=1,
            radius_min_pixels=radius_min_pixels,
            radius_max_pixels=radius_max_pixels,
            pickable=pickable,
        )

    @staticmethod
    def _create_line_layer(
        gdf: gpd.GeoDataFrame,
        line_color="blue",
        opacity: float = 0.5,
        pickable=True,
    ) -> lb.PathLayer:
        # Convert opacity to 0-255 range
        opacity = int(opacity * 255)

        if line_color in gdf.columns:
            get_color = np.array([SnabbKarta._get_color(item) + [opacity] for item in gdf[line_color]], dtype=np.uint8)
        else:
            rgb_color = [int(c * 255) for c in matplotlib.colors.to_rgb(line_color)]
            get_color = np.array([[*rgb_color, opacity]] * len(gdf), dtype=np.uint8)

        return lb.PathLayer.from_geopandas(
            gdf,
            get_color=get_color,
            auto_highlight=True,
            highlight_color=[255, 0, 0],
            width_min_pixels=1,
            width_max_pixels=10,
            pickable=pickable,
        )

    @staticmethod
    def _create_poly_layer(
        gdf: gpd.GeoDataFrame,
        fill_color="red",
        opacity: float = 0.5,
        poly_highlight=True,
        pickable=True,
    ) -> lb.PolygonLayer:
        # Convert opacity to 0-255 range
        opacity = int(opacity * 255)

        if fill_color in gdf.columns:
            get_fill_color = np.array([SnabbKarta._get_color(item) + [opacity] for item in gdf[fill_color]], dtype=np.uint8)
        else:
            rgb_color = [int(c * 255) for c in matplotlib.colors.to_rgb(fill_color)]
            get_fill_color = np.array([[*rgb_color, opacity]] * len(gdf), dtype=np.uint8)

        return lb.PolygonLayer.from_geopandas(
            gdf,
            get_fill_color=get_fill_color,
            auto_highlight=poly_highlight,
            highlight_color=[0, 255, 0, 255],
            get_line_color=[0, 0, 0, 255],
            line_width_min_pixels=1,
            line_width_max_pixels=1,
            pickable=pickable,
        )

    @staticmethod
    def data_to_gdf(
        data: gpd.GeoDataFrame | pd.DataFrame | set,
        geom_type: str | None = None,  #  'h3', 's2', 'geohash', 'osm', 'uprn', 'usrn', 'postcode'
        # Supported types: geospatial cell IDs (geohash, s2, h3),
        # OSM ID, UPRN, USRN, and postcode.
        geom_col: str | list[str] | None = None,
        # e.g. ['northing', 'easting'], 'h3_8', 'osm_id', 'uprn',  'postcode', 'postcode_sec'
        crs: int = 4326,  # CRS of geom_col
        lookup_gdf: pd.DataFrame | gpd.GeoDataFrame | None = None,  # external df containing geometry
        lookup_key: str = None,  # geometry column name in lookup_gdf
        # lat_col: str = "lat",  # Latitude column name in lookup_gdf
        # lon_col: str = "lon",  # Longitude column name in lookup_gdf
        osm_url: str | None = "https://overpass-api.de/api/interpreter",  # OpenStreetMap server URL
    ) -> gpd.GeoDataFrame:
        # if data is a gpd.GeoDataFrame, simply copy it
        if isinstance(data, gpd.GeoDataFrame):
            gdf = data.copy()
            return gdf

        # Convert pd.DataFrame to gpd.GeoDataFrame
        elif isinstance(data, pd.DataFrame):
            if geom_type in ["geohash", "s2", "s2_int", "h3"]:
                df = data.dropna(subset=geom_col)
                cells = df[geom_col].values
                geoms, res = SpatialIndex.cell_poly(cells, cell_type=geom_type)
                gdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")
            elif geom_type in ["uprn", "usrn", "postcode", "osm"]:
                df = data.dropna(subset=geom_col)
                df = df.sort_values(by=geom_col).reset_index(drop=True)
                gdf = lookup_gdf.merge(df, left_on=lookup_key, right_on=geom_col, how="right")
            else:  # if geom_type == None, find lat, lon columns
                df = data.copy()
                if geom_col:  # if geom_col provided
                    x, y = geom_col[1], geom_col[0]
                else:  # if geom_col = None determine lat, lon columns
                    x = [col for col in df.columns if "lon" in col.lower() or "lng" in col.lower()][0]
                    y = [col for col in df.columns if "lat" in col.lower()][0]
                gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df[x], df[y]), crs=crs).to_crs(4326)
            return gdf

        # Convert set to gpd.GeoDataFrame
        elif isinstance(data, set):  # e.g. set of uprn {1, 26, 27, 30, 31}
            data = [item for item in data if item is not None]  # Remove None from the set if exists
            # Convert geospatial cells to Shapely geometries
            if geom_type in ["geohash", "s2", "s2_int", "h3"]:
                cells = data
                geoms, res = SpatialIndex.cell_poly(cells, cell_type=geom_type)
                gdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")
            # Convert other types to Shapely geometries
            elif geom_type in ["uprn", "usrn", "postcode", "osm"]:
                df = pd.DataFrame({geom_type: sorted(data)})  # Sort data for faster join
                gdf = lookup_gdf.merge(df, left_on=lookup_key, right_on=geom_type, how="right")
            return gdf
        # Unsupported input type
        else:
            raise TypeError("Invalid input type: data must be gpd.GeoDataFrame, pd.DataFrame or a set")

    @staticmethod
    def plp(
        data_list: gpd.GeoDataFrame | pd.DataFrame | set | list[gpd.GeoDataFrame | pd.DataFrame | set],
        geom_type: str | None = None,  #  'h3', 's2', 'geohash', 'osm', 'uprn', 'usrn', 'postcode'
        # Supported types: geospatial cell IDs (geohash, s2, h3),
        # OSM ID, UPRN, USRN, and postcode.
        geom_col: str | list[str] | None = None,
        # e.g. ['northing', 'easting'], 'h3_8', 'osm_id', 'uprn',  'postcode', 'postcode_sec'
        crs: int = 4326,  # CRS of geom_col
        lookup_gdf: pd.DataFrame | gpd.GeoDataFrame | None = None,  # external df containing geometry
        lookup_key: str = None,  # geometry column name in lookup_gdf
        # lat_col: str = "lat",  # Latitude column name in lookup_gdf
        # lon_col: str = "lon",  # Longitude column name in lookup_gdf
        osm_url: str | None = "https://overpass-api.de/api/interpreter",  # OpenStreetMap server URL
        # Map tiles
        tiles: str = CartoStyle.Positron,  # DarkMatter
        pitch: int = 30,
        map_height: int = 800,
        # Point
        cluster: bool = False,
        heatmap: bool = False,
        point_color: str = "blue",
        point_opacity: float = 1,
        speed_field: str = "speed",
        speed_limit_field: str = "speedlimit",
        point_radius: int | str = 1,
        radius_min_pixels: int = 1,
        radius_max_pixels: int = 10,
        buffer_radius: int = 0,
        ring_inner_radius: int = 0,
        ring_outer_radius: int = 0,
        # LineString
        line_color: str = "blue",
        line_opacity: float = 0.5,
        # Polygon
        centroid: bool = False,  # if True it shows centroids of polygons on the map.
        fill_color: str = "red",
        highlight_color: str = "green",
        fill_opacity: float = 0.25,
        highlight_opacity: float = 0.5,
        geohash_res: int = 0,
        s2_res: int = -1,
        h3_res: int = -1,
        force_full_cover: bool = True,
        geohash_inner: bool = False,
        compact: bool = False,
    ) -> lb.Map:
        minlat, maxlat, minlon, maxlon = 90, -90, 180, -180
        layers = []

        # Ensure `data_list` is always a list (of gpd.GeoDataFrames, df.DataFrames or dict)
        data_list = data_list if isinstance(data_list, list) else [data_list]
        # Iterate through each set, pd.DataFrame or gpd.GeoDataFrame in the list to add layers to the map
        for data in data_list:
            gdf = SnabbKarta.data_to_gdf(data, geom_type, geom_col, crs, lookup_gdf, lookup_key)
            # Update overall bounding box
            gminlon, gminlat, gmaxlon, gmaxlat = gdf.total_bounds  # gminlon: gdf minlon
            minlat, minlon = min(minlat, gminlat), min(minlon, gminlon)  # minlat: total minlat
            maxlat, maxlon = max(maxlat, gmaxlat), max(maxlon, gmaxlon)

            # Create layers
            for geom in gdf.geometry.type.unique():
                sdf = gdf[gdf.geometry.type == geom]  # subset gdf
                if geom == "Point":
                    layers.append(
                        SnabbKarta._create_point_layer(
                            sdf,
                            color=point_color,
                            opacity=point_opacity,
                            speed_field=speed_field,
                            speed_limit_field=speed_limit_field,
                            get_radius=point_radius,
                        )
                    )
                elif geom in ("LineString", "MultiLineString"):
                    layers.append(SnabbKarta._create_line_layer(sdf, line_color=line_color))
                elif geom in ("Polygon", "MultiPolygon"):
                    layers.append(SnabbKarta._create_poly_layer(sdf, fill_color=fill_color))

                # Show centroids of the geometry if `centroid=True`
                if centroid:
                    cdf = gpd.GeoDataFrame({"geometry": gdf.centroid}, crs=gdf.crs)  # centroid df
                    centroid_layer = SnabbKarta._create_point_layer(cdf, get_radius=1000)
                    layers.append(centroid_layer)

                # Create a buffer layer if `buffer_radius > 0`
                if buffer_radius > 0:
                    bgdf = gdf[["geometry"]]  # buffered gdf: Create a copy of the GeoDataFrame to modify geometries
                    # Apply buffer to geometries using the specified radius in meters
                    bgdf["geometry"] = bgdf.to_crs(GeomUtils.find_proj(geom)).buffer(buffer_radius).to_crs("EPSG:4326")
                    # Add the buffer layer to the map
                    buffer_layer = SnabbKarta._create_poly_layer(bgdf)
                    layers.append(buffer_layer)

                # Create ring layer if `ring_outer_radius > 0`
                if ring_outer_radius > 0:
                    bgdf = gdf[["geometry"]]  # buffered gdf: Create a copy of the GeoDataFrame to modify geometries
                    # Create ring shapes by applying an outer and inner buffer, subtracting the inner from the outer
                    bgdf["geometry"] = (
                        bgdf.to_crs(GeomUtils.find_proj(geom))
                        .buffer(ring_outer_radius)
                        .difference(bgdf.to_crs(GeomUtils.find_proj(geom)).buffer(ring_inner_radius))
                        .to_crs("EPSG:4326")
                    )  # radius in meters
                    # Add the ring-shaped geometries to the map as polygons
                    ring_layer = SnabbKarta._create_poly_layer(bgdf)
                    layers.append(ring_layer)

                # Geohash visualization if `geohash_res > 0`
                if geohash_res > 0:  # inner=False doesn't work if compact=True
                    # Create a polygon for bounding box if input is not a polygon
                    if geom in ("Polygon", "MultiPolygon"):
                        cdf = gdf[["geometry"]]
                    else:
                        # Get the concave hull with a ratio parameter (0-1)
                        # Smaller ratio = tighter fit, larger ratio = more convex
                        tight_polygon = shapely.concave_hull(gdf.geometry.unary_union, ratio=0.1)
                        cdf = gpd.GeoDataFrame(geometry=[tight_polygon], crs=gdf.crs)

                    # Convert geometries to geohash cells and their geometries
                    cells, _ = SpatialIndex.ppoly_cell(cdf, cell_type="geohash", res=geohash_res, compact=compact)
                    geoms, res = SpatialIndex.cell_poly(cells, cell_type="geohash")
                    cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

                    geohash_layer = SnabbKarta._create_poly_layer(cdf, fill_color="green")
                    layers.append(geohash_layer)

                # S2 cell visualization if `s2_res > -1`
                if s2_res > -1:
                    # Create a polygon for bounding box if input is not a polygon
                    if geom in ("Polygon", "MultiPolygon"):
                        cdf = gdf[["geometry"]]
                    else:
                        tight_polygon = shapely.concave_hull(gdf.geometry.unary_union, ratio=0.1)
                        cdf = gpd.GeoDataFrame(geometry=[tight_polygon], crs=gdf.crs)

                    # Convert geometries to S2 cells and their geometries
                    cells, _ = SpatialIndex.ppoly_cell(cdf, cell_type="s2", res=s2_res, compact=compact)
                    geoms, res = SpatialIndex.cell_poly(cells, cell_type="s2")
                    cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

                    s2_layer = SnabbKarta._create_poly_layer(cdf, fill_color="green")
                    layers.append(s2_layer)

                # H3 cell visualization if `h3_res > -1`
                if h3_res > -1:
                    if geom in ("Polygon", "MultiPolygon"):
                        cdf = gdf[["geometry"]]
                    else:
                        tight_polygon = shapely.concave_hull(gdf.geometry.unary_union, ratio=0.1)
                        cdf = gpd.GeoDataFrame(geometry=[tight_polygon], crs=gdf.crs)

                    # Convert geometries to H3 cells and their Shapely hexagons
                    cells, _ = SpatialIndex.ppoly_cell(
                        cdf, cell_type="h3", res=h3_res, force_full_cover=force_full_cover, compact=compact
                    )
                    geoms, res = SpatialIndex.cell_poly(cells, cell_type="h3")
                    cdf = gpd.GeoDataFrame({"id": cells, "res": res, "geometry": geoms}, crs="EPSG:4326")

                    h3_layer = SnabbKarta._create_poly_layer(cdf, fill_color="green")
                    layers.append(h3_layer)
            # for geom in gdf.geometry.type.unique():
        # for gdf in data_list:

        # Create a base map using the bounding box
        sw = [minlat, minlon]  # South West (bottom left corner)
        ne = [maxlat, maxlon]  # North East (top right corner)

        # Calculate center and zoom if not provided
        lat_center, lon_center = (sw[0] + ne[0]) / 2, (sw[1] + ne[1]) / 2
        max_length = max(ne[0] - sw[0], ne[1] - sw[1])  # max(delta_lat, delta_lon)
        # Adjust zoom baseline depending on map extent
        if max_length < 0:  # if no data available, show lat, lon = (0, 0)
            zoom = 5
        elif max_length == 0:  # if only one point available, set the zoom level to 20
            zoom = 20
        elif max_length > 5:  # large area (e.g. whole UK)
            zoom = 12 - math.log(max_length * 2, 1.5)
        else:  # smaller area (e.g. London)
            zoom = 11 - math.log(max_length * 2, 1.5)

        return lb.Map(
            layers=layers,
            basemap_style=tiles,
            height=map_height,
            view_state={
                "longitude": lon_center,
                "latitude": lat_center,
                "zoom": zoom,
                "pitch": pitch,
                "bearing": 0,
            },
        )


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
    def find_proj(geom: BaseGeometry) -> str:
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
    def get_coords(gdf: pd.DataFrame | gpd.GeoDataFrame, coord_cols: list = None) -> list[list[int]]:
        """
        Extract coordinate pairs from spatial data.

        Handles both GeoDataFrames (geometry-based) and DataFrames (column-based).
        For GeoDataFrames: extracts point coordinates or centroids for non-point geometries.
        For DataFrames: extracts coordinates from specified or auto-detected columns.

        Parameters
        ----------
        gdf : pandas.DataFrame or geopandas.GeoDataFrame
            Input data containing spatial information. Must have either:
            - For GeoDataFrame: a 'geometry' column with shapely geometries
            - For DataFrame: columns containing x/y coordinate values

        coord_cols : list of str, optional
            For DataFrames only: list of two column names [lat, lon].
            If None, attempts to auto-detect lat/lon columns.
            Ignored for GeoDataFrames.

        Returns
        -------
        list[list[int]]
            List of coordinate pairs [x, y]. Each sublist contains two numeric values.
            For point geometries: direct coordinates.
            For non-point geometries: centroid coordinates.
            For DataFrames: values from the specified coordinate columns.

        Notes
        -----
        - The return type annotation `list[list[int]]` is nominal; actual return type
          is typically `numpy.ndarray` with float values.
        - For GeoDataFrames: returns centroids for non-point geometries (Polygon,
          LineString, etc.).
        - Auto-detection for DataFrames searches for 'lon'/'lng' and 'lat' substrings
          (case-insensitive).
        """
        # Check if input is a GeoDataFrame (geometry-based)
        if isinstance(gdf, gpd.GeoDataFrame):
            # Use .geom_type to check geometry types efficiently
            geom_types = set(gdf.geometry.geom_type)

            # Ternary-like logic: if only points, use direct coordinates; else use centroids
            if geom_types == {"Point"}:
                # All geometries are points - fastest approach
                return np.column_stack((gdf.geometry.x, gdf.geometry.y))
            else:
                # Mixed or non-point geometries - use centroids
                return np.column_stack((gdf.geometry.centroid.x, gdf.geometry.centroid.y))
        else:
            # DataFrame handling (column-based coordinates)
            if coord_cols:  # if coordinate columns explicitly provided
                x, y = coord_cols[1], coord_cols[0]  # x = second element (longitude), y = first element (latitude)
            else:  # if no columns specified, auto-detect
                # Find longitude column (search for 'lon' or 'lng' in column names)
                x = [col for col in gdf.columns if "lon" in col.lower() or "lng" in col.lower()][0]
                # Find latitude column (search for 'lat' in column names)
                y = [col for col in gdf.columns if "lat" in col.lower()][0]

            # Extract coordinate pairs from DataFrame columns and convert to numpy array
            return gdf[[x, y]].to_numpy()

    @staticmethod
    def geom_stats(
        geom: Polygon | MultiPolygon | None = None, projection: str | None = None, unit: str = "m"
    ) -> list[int | float] | None:
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
    def flatten_3d(geoms: gpd.GeoSeries) -> gpd.GeoSeries:
        """
        Flattens a GeoSeries of 3D geometries into 2D geometries.

        This function removes the z-coordinate from each geometry in the input GeoSeries,
        converting it into its 2D equivalent. Non-3D geometries are preserved.

        Parameters
        ----------
        geom : gpd.GeoSeries
            A GeoSeries containing geometries (with or without z-coordinates).

        Returns
        -------
        gpd.GeoSeries
            A GeoSeries of 2D geometries.

        Notes
        -----
        - Preserves polygon holes (interior rings).
        - Keeps non-3D geometries unchanged.
        - Supports Point, LineString, Polygon, MultiPoint, MultiLineString,
          MultiPolygon, and GeometryCollection.
        """

        def drop_z(geom):
            if geom is None or not getattr(geom, "has_z", False):
                return geom

            if isinstance(geom, Point):
                return Point(geom.x, geom.y)

            elif isinstance(geom, LineString):
                return LineString([(x, y) for x, y, _ in geom.coords])

            elif isinstance(geom, Polygon):
                exterior = [(x, y) for x, y, _ in geom.exterior.coords]
                interiors = [[(x, y) for x, y, _ in ring.coords] for ring in geom.interiors]
                return Polygon(exterior, interiors)

            elif isinstance(geom, MultiPoint):
                return MultiPoint([drop_z(g) for g in geom.geoms])

            elif isinstance(geom, MultiLineString):
                return MultiLineString([drop_z(g) for g in geom.geoms])

            elif isinstance(geom, MultiPolygon):
                return MultiPolygon([drop_z(g) for g in geom.geoms])

            elif isinstance(geom, GeometryCollection):
                return GeometryCollection([drop_z(g) for g in geom.geoms])

            else:
                return geom  # Unknown type: leave unchanged

        return [drop_z(geom) for geom in geoms]

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
    def h3_stats(geom: BaseGeometry, h3_res: int, compact: bool = False) -> tuple[int, float]:
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
        cells = SpatialIndex.poly_cell([geom], cell_type="h3", res=h3_res)
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
    def way_to_geom(way_id: int, osm_url: str = "https://overpass-api.de/api/interpreter") -> LineString | Polygon:
        """
        Converts an OSM way ID into a Shapely Polygon or LineString object.

        This function retrieves the geometry corresponding to the given OSM way ID and
        returns it as a Shapely `Polygon` or `LineString` object based on whether the way
        forms a closed loop or not.

        Parameters
        ----------
        way_id : int
            The OpenStreetMap (OSM) way ID to be retrieved.
        osm_url : str, optional
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
        >>> osm_url = "https://overpass-api.de/api/interpreter"
        >>> geometry = way_to_geom(way_id, osm_url)
        >>> print(geometry)
        POLYGON ((13.3888 52.5170, 13.3976 52.5291, 13.4286 52.5232, 13.3888 52.5170))
        """
        query = f"[out:json][timeout:600][maxsize:4073741824];way({way_id});out geom;"
        response = requests.get(osm_url, params={"data": query}).json()
        response = response["elements"][0]
        geom = response["geometry"]
        coords = [(node["lon"], node["lat"]) for node in geom]
        if geom[0] == geom[-1]:  # Check if the way forms a closed loop
            return Polygon(coords)
        else:
            return LineString(coords)

    @staticmethod
    def ways_to_geom(ids: list[int], osm_url: str = "https://overpass-api.de/api/interpreter") -> list[LineString | Polygon]:
        """
        Converts an array of OpenStreetMap (OSM) way IDs into Shapely geometries.

        This function retrieves the geometries corresponding to the given OSM way IDs and
        returns a list of Shapely `LineString` or `Polygon` objects based on the geometries
        fetched from the OSM API.

        Parameters
        ----------
        ids : list of int
            A list of OSM way IDs to be retrieved.
        osm_url : str, optional
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
        >>> osm_url = "https://overpass-api.de/api/interpreter"
        >>> geometries = ways_to_geom(way_ids, osm_url)
        >>> print(geometries)
        [<shapely.geometry.polygon.Polygon object at 0x...>,
         <shapely.geometry.linestring.LineString object at 0x...>]
        """
        query = "[out:json][timeout:600][maxsize:4073741824];"
        for item in ids:
            query += f"way({item});out geom;"

        response = requests.get(osm_url, params={"data": query}).json()
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
    def map_matching(df: pd.DataFrame, cost: str, url: str, format: str = "osrm") -> dict | None:
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
    point_cell(lats, lons, cell_type, res):
        Converts latitude and longitude coordinates into spatial index representations.

    poly_cell(geoms, cell_type, res, dump_path):
        Converts a list of geometries into a set of unique spatial cells.

    ppoint_cell(lats, lons, cell_type, res):
        Converts latitude and longitude coordinates into spatial index representations in parallel.

    ppoly_cell(mdf, cell_type, res, compact, dump_path, verbose):
        Performs parallelized conversion of geometries in a GeoDataFrame to cell identifiers.

    cell_point(cells, cell_type):
        Converts a list of cell IDs into their corresponding centroids.

    cell_poly(cells, cell_type):
        Converts a list of spatial cells to their corresponding geometries and resolution levels.

    pcell_point(cells, cell_type):
        Converts a list of cell IDs into their corresponding latitude and longitude points in parallel.

    pcell_poly(cells, cell_type):
        Parallelized version of `cell_poly`, converting a list of spatial cells to geometries and resolution levels.
    """

    @staticmethod
    def geohash_to_poly(geohash: str) -> Polygon:
        """
        Convert a geohash string into a Shapely polygon representing its bounding box.

        Parameters
        ----------
        geohash : str
            A valid geohash string to decode.

        Returns
        -------
        polygon : shapely.geometry.Polygon
            A Shapely Polygon representing the bounding box of the geohash.
            The polygon is defined as a rectangle covering the geohash area.

        Examples
        --------
        >>> from shapely.geometry import Polygon
        >>> poly = SpatialIndex.geohash_to_poly("9q8yyzd")
        >>> isinstance(poly, Polygon)
        True
        >>> list(poly.exterior.coords)[0]
        (-122.421875, 37.7734375)
        """
        lat, lon, lat_err, lon_err = pygeohash.decode_exactly(geohash)
        return box(lon - lon_err, lat - lat_err, lon + lon_err, lat + lat_err)

    @staticmethod
    def geohash_adjacent(geohash: str, direction: str) -> str:
        """
        Compute the adjacent geohash in a specified cardinal direction.

        This function calculates the neighboring geohash for a given geohash string
        in one of the four cardinal directions: north (`'n'`), south (`'s'`), east (`'e'`), or west (`'w'`).
        It handles edge cases where the geohash is on the border of a cell.

        Parameters
        ----------
        geohash : str
            A valid geohash string. Case-insensitive.
        direction : str
            The direction in which to find the adjacent geohash. Must be one of:
            `'n'` (north), `'s'` (south), `'e'` (east), `'w'` (west).

        Returns
        -------
        adjacent : str
            The geohash string representing the adjacent cell in the specified direction.

        Raises
        ------
        ValueError
            If the direction is not one of `'n'`, `'s'`, `'e'`, or `'w'`.

        Examples
        --------
        >>> SpatialIndex.geohash_adjacent("u4pruydqqvj", "n")
        'u4pruydqqvm'

        >>> SpatialIndex.geohash_adjacent("ezs42", "e")
        'ezs48'

        >>> SpatialIndex.geohash_adjacent("dr5r", "s")
        'dr5n'
        """
        _base32 = "0123456789bcdefghjkmnpqrstuvwxyz"
        _neighbor = {
            "n": ["p0r21436x8zb9dcf5h7kjnmqesgutwvy", "bc01fg45238967deuvhjyznpkmstqrwx"],
            "s": ["14365h7k9dcfesgujnmqp0r2twvyx8zb", "238967debc01fg45kmstqrwxuvhjyznp"],
            "e": ["bc01fg45238967deuvhjyznpkmstqrwx", "p0r21436x8zb9dcf5h7kjnmqesgutwvy"],
            "w": ["238967debc01fg45kmstqrwxuvhjyznp", "14365h7k9dcfesgujnmqp0r2twvyx8zb"],
        }
        _border = {
            "n": ["prxz", "bcfguvyz"],
            "s": ["028b", "0145hjnp"],
            "e": ["bcfguvyz", "prxz"],
            "w": ["0145hjnp", "028b"],
        }

        if direction not in _neighbor:
            raise ValueError("Direction must be one of: 'n', 's', 'e', 'w'.")

        geohash = geohash.lower()
        last = geohash[-1]
        parent = geohash[:-1]
        type_ = len(geohash) % 2

        if last in _border[direction][type_] and parent:
            parent = SpatialIndex.geohash_adjacent(parent, direction)

        neighbor_index = _neighbor[direction][type_].index(last)
        return parent + _base32[neighbor_index]

    @staticmethod
    def geohash_neighbors(geohash: str) -> list[str]:
        """
        Return the 8 neighboring geohashes surrounding a given geohash.

        The neighbors are ordered as follows:
        [northwest, north, northeast, west, center, east, southwest, south, southeast]

        Parameters
        ----------
        geohash : str
            A valid geohash string.

        Returns
        -------
        neighbors : list of str
            A list of 9 geohash strings, including the input geohash and its 8 neighbors.

        Examples
        --------
        >>> SpatialIndex.geohash_neighbors("u4pruydqqvj")
        ['u4pruydqqvh', 'u4pruydqqvm', 'u4pruydqqvn',
         'u4pruydqqvf', 'u4pruydqqvj', 'u4pruydqqvk',
         'u4pruydqqve', 'u4pruydqqvg', 'u4pruydqqvl']

        Notes
        -----
        This function assumes `geohash_adjacent` is implemented as a static method
        of the `SpatialIndex` class.
        """
        n = SpatialIndex.geohash_adjacent(geohash, "n")
        s = SpatialIndex.geohash_adjacent(geohash, "s")
        e = SpatialIndex.geohash_adjacent(geohash, "e")
        w = SpatialIndex.geohash_adjacent(geohash, "w")
        return [
            SpatialIndex.geohash_adjacent(n, "w"),  # NW
            n,  # N
            SpatialIndex.geohash_adjacent(n, "e"),  # NE
            w,  # W
            geohash,  # Center
            e,  # E
            SpatialIndex.geohash_adjacent(s, "w"),  # SW
            s,  # S
            SpatialIndex.geohash_adjacent(s, "e"),  # SE
        ]

    @staticmethod
    def poly_to_geohashes(poly: Polygon, precision: int, inner: bool = False) -> set[str]:
        """
        Cover a polygon with geohashes using pygeohash.

        Parameters
        ----------
        poly : shapely.geometry.Polygon
            Input polygon.
        precision : int
            Geohash precision (length), higher = finer.
        inner : bool
            If True, only keep geohashes fully inside the polygon.

        Returns
        -------
        Set[str]
            Set of geohash strings covering the polygon.
        """
        visited = set()
        result = set()
        q = deque()

        start = pygeohash.encode(poly.centroid.y, poly.centroid.x, precision)
        q.append(start)

        prepared = prep(poly)

        while q:
            gh = q.popleft()
            if gh in visited:
                continue
            visited.add(gh)

            cell = SpatialIndex.geohash_to_poly(gh)

            if inner:
                if prepared.contains(cell):
                    result.add(gh)
                else:
                    continue
            else:
                if prepared.intersects(cell):
                    result.add(gh)
                else:
                    continue

            for ngh in SpatialIndex.geohash_neighbors(gh):
                if ngh not in visited:
                    q.append(ngh)

        return result

    @staticmethod
    def point_cell(lats: list[float], lons: list[float], cell_type: str, res: int) -> list:
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
        >>> point_cell(lats, lons, "geohash", 6)
        ['9q8yy', '9qh0b']
        """
        if cell_type == "geohash":
            return [
                pygeohash.encode(lat, lon, res) if lat is not math.isnan(lat) and not math.isnan(lon) else None
                for lat, lon in zip(lats, lons)
            ]
        elif cell_type == "s2":  # string
            return [
                s2.geo_to_s2(lat, lon, res) if lat is not math.isnan(lat) and not math.isnan(lon) else None
                for lat, lon in zip(lats, lons)
            ]
        elif cell_type == "s2_int":  # int data type requires less memory
            return [
                int(s2.geo_to_s2(lat, lon, res), 16) if lat is not math.isnan(lat) and not math.isnan(lon) else None
                for lat, lon in zip(lats, lons)
            ]
        elif cell_type == "h3":
            return [
                h3.latlng_to_cell(lat, lon, res) if lat is not math.isnan(lat) and not math.isnan(lon) else None
                for lat, lon in zip(lats, lons)
            ]
        else:
            raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', 's2_int' or 'h3'.")

    @staticmethod
    def poly_cell(
        geoms: list[Polygon | MultiPolygon],
        cell_type: str,
        res: int,
        force_full_cover: bool = True,
        dump_path: str = None,
    ) -> list[str] | None:
        """
        Converts a list of geometries into a set of unique spatial cells based on the specified cell type and resolution.

        This function takes a list of Shapely geometries (e.g., Polygon, MultiPolygon) and converts them into spatial cells
        using one of the supported cell systems: Geohash, S2, or H3. The resulting cells are returned as a list of unique
        cell IDs. If `dump_path` is set to a valid directory path, the cells are saved to a file in that directory, instead of being returned.

        Parameters
        ----------
        geoms : list of shapely.geometry.Polygon or shapely.geometry.MultiPolygon
            A list of Shapely geometry objects (Polygon or MultiPolygon).

        cell_type : str
            The type of spatial cell system to use. Supported values are "geohash", "s2", or "h3".

        res : int
            The resolution level for the spatial cells. The resolution parameter determines the granularity of the cells.

        force_full_cover : bool, default=True
            If True, ensures that the entire polygon is fully covered by H3 cells, possibly including cells that extend beyond the polygon boundary.

        dump_path : str, optional
            If set to a valid directory path (string), the cells are saved to a file in the specified folder.
            The file will be saved in a subdirectory structure following the pattern: `/path/to/dir/cell_type/res/`.
            If `dump_path` is None, the function returns the list of cell IDs. Default is None.

        Returns
        -------
        list of str or None
            If `dump_path` is None, a list of unique cell IDs is returned.
            If `dump_path` is provided, None is returned after saving the cells to a file.

        Raises
        ------
        ValueError
            If `cell_type` is not one of the supported values ("geohash", "s2", "h3").

        Examples
        --------
        >>> from shapely.geometry import Polygon, MultiPolygon
        >>> geometries = [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]), MultiPolygon([...])]
        >>> # Convert geometries to H3 cells at resolution 9
        >>> h3_cells = poly_cell(geometries, cell_type="h3", res=9)

        >>> # Convert geometries to S2 cells and save to a directory
        >>> poly_cell(geometries, cell_type="s2", res=10, dump_path="~/Desktop/s2")
        """
        if cell_type == "geohash":
            cell_ids = list(
                {geohash for geom in geoms for geohash in SpatialIndex.poly_to_geohashes(geom, precision=res, inner=False)}
            )
        elif cell_type == "s2":
            # Append GeoJSON-like representations of Shapely geometries (Polygon or MultiPolygon).
            # getattr(geom, 'geoms', [geom]) checks if the geometry has a .geoms attribute (MultiPolygon case)
            # If yes: returns geom.geoms (iterable of Polygons)
            # If no: wraps single Polygon in a list [geom]
            geojson_polys = [g.__geo_interface__ for geom in geoms for g in getattr(geom, "geoms", [geom])]
            cell_ids = list(
                {
                    cell["id"]
                    for geojson_poly in geojson_polys
                    for cell in s2.polyfill(geojson_poly, res, geo_json_conformant=True, with_id=True)
                }
            )
        elif cell_type == "h3":
            if not force_full_cover:
                geojson_polys = [g.__geo_interface__ for geom in geoms for g in getattr(geom, "geoms", [geom])]
                cell_ids = [
                    cell_id for geojson_poly in geojson_polys for cell_id in h3.geo_to_cells(geojson_poly, res)
                ]  # array of h3 IDs covering geoms
            else:  # use the buffered bbox of geom instead.
                mpoly = unary_union(geoms)
                geojson_polys = []
                # edge_length = 1290 * 0.376**resolution  Approximate H3 edge length (in km) using exponential decay model
                # diameter_km = 2 * edge_length           Diameter of a regular hexagon is approximately twice the edge length
                # diameter_deg = diameter_km / 100        Convert diameter from kilometers to degrees (approx. 1 degree ≈ 100 km)
                # Use diameter_deg to buffer the geometry (roughly the width of an H3 hexagon at this resolution)
                buffer_radius = 26 * 0.376**res
                for geom in geoms:
                    minx, miny, maxx, maxy = geom.bounds
                    geom = box(minx, miny, maxx, maxy).buffer(buffer_radius)
                    geojson_polys += [
                        g.__geo_interface__ for g in (getattr(geom, "geoms", [geom]))
                    ]  # array of GeoJson representation of geoms
                cell_ids = [cell_id for geojson_poly in geojson_polys for cell_id in h3.geo_to_cells(geojson_poly, res)]
                cell_ids = list(set(cell_ids))

                hex_polys, _ = SpatialIndex.cell_poly(
                    cell_ids, cell_type="h3"
                )  # Convert H3 cell IDs to Shapely Polygons for Intersection Checks
                inter_hex_polys = [
                    hex_poly for hex_poly in hex_polys if hex_poly.intersects(mpoly)
                ]  # Find the hexagons which intersects with mpoly
                cell_ids = [h3.latlng_to_cell(hex_poly.centroid.y, hex_poly.centroid.x, res) for hex_poly in inter_hex_polys]
        else:
            raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', or 'h3'.")

        if not dump_path:
            return cell_ids
        else:
            # Create the directories if they don't exist
            cells_path = os.path.expanduser(f"{dump_path}/{cell_type}/{res}")
            os.makedirs(cells_path, exist_ok=True)
            with open(f"{cells_path}/{datetime.now()}.txt", "w") as json_file:
                json.dump(cell_ids, json_file)
            return None

    @staticmethod
    def ppoint_cell(lats: list[float], lons: list[float], cell_type: str, res: int) -> list:
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
        >>> ppoint_cell(lats, lons, cell_type, res)
        ['8928308280fffff', '8a28308280fffff']
        """
        n_cores = cpu_count()

        # Prepare arguments for parallel processing
        lat_chunks = np.array_split(lats, 4 * n_cores)
        lon_chunks = np.array_split(lons, 4 * n_cores)
        args = zip(lat_chunks, lon_chunks, [cell_type] * 4 * n_cores, [res] * 4 * n_cores)

        # Parallelize the conversion using Pool.starmap
        with Pool(n_cores) as pool:
            cells = pool.starmap(SpatialIndex.point_cell, args)
        cells = [item for sublist in cells for item in sublist]  # Flatten the list of cells

        return cells

    @staticmethod
    def ppoly_cell(
        mdf: gpd.GeoDataFrame,
        cell_type: str,
        res: int,
        force_full_cover: bool = True,
        compact: bool = False,
        dump_path: str = None,
        verbose: bool = False,
    ) -> tuple[list[str], int]:
        """
        Performs a parallelised conversion of geometries in a GeoDataFrame to cell identifiers of a specified type
        (e.g., Geohash, S2, or H3), optionally compacting the result to reduce the number of cells.

        This function first divides the bounding box of the input GeoDataFrame into smaller grid cells, then calculates
        the intersection between these grid cells and the input geometries. The resulting geometries are processed in
        parallel to generate cell identifiers according to the specified `cell_type` and `res` (resolution). The result
        can be compacted to reduce the number of cells. Optionally, if `dump_path` is provided, the results are saved in multiple
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

        force_full_cover : bool, default=True
            If True, ensures that the entire polygon is fully covered by H3 cells, possibly including cells that extend beyond the polygon boundary.

        compact : bool, optional, default=False
            If True, compact the resulting cells to reduce their number. This is typically applicable for S2 and H3 cells.

        dump_path : str, optional
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
        >>> cells, count = ppoly_cell(mdf, cell_type="h3", res=10, compact=True, dump_path="~/Desktop/h3", verbose=True)
        >>> print(f"Generated {count} cells: {cells}")
        """
        # Determine the number of slices and grid cells based on CPU cores
        n_cores = cpu_count()
        slices = 128 * n_cores

        if dump_path and not verbose:
            cells_path = os.path.abspath(os.path.expanduser(f"{dump_path}/{cell_type}/{res}"))
            print(f"Writing cell IDs in parallel to {4 * n_cores} files in {cells_path} directory ...")

        if verbose:
            print(datetime.now())
            print("\nSlicing the bounding box of polygons ...")
            start_time = time()

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
            print(f"    {slices:,} slices created")
            print(f"    {elapsed_time} seconds")

            print("Performing intersection between grid and polygons ...")
            start_time = time()

        # Perform intersection between input geometries and grid cells
        gmdf = gpd.overlay(mdf, gmdf, how="intersection")  # grid mdf

        if verbose:
            elapsed_time = round(time() - start_time)
            print(f"    {len(gmdf):,} intersected slices")
            print(f"    {elapsed_time} seconds")

            if dump_path:
                cells_path = os.path.abspath(os.path.expanduser(f"{dump_path}/{cell_type}/{res}"))
                print(f"Writing cell IDs in parallel to {4 * n_cores} files in {cells_path} directory ...")
            else:
                print("Calculating cell IDs in parallel ...")
            start_time = time()

        # Shuffle geometries for even load distribution across chunks
        gmdf = gmdf.sample(frac=1)
        geom_chunks = np.array_split(list(gmdf.geometry), 4 * n_cores)
        args = zip(
            geom_chunks,
            [cell_type] * 4 * n_cores,
            [res] * 4 * n_cores,
            [force_full_cover] * 4 * n_cores,
            [dump_path] * 4 * n_cores,
        )

        # Parallel processing to generate cells
        if dump_path:
            with Pool(n_cores) as pool:
                pool.starmap(SpatialIndex.poly_cell, args)
            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"    {elapsed_time} seconds")
            return
        else:
            with Pool(n_cores) as pool:
                cells = pool.starmap(SpatialIndex.poly_cell, args)
            cells = [item for sublist in cells for item in sublist]  # Flatten the list of cells

            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"    {len(cells):,} cells")
                print(f"    {elapsed_time} seconds")

            # Remove duplicates based on cell type
            if verbose:
                print("Removing duplicate cells ...")
                start_time = time()
            cells = list(set(cells))  # Remove duplicate cells
            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"    {len(cells):,} cells")
                print(f"    {elapsed_time} seconds")

            cell_counts = len(cells)  # Total unique cell count

            # Compact the cells if needed
            if compact:
                if verbose:
                    print("Compacting cells ...")
                    start_time = time()
                cells = CellUtils.compact_cells(cells, cell_type)
                if verbose:
                    elapsed_time = round(time() - start_time)
                    print(f"    {len(cells):,} cells")
                    print(f"    {elapsed_time} seconds")

            return cells, cell_counts

    @staticmethod
    def cell_point(cells: list[str | int], cell_type: str) -> list[tuple[float, float]]:
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
        >>> cell_point(["ezs42", "u4pruydqqvj"], cell_type="geohash")
        [(42.6, -5.6), (57.64911, 10.40744)]

        >>> cell_point(["8928308280fffff"], cell_type="h3")
        [(37.775938728915946, -122.41795063018799)]

        >>> cell_point([9744573459660040192], cell_type="s2_int")
        [(37.7749, -122.4194)]

        >>> cell_point(["89c25c"], cell_type="s2")
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
    def cell_poly(cells: list, cell_type: str) -> tuple:
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
            - `geoms` : list of shapely.geometry.Polygon
                A list of Polygon geometries representing the spatial boundaries of the input cells.
            - `res` : list of int
                A list of resolution levels corresponding to each cell in the input.

        Raises
        ------
        ValueError
            If `cell_type` is not one of "geohash", "h3", or "s2".

        Example
        -------
        >>> from shapely.geometry import Polygon
        >>> cells = ["ezs42", "ezs43"]  # Geohash cells
        >>> cell_type = "geohash"
        >>> geoms, res = cell_poly(cells, cell_type)
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
        if cell_type not in {"geohash", "s2", "s2_int", "h3"}:
            raise ValueError(f"Invalid cell_type '{cell_type}'. Accepted values are: 'geohash', 'h3', 's2' and 's2_int'.")

        if cell_type == "s2_int":
            cells = [hex(cell)[2:] if cell is not None else None for cell in cells]

        # Create geometry objects based on cell type
        geoms = [
            None
            if cell is None
            else SpatialIndex.geohash_to_poly(cell)
            if cell_type == "geohash"
            # Shapely expects (lng, lat) format, so we reverse the coordinates returned by cell_to_boundary
            else Polygon([(lng, lat) for lat, lng in h3.cell_to_boundary(cell)])
            if cell_type == "h3"
            else Polygon(s2.s2_to_geo_boundary(cell, geo_json_conformant=True))
            for cell in cells
        ]

        # Determine resolution level based on cell type
        res = [
            np.nan
            if cell is None
            else len(cell)
            if cell_type == "geohash"
            else int(cell[1], 16)
            if cell_type == "h3"
            else s2.CellId.from_token(cell).level()
            for cell in cells
        ]

        return geoms, res

    @staticmethod
    def pcell_point(cells: list[str | int], cell_type: str) -> list[tuple[float, float]]:
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
            points = pool.starmap(SpatialIndex.cell_point, args)
        points = [item for sublist in points for item in sublist]  # Flatten the list of cells

        return points

    @staticmethod
    def pcell_poly(cells: list[str | int], cell_type: str) -> tuple:
        """
        Parallelized version of `cell_poly`, converting a list of spatial cells to geometries and resolution levels.

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
            - `geoms` : list of shapely.geometry.Polygon
                Polygon geometries representing the boundaries of input cells.
            - `res` : list of int
                Resolution levels for each cell in the input.
        """
        n_cores = cpu_count()

        cell_chunks = np.array_split(cells, 4 * n_cores)
        # Convert each numpy array to a list which converts numpy.str_ to str
        cell_chunks = [arr.tolist() for arr in cell_chunks]
        args = zip(cell_chunks, [cell_type] * 4 * n_cores)

        with Pool(n_cores) as pool:
            results = pool.starmap(SpatialIndex.cell_poly, args)

        # Unpack `res` and `geoms` from the result tuples
        geoms = [g for result in results for g in result[0]]
        res = [r for result in results for r in result[1]]

        return geoms, res


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

    geom_in_poly(gdf: gpd.GeoDataFrame, mdf: gpd.GeoDataFrame, mdf_col: Optional[str] = None) -> Union[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, np.ndarray]]
        Performs a spatial intersection between two GeoDataFrames and returns the intersecting subset.

    poverlay(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame, how: str = "intersection", keep_geom_type: bool = False) -> gpd.GeoDataFrame
        Executes a parallelized spatial overlay operation between two GeoDataFrames.

    Notes
    -----
    - The class relies on GeoPandas and Shapely for spatial operations and multiprocessing for parallelization.
    - Ensure that input GeoDataFrames have the same coordinate reference system (CRS) for accurate results.

    Examples
    --------
    >>> # Example for flatten_3d
    >>> gdf_2d = SpatialOps.flatten_3d(gdf.geometry)
    >>> print(gdf_2d)

    >>> # Example for line_to_points
    >>> point_gdf = SpatialOps.line_to_points(line_gdf.iloc[0])
    >>> print(point_gdf)

    >>> # Example for geom_in_poly
    >>> gdf_mask, geom_counts, poly_ids = SpatialOps.geom_in_poly(gdf, mdf, mdf_col="region_id")
    >>>

    >>> # Example for poverlay
    >>> result_gdf = SpatialOps.poverlay(gdf1, gdf2, how="intersection")
    >>> print(result_gdf)
    """

    @staticmethod
    def geom_in_poly(
        gdf: gpd.GeoDataFrame, mdf: gpd.GeoDataFrame, mdf_col: str | None = None
    ) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute spatial intersections between two GeoDataFrames.

        For each geometry in `gdf`, determine if it intersects with any polygon in `mdf` (which must
        contain MultiPolygons or Polygons). Optionally track which `mdf` polygon contains each
        intersecting `gdf` geometry.

        Parameters
        ----------
        gdf : gpd.GeoDataFrame
            GeoDataFrame containing geometries to test for intersection with `mdf`.
            Can be any geometry type (Points, LineStrings, Polygons, etc.).
        mdf : gpd.GeoDataFrame
            GeoDataFrame containing MultiPolygon/Polygon geometries to use as spatial filter.
            All geometries must be polygonal (will raise error if lines/points are found).
        mdf_col : str, optional
            Column name in `mdf` containing unique polygon identifiers.
            Required if you need to track which polygon contains each intersecting geometry.

        Returns
        -------
        gdf_mask : np.ndarray[bool]
            Boolean array where True indicates the corresponding `gdf` geometry
            intersects at least one polygon in `mdf`.
            Shape: (len(gdf),). Use as: `filtered_gdf = gdf[gdf_mask]`
        geom_counts : np.ndarray[int]
            Array counting how many `gdf` geometries intersect each `mdf` polygon.
            Shape: (len(mdf),). Use as: `mdf['count'] = geom_counts`
        poly_ids : np.ndarray[object], optional
            Only returned when `mdf_col` is specified.
            Contains the `mdf_col` value of the containing polygon for each intersecting
            `gdf` geometry (None where no intersection exists).
            Shape: (len(gdf),). Use as: `gdf['poly_id'] = poly_ids`

        Notes
        -----
        1. This function is faster than `gpd.sjoin` only for small datasets (typically fewer than 2,000 features).
           For larger datasets, `gpd.sjoin` is generally more efficient and scalable.
        2. For `gdf` geometries intersecting multiple `mdf` polygons, only the last
           encountered polygon ID (from iteration order) is kept in `poly_ids`.

        Examples
        --------
        >>> import geopandas as gpd

        # Basic intersection check
        >>> points = gpd.read_file("points.shp")
        >>> zones = gpd.read_file("zones.shp")
        >>> gdf_mask, geom_counts = geom_in_poly(points, zones)

        # With polygon ID tracking
        >>> gdf_mask, geom_counts, poly_ids = geom_in_poly(points, zones, mdf_col="zone_id")
        >>> # Get only points within zones
        >>> points_in_zones = points[gdf_mask]
        >>> # Add zone IDs to the points
        >>> points["zone_id"] = poly_ids
        """
        gdf_mask = np.zeros(len(gdf), dtype=bool)
        geom_counts = np.zeros(len(mdf), dtype=int)

        if mdf_col is not None:
            poly_ids = np.full(len(gdf), None, dtype=object)

        for idx, geom in enumerate(mdf.geometry):
            mask = gdf.intersects(geom)
            gdf_mask |= mask
            geom_counts[idx] = mask.sum()

            if mdf_col is not None and mask.any():
                poly_ids[mask] = mdf.at[idx, mdf_col]

        return (gdf_mask, geom_counts, poly_ids) if mdf_col is not None else (gdf_mask, geom_counts)

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
    def proximity_counts(coords: np.ndarray, coords_crs: int | str | pyproj.CRS = 4326, radius: int = 100) -> list[int]:
        """
        Calculate the number of neighboring points within a specified radius for each point.

        Transforms coordinates to British National Grid (EPSG:27700), builds a KD-Tree,
        and performs radius queries to count neighbors within the given radius.

        Parameters
        ----------
        coords : np.ndarray
            Array of point coordinates with shape (n_points, 2).
            Coordinates should be in longitude/latitude order (x, y) if using WGS84.
            Expected dtype: float (np.float32 or np.float64).

        coords_crs : Union[int, str, pyproj.CRS], optional
            Coordinate Reference System of input coordinates.
            Default: 4326 (WGS84 - longitude/latitude).
            Can be EPSG code (int or str), PROJ string, or pyproj.CRS object.

        radius : float, optional
            Search radius in meters. All points within this distance will be counted
            as neighbors (excluding the point itself). Default: 100.0 meters.

        Returns
        -------
        List[int]
            List of neighbor counts for each input point, with length equal to n_points.
            Each value represents the number of other points within the specified radius
            (self-count is excluded).

        Raises
        ------
        ValueError
            If `coords` is empty or has incorrect shape (not (n, 2)).
            If `radius` is negative or zero.

        TypeError
            If `coords` is not a numpy array.

        pyproj.exceptions.CRSError
            If `coords_crs` cannot be parsed to a valid coordinate reference system.

        Notes
        -----
        1. The function transforms coordinates to EPSG:27700 (British National Grid)
           before distance calculations, ensuring metric distances are accurate.
        2. Self-counts are excluded from neighbor counts.
        3. The KD-Tree query uses Euclidean distance in the projected CRS (EPSG:27700),
           which is accurate for distances up to tens of kilometers in the UK.
        4. For large datasets (>50,000 points), consider using approximate nearest
           neighbor methods or spatial indexing for better performance.

        Examples
        --------
        >>> import numpy as np
        >>> coords = np.array(
        ...     [
        ...         [-0.1278, 51.5074],  # London
        ...         [-1.8904, 52.4862],  # Birmingham
        ...         [-2.2426, 53.4808],  # Manchester
        ...     ]
        ... )

        >>> # Count neighbors within 50km
        >>> counts = proximity_counts(coords, radius=50000)
        >>> counts
        [0, 0, 0]  # No points within 50km of each other

        >>> # Dense cluster of points
        >>> cluster = np.array(
        ...     [
        ...         [530000, 180000],  # BNG coordinates (London area)
        ...         [530100, 180100],
        ...         [530050, 180050],
        ...     ]
        ... )
        >>> counts = proximity_counts(cluster, coords_crs=27700, radius=150)
        >>> counts
        [2, 1, 2]  # Each point has 1-2 neighbors within 150 meters

        >>> # Edge case: radius too small
        >>> counts = proximity_counts(coords, radius=0.1)
        >>> counts
        [0, 0, 0]  # No points within 10cm

        See Also
        --------
        sklearn.neighbors.KDTree : KD-Tree implementation used for spatial queries.
        pyproj.Transformer : Coordinate transformation utility.
        scipy.spatial.cKDTree : Alternative KD-Tree implementation (faster for pure queries).

        Warnings
        --------
        - Using geographic coordinates (lon/lat) with metric radius without projection
          will give incorrect results. This function handles projection automatically.
        - For global datasets, consider using Great Circle distance or Haversine formula
          instead of projecting to a local CRS like EPSG:27700.
        - Memory usage scales with O(n²) for dense clusters with large radii.
        """
        # Validate inputs
        if not isinstance(coords, np.ndarray):
            raise TypeError(f"coords must be numpy array, got {type(coords)}")

        if coords.ndim != 2 or coords.shape[1] != 2:
            raise ValueError(f"coords must have shape (n, 2), got {coords.shape}")

        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")

        if len(coords) == 0:
            return []

        # Transform coordinates from WGS84 to British National Grid (27700)
        transformer = pyproj.Transformer.from_crs(coords_crs, 27700, always_xy=True).transform
        coords = np.array(transformer(coords[:, 0], coords[:, 1])).T
        # Build KDTree and query neighbors
        tree = KDTree(coords)
        indices = tree.query_ball_tree(tree, r=radius)
        # Count neighbors excluding self
        neighbor_counts = [len(nbrs) - 1 for nbrs in indices]
        return neighbor_counts

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


plp = Karta.plp
