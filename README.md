<h1 align="center"><img align="center" src="https://geoparse.io/graphics/geoparse_logo.png" alt="GeoParse Logo" width="200"/></h1>
<h1 align="center">GeoParse</h1>
<h3 align="center">All About Points <img src="https://geoparse.io/graphics/point.png" width="10"/> Lines <img src="https://geoparse.io/graphics/line.png" width="40"/> and Polygons <img src="https://geoparse.io/graphics/polygon.png" width="30"/></h3>

---

[![GeoParse](https://img.shields.io/badge/GeoParse-008000.svg)](https://geoparse.io)
[![Tutorials](https://img.shields.io/badge/Tutorials-Start%20Learning-orange?logo=book&logoColor=white)](https://geoparse.io/tutorials)
[![Docs](https://img.shields.io/badge/Docs-Read%20Now-blue?logo=readthedocs&logoColor=white)](https://geoparse.io)
[![License](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative&logoColor=white)](https://github.com/geoparse/geospatial/blob/main/LICENSE)
[![PythonVersion](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue?style=flat&logo=python&logoColor=white&labelColor=gray)](https://www.python.org/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-Ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Contributors](https://img.shields.io/github/contributors/geoparse/geospatial)](https://github.com/geoparse/geospatial/graphs/contributors)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

GeoParse is an open-source project that provides tools for the visualization, analysis, and manipulation of vector geospatial data.
It supports a wide range of applications, including
telematics analysis (e.g., vehicle trajectories and movement patterns),
infrastructure mapping (e.g., roads, rivers, water streams, buildings, and lakes),
and
footfall analytics (e.g., pedestrian density and mobility trends).

GeoParse builds on top of popular tools and libraries like 
[GDAL](https://gdal.org/en/stable/), 
[GeoPandas](https://geopandas.org/en/stable/), 
[Folium](https://python-visualization.github.io/folium/latest/), 
and [Overpass API](https://wiki.openstreetmap.org/wiki/Overpass_API) 
providing a powerful toolkit for working with geospatial data.
It focuses on straightforward visualization, efficient spatial indexing, advanced geometry manipulations, and user-friendly utilities for handling OpenStreetMap data effortlessly.

With GeoParse, you can:

1. **Create Interactive Maps** – Generate interactive maps with multiple tile layers (e.g., OpenStreetMap, Satellite, Dark Mode).
2. **Visualize Multiple Formats** – Display data from CSV and GIS formats like Shapefile, GPKG, GeoJSON, and GeoParquet.
3. **Generate Heatmaps and Clusters** – Create heatmaps and clusters from point data to visualize density and patterns.
4. **Build Choropleth Maps** – Create thematic maps to visualize data distributions across regions.
5. **Use Spatial Indexing** – Convert geometries into spatial grids like [Geohash](https://en.wikipedia.org/wiki/Geohash), [S2](https://github.com/google/s2geometry), and [H3](https://github.com/uber/h3) for efficient spatial queries, analysis and grid visualization.
6. **Run Parallelized Spatial Operations** – Perform spatial operations in parallel to handle large datasets efficiently.
7. **Leverage Geospatial Utilities** – Perform advanced spatial operations such as intersection, union, and distance calculations (e.g., Haversine, Vincenty).
8. **Integrate with OpenStreetMap (OSM)** – Fetch and convert [OpenStreetMap (OSM)](https://www.openstreetmap.org/about) way IDs into Shapely geometries for visualization and analysis.
9. **Geocode and Reverse Geocode** – Use the Google Maps API to convert addresses to coordinates and vice versa.
10. **Match Maps** – Align GPS points to road networks using the [Valhalla](https://github.com/valhalla/valhalla) routing engine.

# Ducumentation
The documentation is available online at [geoparse.io](https://geoparse.io/).

# Install
```sh
pip install git+https://github.com/geoparse/geoparse.git
```

---

# Tutorials

GeoParse includes the following classes for visualizing and analyzing geospatial data.

1. [Karta Class](https://geoparse.io/tutorials/karta.html) named after the Swedish word for *map*, is designed for creating and customizing interactive maps.

3. **GeomUtils Class** stands for *Geometry Utilities* and provides utility functions for working with geometry objects. Key functionalities include determining UTM projections, transforming geometries between coordinate reference systems (CRS), and calculating geometric statistics like area and perimeter.

4. **CellUtils Class** stands for *Cell Utilities* and provides utility functions for compacting and uncompacting spatial cells like H3, S2, and Geohash. It also supports statistical analysis of spatial cells, such as calculating the number of cells and their area for a given geometry.

5. **OSMUtils Class** stands for *OSM Utilities* and provides utility functions for working with OpenStreetMap (OSM) data and routing engines built on top of OSM. It allows users to retrieve OSM way geometries (either polygons or lines) and perform map matching of GPS coordinates to road networks. 

6. **SpatialIndex Class** provides methods for *spatial indexing* operations. It converts geographic coordinates and geometries into spatial index representations and vice versa, utilizing popular encoding systems like H3, S2, and Geohash. It also accelerates computations through parallel processing, making it useful for efficient spatial queries and handling large datasets.

7. **SpatialOps Class** stands for *Spatial Operations* and provides methods for handling 3D geometries, converting LineStrings to Points, performing spatial intersections, and executing parallelized spatial overlay operations. It also includes utilities for calculating distances (e.g., haversine and Vincenty formulas) and geocoding addresses using the Google Geocoding API.

## Geospatial Data Visualization

Karta class provides methods for adding points, lines, polygons, and choropleth layers to a map. The primary function, `plp` (point, line, polygon), supports various visualization styles and configurations, including clustering, heatmaps, and cell-based layers (e.g., H3, S2, geohash). `plp` can visualize the vector data on a map with the following tile layers.
| Light | Dark |
| ----- | ---- | 
| A minimalist, light-colored basemap that serves as a subtle background, emphasizing overlaid data. | A high-contrast, dark-themed map ideal for vibrant data overlays and nighttime aesthetics. |
| <img src="https://geoparse.io/graphics/layer_light_london.png?" height="400">          | <img src="https://geoparse.io/graphics/layer_dark_london.png?" height="400">   |

| Outdoors | Satellite |
| -------- | --------- |
| Designed for outdoor enthusiasts, featuring hiking trails, biking paths, natural landmarks, and elevation contours. | A basemap displaying satellite imagery of the Earth's surface, useful for real-world context and analyses requiring detailed imagery. |
| <img src="https://geoparse.io/graphics/layer_outdoors_kiasar.png?" height="400"> | <img src="https://geoparse.io/graphics/layer_satellite_venice.png?" height="400"> |

| OSM |
| --- |
| A general-purpose map powered by OpenStreetMap, showcasing roads, buildings, and points of interest (POIs).   |
| <div align="center"><img src="https://geoparse.io/graphics/layer_osm_sf.png?" height="400"></div> |

`plp` accepts either a pandas `DataFrame` or a GeoPandas `GeoDataFrame` to render geometry data. For a `DataFrame`, the `plp` method in `Karta` class automatically identifies columns with names containing "lat" and "lon" (case-insensitive) to use as latitude and longitude for plotting points on the map. If no columns contain these keywords, or if more than two columns contain these keywords, you must explicitly specify the latitude and longitude using the `y` and `x` parameters, respectively, e.g., `plp(df, x="easting", y="northing")`. Note that plp assumes all data is in the [EPSG:4326](https://epsg.io/4326) projection. For a `GeoDataFrame`, the `plp` function can render Shapely objects such as `Point`, `LineString`, `Polygon`, and `MultiPolygon`.

### Point
In the following example, we demonstrate how to display points from a CSV file, customize the map with point colors and popups, and add layers such as heatmaps and clusters.

```python
import pandas as pd
from geoparse.geoparse import plp

df = pd.read_csv("https://geoparse.io/tutorials/data/fatal_crash_great_britain_2023.csv")
df.head()
```
| date       | time  | latitude  | longitude  | number_of_casualties | speed_limit | highway  |
|------------|-------|-----------|------------|----------------------|-------------|----------|
| 03/01/2023 | 19:12 | 51.356551 | -0.097759  | 1                    | 30          | primary  |
| 07/01/2023 | 10:05 | 51.593701 | 0.022379   | 1                    | 30          | primary  |
| 14/01/2023 | 16:15 | 51.466689 | -0.011289  | 1                    | 20          | primary  |
| 15/01/2023 | 19:51 | 51.671577 | -0.037543  | 1                    | 30          | tertiary |
| 16/01/2023 | 19:22 | 51.447944 | 0.117279   | 1                    | 30          | primary  |

After loading the data, you can easily visualize it on a map using `plp(df)`.
By default, the `plp` method displays points in blue (left figure), but you can change the point color using the `point_color` argument (center figure).
To apply custom colors, you can also use RGB hex codes (right figure).

<table width="100%">
<table>
  <tr>
    <td style="vertical-align: top;"><pre><code>plp(df)                  </code></pre></td>
    <td style="vertical-align: top;">
      <pre><code>plp(
    df,
    point_color="purple"   
)</code></pre>
    </td>
    <td style="vertical-align: top;">
      <pre><code>plp(
    df,
    point_color="#cc5500" 
)</code></pre>
    </td>
  </tr>
  <tr>
    <td width="33%"><img src="https://geoparse.io/graphics/casualty_map.png?"/></td>
    <td width="33%"><img src="https://geoparse.io/graphics/casualty_map_purple.png?"/></td>
    <td width="33%"><img src="https://geoparse.io/graphics/casualty_map_brown.png?"/></td>
  </tr>
</table>


To enhance interpretability, `plp` allows grouping and coloring points based on the values of a specific feature such as `number_of_casualties`. 
In the visualizations below, crash points are color-coded: blue for 1 casualty, green for 2, and red for 3. 
Additionally, by using the `point_popup` argument, hovering over a point reveals a popup with detailed contextual information, 
as illustrated in the image on the right.

<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(
    df,
    point_color="number_of_casualties"
)
      </code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(
    df,
    point_color="number_of_casualties",
    point_popup={
        "Date": "date",
        "Casualties": "number_of_casualties"
    }
)
      </code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <div align="center"><img src="https://geoparse.io/graphics/casualty_colored.png?" height="400"></div>
    </td>
    <td>
      <div align="center"><img src="https://geoparse.io/graphics/casualty_colored_popup.png?" height="400"></div>
    </td>
  </tr>
</table>

#### Heatmap and cluster
The `plp` function can also add heatmap and cluster layers to the map. In the left image, we see a heatmap of fatal road crashes in Great Britain. The center image displays the cluster view. Both layers can also be shown together on a single map, as demonstrated in the right image.

<table>
  <tr>
    <td style="vertical-align: bottom; text-align: center;">
      <pre><code>
plp(df, heatmap=True)
      </code></pre>
    </td>
    <td style="vertical-align: bottom; text-align: center;">
      <pre><code>
plp(df, cluster=True)
      </code></pre>
    </td>
    <td style="vertical-align: bottom; text-align: center;">
      <pre><code>
plp(
    df, 
    heatmap=True,
    cluster=True
)
      </code></pre>
    </td>
  </tr>
  <tr>
    <td style="text-align: center;">
      <img src="https://geoparse.io/graphics/casualty_heatmap.png" height="360" alt="Default plot">
    </td>
    <td style="text-align: center;">
      <img src="https://geoparse.io/graphics/casualty_cluster.png" height="360" alt="Purple points plot">
    </td>
    <td style="text-align: center;">
      <img src="https://geoparse.io/graphics/casualty_heatmap_cluster.png" height="360" alt="Brown points plot">
    </td>
  </tr>
</table>


#### Buffer and ring
`plp` can create a buffer zone around each point, forming a circular area centered on the point. This is useful for visualizing spatial influence or performing proximity-based analysis. For example, it can help identify features within 100 meters of a crash site, as shown in the left figure. `plp` can also generate a ring-shaped buffer, sometimes called a "donut buffer," around each point. Each ring is defined by an inner and outer radius. In the example shown in the right cell, the ring starts 100 meters from each point and extends to 200 meters. This approach is useful when you want to exclude the immediate area around a point and focus on a specific surrounding zone.



<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(df, buffer_radius=100)
      </code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(
    df, 
    ring_inner_radius=100,
    ring_outer_radius=200
)
      </code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://geoparse.io/graphics/casualty_buffer.png" height="480">
    </td>
    <td>
      <img src="https://geoparse.io/graphics/casualty_ring.png" height="480">
    </td>
  </tr>
</table>

#### Trajectory - Visualization
`plp` processes a trajectory as an ordered set of points and renders it using point-based visualization, as previously described. By default, the points are colored blue, but this can be changed using the `point_color` parameter. As with other point visualizations mentioned earlier, you can customize the color based on a feature value, e.g., `vin` (vehicle identification number) in the right-hand map, to distinguish between different journeys.

```python
df = pd.read_csv("https://geoparse.io/tutorials/data/trajectory.csv", dtype={"speedlimit": "Int8"})
df.head()
```
| vin | lat       | lon       | dt                  | speed_mph | highway      | name           | ref  | speedlimit_mph |
|-----|-----------|-----------|---------------------|-----------|--------------|----------------|------|----------------|
| 13  | 53.420805 | -2.212972 | 2023-11-18 23:46:35 | 0.0       | residential  | Lane End Road  | NaN  | \<NA\>         |
| 13  | 53.420618 | -2.213305 | 2023-11-18 23:47:35 | 19.0      | trunk        | Kingsway       | A34  | 40             |
| 13  | 53.418982 | -2.214270 | 2023-11-18 23:48:35 | 32.6      | trunk        | Kingsway       | A34  | 40             |
| 13  | 53.414234 | -2.217178 | 2023-11-18 23:49:35 | 30.1      | trunk        | Kingsway       | A34  | 40             |
| 13  | 53.407416 | -2.221233 | 2023-11-18 23:50:35 | 37.6      | trunk        | Kingsway       | A34  | 40             |


<table>
  <tr>
    <td width="410" style="vertical-align: top;">
        <pre><code>plp(df)</code></pre>
    </td>
    <td width="447" style="vertical-align: top;">
      <pre><code>plp([df, df2],
    point_color="vin"
)</code></pre>
    </td>
  </tr>
  <tr>
    <td><img src="https://geoparse.io/graphics/traj_points.png?"></td>
    <td><img src="https://geoparse.io/graphics/trajs.png?"></td>
  </tr>
</table>


In addition, specifically for trajectories, `plp` can display lines that connect points to form continuous paths. It also supports animated ant paths for line geometries, allowing users to visualize the direction of movement.

<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>plp(df, line=True)</code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>plp(df, antpath=True)</code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://geoparse.io/graphics/traj_line.png?" height="480">
    </td>
    <td>
      <img src="https://geoparse.io/graphics/traj_antpath.gif?" height="480">
    </td>
  </tr>
</table>

#### Trajectory - Speeding detection

By passing the name of the speed column (`speed_mph` in the above example) to the `point_color` parameter, and specifying both `speed_field` and `speed_limit_field` in the `plp` function, each GPS point is plotted and color-coded based on how the vehicle's speed compares to the road's speed limit. The `speed_field` and `speed_limit_field` parameters indicate which columns contain the actual vehicle speed and the legal speed limit, respectively. To convey speeding severity, a color-coded band is applied according to the percentage by which the vehicle exceeds the speed limit. For instance, a speeding value of 10% corresponds to driving at 55 miles per hour on a road with a 50 mile per hour limit.

    Blue: No speeding
    Green: 0% < Speeding < 10%
    Yellow: 10% ≤ Speeding < 20%
    Orange: 20% ≤ Speeding < 30%
    Red: 30% ≤ Speeding < 40%
    Black: Speeding ≥ 40%
    Purple: Speed limit unavailable
    
Using the `point_popup` argument, when a user hovers over a point on the map, a popup will display relevant information for that location. This feature supports exploratory analysis of driving behavior in relation to road regulations. Combined with the visual encoding of speeding severity, it allows users to quickly identify areas of excessive speeding, compliance, or missing speed limit data.
<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(
    df,
    speed_field="speed_mph",
    speed_limit_field="speedlimit_mph",
    point_color="speed_mph",
)
</code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(
    df,
    speed_field="speed_mph",
    speed_limit_field="speedlimit_mph",
    point_color="speed_mph",
    point_popup={
        "Time": "dt", 
        "Speed (mile/h)": "speed_mph", 
        "Speed Limit (mile/h)": "speedlimit_mph"
    }
)
</code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://geoparse.io/graphics/traj_speeding.png?" height="480">
    </td>
    <td>
      <img src="https://geoparse.io/graphics/traj_speeding_popup.png?" height="480">
    </td>
  </tr>
</table>

It's worth mentioning that the default values for `speed_field` and `speed_limit_field` are `speed` and `speedlimit`, respectively. 
So in this case, there's no need to specify these arguments.

```python
df = df.rename(columns={"speed_mph": "speed", "speedlimit_mph": "speedlimit"})
plp(df, point_color="speed")
```

### Line and Polygon
Using `GeoPandas`, we can read a geospatial file and display its contents using `plp` function. The left image illustrates the border of Luxembourg, represented as a Shapely `Polygon` object. The center image depicts the main roads in Luxembourg, represented as Shapely `LineString` objects. Additionally, `plp` can accept two `GeoDataFrame` objects as a list and display both of them on a single map, as shown in the right image.

| Polygn                                       | LineString                                        | Both                                                            |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf)`                            | `plp(road_gdf, line_color='gold', line_weight=1)` | `plp([border_gdf, road_gdf], line_color='gold', line_weight=1)` |
|![](https://geoparse.io/graphics/luxembourg_border.png) | ![](https://geoparse.io/graphics/luxembourg_roads.png)      | ![](https://geoparse.io/graphics/luxembourg_border_roads.png)             |


Using the `plp` function, we can also add spatial index polygonal layers such as `GeoHash`, Google `S2`, and Uber `H3`. 


| Geohash                                      | S2                                                | H3                                                              |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf, geohash_res=5)`             | `plp(border_gdf, s2_res=11)`                      | `plp(border_gdf, h3_res=6)`                                     |
|![](https://geoparse.io/graphics/geohash_5.png)         | ![](https://geoparse.io/graphics/s2_11.png)                 | ![](https://geoparse.io/graphics/h3_6.png)                                |
 
If the `compact` parameter is set to True, `plp` calculates the parent cell IDs to create a compact representation of the cells.

| Geohash                                      | S2                                                | H3                                                              |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf, geohash_res=7, compact=True)`             | `plp(border_gdf, s2_res=15, compact=True)`                      | `plp(border_gdf, h3_res=10, compact=True)`                                     |
|![](https://geoparse.io/graphics/geohash_compact.png)         | ![](https://geoparse.io/graphics/s2_compact.png)                 | ![](https://geoparse.io/graphics/h3_compact.gif)                                |
 

### OSM Ways

The `plp` function can also accept `OpenStreetMap` (OSM) Way IDs instead of `DataFrame` or `GeoDataFrame` objects and visualize them as `LineString` or `Polygon` geometries. The left image illustrates two ways of the Tokyo Metro Line represented as `LineString` geometries, while the right image depicts three ways in the Louvre Museum visualized as `Polygon` geometries.

|                                         |                                                                            | 
| --------------------------------------- | ---------------------------------------------------------------------------|
| `plp(osm_ways=[893074361, 666201240])`  | `plp(osm_ways=[335265936, 53820456, 1117218957], s2_res=22, compact=True)` | 
|![](https://geoparse.io/graphics/osm_way_line.png) | ![](https://geoparse.io/graphics/osm_way_polygon.gif)                                | 



We recommend starting your GeoParse journey with the HTML versions of the tutorial notebooks, available [here](https://geoparse.io/tutorials/karta.html).
The notebooks cover a wide range of use cases to help you get familiar with GeoParse’s capabilities. 
You can find them on [GitHub](https://github.com/geoparse/geoparse/tree/main/tutorials). 
To run the notebooks, use the following command after setting up the environment described in the 
[contributing](https://github.com/geoparse/geoparse/blob/main/CONTRIBUTING.md#environment-setup)
guide on GitHub:

```sh
jupyter lab 
```

# Contributing to GeoParse
We appreciate all contributions, including bug reports, fixes, documentation improvements, feature enhancements, and ideas. 
A detailed overview on how to contribute can be found in the [contributing](https://github.com/geoparse/geoparse/blob/main/CONTRIBUTING.md) guide on GitHub.

# License and Credits
GeoParse is licensed under the [MIT license](https://github.com/geoparse/geoparse/blob/main/LICENSE) and is written and maintained by Abbas Kiasari (abbas@geoparse.io).
