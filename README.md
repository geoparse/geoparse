
<h1 align="center"><img align="center" src="tutorials/data/geoparse_logo.png" alt="GeoParse Logo" width="200"/></h1>
<h1 align="center">GeoParse</h1>
<h3 align="center">It's all about points <img src="tutorials/data/point.png" width="10"/> lines <img src="tutorials/data/line.png" width="40"/> and polygons <img src="tutorials/data/polygon.png" width="30"/></h3>

---

[![GeoParse](https://img.shields.io/badge/GeoParse-008000.svg)](https://geoparse.io)
[![Docs](https://img.shields.io/badge/Docs-Read%20Now-blue?logo=readthedocs&logoColor=white)](https://geoparse.io)
[![License](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative&logoColor=white)](https://github.com/geoparse/geospatial/blob/main/LICENSE)
[![PythonVersion](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue?style=flat&logo=python&logoColor=white&labelColor=gray)](https://www.python.org/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-Ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Contributors](https://img.shields.io/github/contributors/geoparse/geospatial)](https://github.com/geoparse/geospatial/graphs/contributors)

GeoParse is a Python library designed for the visualization, analysis, and manipulation of vector geospatial data.
It supports a wide range of applications, including telematics analysis (e.g., vehicle trajectories and movement patterns),
footfall analytics (e.g., pedestrian density and mobility trends), and infrastructure mapping (e.g., roads, rivers, water streams, buildings, and lakes).

GeoParse builds on top of popular libraries like GeoPandas and Folium, providing a powerful toolkit for working with geospatial data. It focuses on straightforward
visualization, efficient spatial indexing, advanced geometry manipulations, and user-friendly utilities for handling OpenStreetMap data effortlessly.

---

## ✨ Key Features
* Efficient geospatial indexing using grid-based systems (H3, S2, Geohash)
* Data visualization using Folium maps
* Utilities for working with OpenStreetMap (OSM) data
* Geometry manipulations and conversions between formats
  
---
## ✅ Prerequisites
This repository uses [uv](https://docs.astral.sh/uv/) for Python package and project management.
If you don't have it installed, you need to install it first.
```bash
# On Linux and macOS
curl -LsSf https://astral.sh/uv/install.sh | sh

```

```bash
# On Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

```

You can also update it using the following command if you have `uv` installed:
```
uv self update
```

---
## 🛠️ Installation


`pip install git+https://github.com/geoparse/geoparse.git`

---
## 📖 Ducumentation
The documentation HTML pages are located in `docs/_build/html/`. Open `index.html` to access the documentation, which includes comprehensive descriptions and working examples for each class and function. Additionally, the documentation is available online at [geoparse.io](https://geoparse.io/).

We recommend starting your GeoParse journey with the [tutorial notebooks](https://github.com/geoparse/geoparse/tree/main/tutorials).

---
## 🤝 Contributions Welcome!  

We appreciate contributions from the community! Before submitting a pull request, please:  
✅ Ensure your code passes all tests `(uv run pytest --cov)`.  
✅ Add a new test for any new functionality.  

Thank you for helping improve this project! 🚀  


---
## 🔳 Tile Layers

GeoParse can visualize the vector data on a map with the following tile layers.

| Light                                                                                              | Dark                                                                                       |
| -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | 
| A minimalist, light-colored basemap that serves as a subtle background, emphasizing overlaid data. | A high-contrast, dark-themed map ideal for vibrant data overlays and nighttime aesthetics. |
| <img src="tutorials/graphics/light.png" width="400" height="400">                                  | <img src="tutorials/graphics/dark.png" width="400" height="400">                           |

| Outdoors                                                                                          | Satellite                                                                                       |
| ------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- |
| Designed for outdoor enthusiasts, featuring hiking trails, biking paths, natural landmarks, and elevation contours. | A basemap displaying satellite imagery of the Earth's surface, useful for real-world context and analyses requiring detailed imagery. |
| <img src="tutorials/graphics/outdoors.png" width="400" height="400">                              | <img src="tutorials/graphics/satellite.png" width="400" height="400">                           |

| OSM |
| --- |
| A general-purpose map powered by OpenStreetMap, showcasing roads, buildings, and points of interest. |
| <img src="tutorials/graphics/osm.png" width="400" height="400"> |

---
## GeoParse Classes

GeoParse contains several classes and utility functions designed for geospatial data processing, visualization, and analysis. 
It allows users to handle basic tasks like geocoding and visualization as well as more advanced features like spatial indexing and OSM-based analysis. 
Below is a summary of the main classes and their functionalities:

* **Karta Class**, named after the Swedish word for "map," is designed for creating and customizing interactive maps. It provides methods for adding points, lines, polygons, and choropleth layers to a map. The primary function, `plp` (point, line, polygon), supports various visualization styles and configurations, including clustering, heatmaps, and cell-based layers (e.g., H3, S2, geohash).

* **GeomUtils Class** stands for Geospatial Utilities and provides utility functions for working with geometry objects. Key functionalities include determining UTM projections, transforming geometries between coordinate reference systems (CRS), and calculating geometric statistics like area and perimeter.

* **CellUtils Class** stands for Cell Utilities and provides utility functions for compacting and uncompacting spatial cells like H3, S2, and Geohash. It also supports statistical analysis of spatial cells, such as calculating the number of cells and their area for a given geometry.

* **OSMUtils Class** stands for OSM Utilities and provides utility functions for working with OpenStreetMap (OSM) data and routing engines built on top of OSM. It allows users to retrieve OSM way geometries (either polygons or lines) and perform map matching of GPS coordinates to road networks. 

* **SpatialIndex Class** provides methods for spatial indexing operations. It converts geographic coordinates and geometries into spatial index representations and vice versa, utilizing popular encoding systems like H3, S2, and Geohash. It also accelerates computations through parallel processing, making it useful for efficient spatial queries and handling large datasets.

* **SpatialOps Class** stands for Spatial Operations and provides methods for handling 3D geometries, converting LineStrings to Points, performing spatial intersections, and executing parallelized spatial overlay operations. It also includes utilities for calculating distances (e.g., haversine and Vincenty formulas) and geocoding addresses using the Google Geocoding API.

### Key Features

* **Interactive Map Creation:** Easily create and customize maps with points, lines, polygons, and choropleth layers.
* **Parallel Processing:** Accelerate computations for large datasets using parallelized operations.
* **Geometric Operations:** Perform advanced geometric operations like UTM projection, CRS transformation, and geometric statistics.
* **Spatial Indexing:** Convert coordinates and geometries into spatial cells (Geohash, H3, S2) and vice versa.
* **OSM Integration:** Retrieve and analyze OpenStreetMap data, including way geometries and map matching.
* **Geocoding:** Convert addresses to geographic coordinates and vice versa using the Google Geocoding API.

---
## Examples

You can run [GeoParse examples](https://github.com/geoparse/geoparse/tree/main/tutorials) on MyBinder. No installation required. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/geoparse/geoparse/main?labpath=tutorials%2F00_visualization.ipynb)

### Karta Class

`Karta` is used for visualization and accepts either a pandas `DataFrame` or a GeoPandas `GeoDataFrame` to render geometry data. For a `DataFrame`, the `plp` function in `Karta` automatically identifies columns with names containing "lat" and "lon" (case-insensitive) to use as latitude and longitude for plotting points on the map. If no columns contain these keywords, or if more than two columns contain these keywords, you must explicitly specify the latitude and longitude using the `y` and `x` parameters, respectively, e.g., `plp(df, x="easting", y="northing")`. Note that plp assumes all data is in the [EPSG:4326](https://epsg.io/4326) projection. For a `GeoDataFrame`, the `plp` function can render Shapely objects such as `Point`, `LineString`, `Polygon`, and `MultiPolygon`.


#### Point
In the following example, we demonstrate how to display points from a CSV file, customize the map with point colors and popups, and add layers such as heatmaps and clusters.

```python
import pandas as pd
from geoparse.geoparse import Karta

plp = Karta.plp

df = pd.read_csv("data/great_britain_road_casualties-2023.csv")
df.head()
```
| date       | time  | latitude  | longitude  | number_of_vehicles | number_of_casualties | speed_limit |
|------------|-------|-----------|------------|--------------------|----------------------|-------------|
| 03/01/2023 | 19:12 | 51.356551 | -0.097759  | 1                  | 1                    | 30          |
| 07/01/2023 | 10:05 | 51.593701 | 0.022379   | 1                  | 1                    | 30          |
| 14/01/2023 | 16:15 | 51.466689 | -0.011289  | 1                  | 1                    | 20          |
| 15/01/2023 | 19:51 | 51.671577 | -0.037543  | 1                  | 1                    | 30          |
| 16/01/2023 | 19:22 | 51.447944 | 0.117279   | 2                  | 1                    | 30          |

After loading the data, we can easily display it on a map using `plp(df)`. For a more advanced visualization, we can customize the color of the points based on a feature value (e.g., `speed_limit` in the right map) and create HTML popups that display the attributes of each point.

<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>plp(df)</code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(df, 
    point_color = 'speed_limit', 
    point_popup = {
        'Date': 'date',
        'Number of Casualties': 'number_of_casualties'
    }
)</code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <img src="tutorials/graphics/casualty_map.gif"  height="480">
    </td>
    <td>
      <img src="tutorials/graphics/casualty_popup.gif" height="480">
    </td>
  </tr>
</table>

`plp` can also add heatmap and cluster layers to the map. In the left image, we see the clusters and heatmap of fatal road crashes in Great Britain. If you are working with trajectory data, `plp` can display the direction of movement using the antpath parameter, as shown in the right image.

<table>
  <tr>
    <td style="vertical-align: bottom;">
      <pre><code>plp(df, heatmap=True, cluster=True)</code></pre>
    </td>
    <td style="vertical-align: bottom;">
      <pre><code>
plp(df, antpath=True, line=True)
      </code></pre>
    </td>
  </tr>
  <tr>
    <td>
      <img src="tutorials/graphics/casualty_heatmap_cluster.gif" height="480">
    </td>
    <td>
      <img src="tutorials/graphics/trajectory.gif" height="480">
    </td>
  </tr>
</table>

### Line and Polygon

Using `GeoPandas`, we can read a geospatial file and display its contents using `plp` function. The left image illustrates the border of Luxembourg, represented as a Shapely `Polygon` object. The center image depicts the main roads in Luxembourg, represented as Shapely `LineString` objects. Additionally, `plp` can accept two `GeoDataFrame` objects as a list and display both of them on a single map, as shown in the right image.

| Polygn                                       | LineString                                        | Both                                                            |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf)`                            | `plp(road_gdf, line_color='gold', line_weight=1)` | `plp([border_gdf, road_gdf], line_color='gold', line_weight=1)` |
|![](tutorials/graphics/luxembourg_border.png) | ![](tutorials/graphics/luxembourg_roads.png)      | ![](tutorials/graphics/luxembourg_border_roads.png)             |


Using the `plp` function, we can also add spatial index polygonal layers such as `GeoHash`, Google `S2`, and Uber `H3`. 


| Geohash                                      | S2                                                | H3                                                              |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf, geohash_res=5)`             | `plp(border_gdf, s2_res=11)`                      | `plp(border_gdf, h3_res=6)`                                     |
|![](tutorials/graphics/geohash_5.png)         | ![](tutorials/graphics/s2_11.png)                 | ![](tutorials/graphics/h3_6.png)                                |
 
If the `compact` parameter is set to True, `plp` calculates the parent cell IDs to create a compact representation of the cells.

| Geohash                                      | S2                                                | H3                                                              |
| -------------------------------------------- | --------------------------------------------------|---------------------------------------------------------------- | 
| `plp(border_gdf, geohash_res=7, compact=True)`             | `plp(border_gdf, s2_res=15, compact=True)`                      | `plp(border_gdf, h3_res=10, compact=True)`                                     |
|![](tutorials/graphics/geohash_compact.png)         | ![](tutorials/graphics/s2_compact.png)                 | ![](tutorials/graphics/h3_compact.png)                                |
 

---
### OSM Ways

The `plp` function can also accept `OpenStreetMap` (OSM) Way IDs instead of `DataFrame` or `GeoDataFrame` objects and visualize them as `LineString` or `Polygon` geometries. The left image illustrates two ways of the Tokyo Metro Line represented as `LineString` geometries, while the right image depicts three ways in the Louvre Museum visualized as `Polygon` geometries.

|                                         |                                                                            | 
| --------------------------------------- | ---------------------------------------------------------------------------|
| `plp(osm_ways=[893074361, 666201240])`  | `plp(osm_ways=[335265936, 53820456, 1117218957], s2_res=22, compact=True)` | 
|![](tutorials/graphics/osm_way_line.png) | ![](tutorials/graphics/osm_way_polygon.gif)                                | 


---
