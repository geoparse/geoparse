
<h1 align="center"><img align="center" src="tutorials/data/geoparse_logo.png" alt="GeoParse Logo" width="200"/></h1>
<h1 align="center">GeoParse</h1>
<h3 align="center">It's all about points <img src="tutorials/data/point.png" width="10"/> lines <img src="tutorials/data/line.png" width="40"/> and polygons <img src="tutorials/data/polygon.png" width="30"/></h3>

---

[![GeoParse](https://img.shields.io/badge/GeoParse-008000.svg)](https://geoparse.io)
[![Docs](https://img.shields.io/badge/Docs-Read%20Now-blue?logo=readthedocs&logoColor=white)](https://geo-parse.readthedocs.io/en/latest/)
[![License](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative&logoColor=white)](https://github.com/geoparse/geospatial/blob/main/LICENSE)
[![PythonVersion](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue?style=flat&logo=python&logoColor=white&labelColor=gray)](https://www.python.org/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-Ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Contributors](https://img.shields.io/github/contributors/geoparse/geospatial)](https://github.com/geoparse/geospatial/graphs/contributors)


GeoParse is a Python library designed for the visualization, analysis, and manipulation of vector geospatial data. It builds on top of popular libraries like GeoPandas and Folium, providing a powerful toolkit for working with geospatial data. GeoParse focuses on efficient geospatial indexing, geometry manipulations, and utilities to handle OpenStreetMap data with ease.

---

## Key Features
* Efficient geospatial indexing using grid-based systems (H3, S2, Geohash)
* Data visualization using Folium maps
* Utilities for working with OpenStreetMap (OSM) data
* Geometry manipulations and conversions between formats
  
---

## Installation


`pip install git+https://github.com/geoparse/geospatial.git`

---

## Ducumentation
We recommend starting your GeoParse journey with the [tutorial notebooks](https://github.com/geoparse/geoparse/tree/main/tutorials).

The official API documentation is hosted on [ReadTheDocs](https://geo-parse.readthedocs.io/en/latest/)

---
## Tile Layers

GeoParse can visualize the vector data on a map with the following tile layers.

| Light                                                                                               | Dark                                                                                                |
| --------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- | 
| A minimalist, light-colored basemap that serves as a subtle background, emphasizing overlaid data. | A high-contrast, dark-themed map ideal for vibrant data overlays and nighttime aesthetics.         |
| <img src="tutorials/graphics/light.png" width="400" height="400">                                  | <img src="tutorials/graphics/dark.png" width="400" height="400">                                  |

| Outdoors                                                                                          | Satellite                                                                                       |
| ------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| Designed for outdoor enthusiasts, featuring hiking trails, biking paths, natural landmarks, and elevation contours. | A basemap displaying satellite imagery of the Earth's surface, useful for real-world context and analyses requiring detailed imagery. |
| <img src="tutorials/graphics/outdoors.png" width="400" height="400">                              | <img src="tutorials/graphics/satellite.png" width="400" height="400">                           |

| OSM |
| --- |
| A general-purpose map powered by OpenStreetMap, showcasing roads, buildings, and points of interest. |
| <img src="tutorials/graphics/osm.png" width="400" height="400"> |

---
## GeoParse Modules

GeoParse provides a powerful suite of tools for anyone working with geospatial data. It allows users to handle basic tasks like geocoding and visualization as well as more advanced features like spatial indexing and OSM-based analysis. It includes four key modules, each focusing on different aspects of geospatial tasks.

* **karta**: The [karta](https://github.com/geoparse/geoparse/blob/main/geoparse/karta.py) module, named after the Swedish word for "map," is used to visualize vector data (points, lines, and polygons) on a map. The primary function, `plp` (point, line, polygon), supports various visualization styles and configurations, including clustering, heatmaps, choropleth maps, and cell-based layers (e.g., H3, S2, geohash). `karta` is ideal for creating interactive visual maps from geospatial data.

* **gutils**: The [gutils](https://github.com/geoparse/geoparse/blob/main/geoparse/gutils.py) stands for Geospatial Utilities and provides a wide range of geospatial utilities that serve various tasks such as geometric statistics, transformations, geocoding (address to coordinates), and distance or area calculations.

* **gindex**: The [gindex](https://github.com/geoparse/geoparse/blob/main/geoparse/gindex.py) module stands for Geospatial Index and specializes in spatial indexing systems. It converts geographic coordinates and geometries into spatial index representations and vice versa, utilizing popular encoding systems like Geohash, H3, and S2. It also accelerates computations through parallel processing, making it useful for efficient spatial queries and handling large datasets.

* **osmtools**: The [osmtools](https://github.com/geoparse/geoparse/blob/main/geoparse/osmtools.py) stands for OpenStreetMap Tools and offers tools to work with OpenStreetMap (OSM) data and routing engines built on top of OSM. It allows users to retrieve OSM way geometries (either polygons or lines), decode encoded polyline strings (such as those used by Valhalla and similar services), and perform map matching of GPS coordinates to road networks. This module is essential for extracting, analyzing, and visualizing OSM data as well as processing GPS traces with routing services.


---
## Examples

You can run [GeoParse examples](https://github.com/geoparse/geoparse/tree/main/tutorials) on MyBinder. No installation required. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/geoparse/geoparse/main?labpath=tutorials%2F00_visualization.ipynb)

### karta

`karta` is used for visualization and accepts either a pandas `DataFrame` or a GeoPandas `GeoDataFrame` to render geometry data. For a `DataFrame`, the `plp` function in `karta` identifies columns with names containing `lat` and `lon` to use as latitude and longitude for displaying points on the map. For a `GeoDataFrame`, the `plp` function can render Shapely objects such as `Point`, `LineString`, `Polygon`, and `MultiPolygon`.


#### Point
```python
import pandas as pd
from geoparse.karta import plp

df = pd.read_csv("data/great_britain_road_casualties-2023.csv")
df.head()
```
| date       | time   | latitude   | longitude   | number_of_vehicles  | number_of_casualties  | speed_limit |
|------------|--------|------------|-------------|---------------------|-----------------------|-------------|
| 03/01/2023 | 19:12  | 51.356551  | -0.097759   | 1                   | 1                     | 30          |
| 07/01/2023 | 10:05  | 51.593701  | 0.022379    | 1                   | 1                     | 30          |
| 14/01/2023 | 16:15  | 51.466689  | -0.011289   | 1                   | 1                     | 20          |
| 15/01/2023 | 19:51  | 51.671577  | -0.037543   | 1                   | 1                     | 30          |
| 16/01/2023 | 19:22  | 51.447944  | 0.117279    | 2                   | 1                     | 30          |

After loading the data, we can easily display it on a map using plp(df). For a more advanced visualization, we can customize the color of the points based on a feature value (e.g., `speed_limit` in the right map) and create HTML popups that display the attributes of each point.

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
      <img src="tutorials/graphics/casualty_map.gif"  width="314" height="480">
    </td>
    <td>
      <img src="tutorials/graphics/casualty_popup.gif" height="480">
    </td>
  </tr>
</table>


### Polygon

```python
gdf = gpd.read_file("data/london.geojson")
plp(gdf)
```
![](tutorials/graphics/london.png)

---
### OSM Ways
```python
plp(osm_ways=[335265936, 53820456, 1117218957], s2_res=22, compact=True)
```
![London](tutorials/graphics/osm_ways_2.png)

---
