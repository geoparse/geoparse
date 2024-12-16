
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

- **Light**  
  A minimalist, light-colored basemap that serves as a subtle background, emphasizing overlaid data on dashboards and visualizations.
  
- **Dark**  
  A high-contrast, dark-themed map that is ideal for vibrant data overlays and nighttime aesthetics.

- **Outdoors**  
  Designed for outdoor enthusiasts, it features hiking trails, biking paths, natural landmarks, and elevation contours.

- **Satellite**  
  A basemap displaying satellite imagery of the Earth's surface, useful for real-world visual context, remote sensing, or analyses requiring detailed imagery.

  
- **OSM (OpenStreetMap)**  
  A general-purpose map powered by OpenStreetMap, showcasing roads, buildings, and points of interest.

| Light | Dark |
| ----- | ---- | 
| A minimalist, light-colored basemap that serves as a subtle background, emphasizing overlaid data on dashboards and visualizations. | A high-contrast, dark-themed map that is ideal for vibrant data overlays and nighttime aesthetics. |
| <img src="tutorials/graphics/light.png" width="400" height="400"> | <img src="tutorials/graphics/dark.png" width="400" height="400"> |

| Outdoors | Satellite |
| -------- | --------- | 
| Designed for outdoor enthusiasts, it features hiking trails, biking paths, natural landmarks, and elevation contours. |  A basemap displaying satellite imagery of the Earth's surface, useful for real-world visual context, remote sensing, or analyses requiring detailed imagery. |
| <img src="tutorials/graphics/outdoors.png" width="400" height="400"> | <img src="tutorials/graphics/satellite.png" width="400" height="400"> |

| OSM |
| --- |
| A general-purpose map powered by OpenStreetMap, showcasing roads, buildings, and points of interest. |
| <img src="tutorials/graphics/osm.png" width="400" height="400"> |


---

## Examples

You can run [GeoParse examples](https://github.com/geoparse/geoparse/tree/main/tutorials) on MyBinder. No installation required. [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/geoparse/geoparse/main?labpath=tutorials%2F00_visualization.ipynb)

To try the cutting-edge dev version, use this MyBinder link.
### Point
```python
df = pd.read_csv("data/great_britain_road_casualties-2023.csv")
df.head()
```
![](tutorials/graphics/casualty_df.png)

| | |
|-|-|
| `plp(df)` | `plp(df, heatmap=True, cluster=True)` |
|<img src="tutorials/graphics/casualty_map.gif" width="400" height="480"> | <img src="tutorials/graphics/casualty_heatmap_cluster.gif" width="400" height="480"> |


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
