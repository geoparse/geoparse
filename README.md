
<h1 align="center"><img align="center" src="https://geoparse.io/graphics/geoparse_logo.png" alt="GeoParse Logo" width="200"/></h1>
<h1 align="center">GeoParse</h1>
<h3 align="center">All About Points <img src="https://geoparse.io/graphics/point.png" width="10"/> Lines <img src="https://geoparse.io/graphics/line.png" width="40"/> and Polygons <img src="https://geoparse.io/graphics/polygon.png" width="30"/></h3>

---

[![GeoParse](https://img.shields.io/badge/GeoParse-008000.svg)](https://geoparse.io)
[![Docs](https://img.shields.io/badge/Docs-Read%20Now-blue?logo=readthedocs&logoColor=white)](https://geoparse.io)
[![License](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative&logoColor=white)](https://github.com/geoparse/geospatial/blob/main/LICENSE)
[![PythonVersion](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue?style=flat&logo=python&logoColor=white&labelColor=gray)](https://www.python.org/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-Ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Contributors](https://img.shields.io/github/contributors/geoparse/geospatial)](https://github.com/geoparse/geospatial/graphs/contributors)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

GeoParse is an open-source project that provides tools for the visualization, analysis, and manipulation of vector geospatial data.
It supports a wide range of applications, including
footfall analytics (e.g., pedestrian density and mobility trends),
telematics analysis (e.g., vehicle trajectories and movement patterns), and
infrastructure mapping (e.g., roads, rivers, water streams, buildings, and lakes).

GeoParse builds on top of popular libraries like GeoPandas and folium, providing a powerful toolkit for working with geospatial data.
It focuses on straightforward visualization, efficient spatial indexing, advanced geometry manipulations, and user-friendly utilities for handling OpenStreetMap data effortlessly.

With GeoParse, you can:

1. **Create Interactive Maps** – Generate interactive maps with multiple tile layers (e.g., OpenStreetMap, Satellite, Dark Mode) using [folium](https://github.com/python-visualization/folium).
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

# Tutorials
We recommend starting your GeoParse journey with the HTML versions of the tutorial notebooks, available [here](https://geoparse.io/tutorials/karta.html).
The notebooks cover a wide range of use cases to help you get familiar with GeoParse’s capabilities. 
You can find them on [GitHub](https://github.com/geoparse/geoparse/tree/main/tutorials). 
To run the notebooks, use the following command after setting up the environment described in the [contributing](https://github.com/geoparse/geoparse/blob/main/CONTRIBUTING.md) guide on GitHub:
```sh
jupyter lab 
```

# Contributing to GeoParse
We appreciate all contributions, including bug reports, fixes, documentation improvements, feature enhancements, and ideas. 
A detailed overview on how to contribute can be found in the [contributing](https://github.com/geoparse/geoparse/blob/main/CONTRIBUTING.md) guide on GitHub.

# License and Credits
GeoParse is licensed under the [MIT license](https://github.com/geoparse/geoparse/blob/main/LICENSE) and is written and maintained by Abbas Kiasari (abbas@geoparse.io).
