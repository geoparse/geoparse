import collections
import json
import os
from datetime import datetime
from math import sqrt
from multiprocessing import Pool, cpu_count
from time import time
from typing import List, Tuple, Union

import geopandas as gpd
import numpy as np
import pygeohash
from h3 import h3
from polygon_geohasher.polygon_geohasher import geohash_to_polygon, polygon_to_geohashes
from s2 import s2
from shapely.geometry import MultiPolygon, Polygon
from shapely.geometry.base import BaseGeometry


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
        return [h3.geo_to_h3(lat, lon, res) for lat, lon in zip(lats, lons)]
    else:
        raise ValueError(f"Unsupported cell type: {cell_type}. Choose 'geohash', 's2', 's2_int', or 'h3'.")


def cellpoint(cells, cell_type):
    if cell_type == "s2":
        return [(s2.CellId(cell).to_lat_lng().lat().degrees, s2.CellId(cell).to_lat_lng().lng().degrees) for cell in cells]


def polycell(geoms: List[Union[Polygon, MultiPolygon]], cell_type: str, res: int, dump: str = None) -> Union[List[str], None]:
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

    cells = []
    if cell_type == "geohash":
        cells = set()
        for geom in geoms:
            cells |= polygon_to_geohashes(geom, precision=res, inner=False)  # Collect Geohashes for each Polygon
        cells = list(cells)

    elif cell_type == "s2":
        for poly in polys:
            cells += s2.polyfill(poly, res, geo_json_conformant=True, with_id=True)  # Use S2 to fill each Polygon
        cells = [item["id"] for item in cells]  # Keep only the cell IDs
        cells = list(set(cells))  # Remove duplicates

    elif cell_type == "h3":
        for poly in polys:
            cells += h3.polyfill(poly, res, geo_json_conformant=True)  # Use H3 to fill each Polygon

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
    inputs = zip(geom_chunks, [cell_type] * 4 * n_cores, [res] * 4 * n_cores, [dump] * 4 * n_cores)

    # Parallel processing to generate cells
    if dump:
        with Pool(n_cores) as pool:
            pool.starmap(polycell, inputs)
        if verbose:
            elapsed_time = round(time() - start_time)
            print(f"{elapsed_time} seconds.")
        return
    else:
        with Pool(n_cores) as pool:
            cells = pool.starmap(polycell, inputs)
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
            cells = compact_cells(cells, cell_type)
            if verbose:
                elapsed_time = round(time() - start_time)
                print(f"{elapsed_time} seconds.")

        return cells, cell_counts


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
        else cell[1]
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
        else Polygon(h3.h3_to_geo_boundary(cell, geo_json=True))
        for cell in cells
    ]

    return res, geoms


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
        return list(h3.compact(cells))
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
    cells = polycell(geom, cell="h3", res=h3_res)
    area = h3.hex_area(h3_res, unit="km^2")
    if compact:
        cells = h3.compact(cells)
    return len(cells), area
