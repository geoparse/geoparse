"""
GeoParse Viewer - A Streamlit application for visualizing geospatial data.

Supports GeoParquet, GeoPackage, GeoJSON, and Shapefile formats.



You need to install geoparse as an editable (development) package, so that changes to the source code are immediately reflected without reinstalling.

uv pip install -e .
"""

import base64
import logging
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd
import requests
import streamlit as st
from streamlit.runtime.uploaded_file_manager import UploadedFile

from geoparse.geoparse import GeomUtils, plp

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
TARGET_CRS = "EPSG:4326"
SUPPORTED_EXTENSIONS = (".parquet", ".gpkg", ".geojson", ".shp", ".csv", ".csv.gz")
MAP_HEIGHT = 600
PREVIEW_ROWS = 3

BACKEND_URL = "http://localhost:8000"


def configure_page() -> None:
    """Configure Streamlit page settings."""
    st.set_page_config(
        layout="wide",
        page_title="GeoParse Viewer",
        page_icon="graphics/geoparse_logo.png",
    )


def inject_custom_css() -> None:
    """Inject custom CSS styles into the Streamlit app."""
    st.markdown(
        """
        <style>
            .metric-box {
                text-align: center;
                padding: 15px;
                border-radius: 10px;
                margin: -15px 0;
            }
            .metric-label {
                font-size: 16px;
                color: #000;
            }
            .metric-value {
                font-size: 24px;
                font-weight: bold;
                color: #000;
            }
            .spacer {
                height: 75px;
            }
            .centered-header {
                text-align: center;
                margin-bottom: 2rem;
            }
            .centered-header img {
                width: 60px;
                height: 60px;
                display: block;
                margin: 0 auto 10px auto;
            }
            .centered-header h1 {
                margin: 0;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )


def render_header(logo_path: Path) -> None:
    """Render the centered logo and title header."""
    logo_html = ""
    if logo_path.exists():
        try:
            with open(logo_path, "rb") as img_file:
                b64 = base64.b64encode(img_file.read()).decode()
            logo_html = f'<img src="data:image/png;base64,{b64}" alt="GeoParse Logo">'
        except IOError as e:
            logger.warning(f"Failed to load logo: {e}")

    st.markdown(
        f"""
        <div class="centered-header">
            {logo_html}
            <h1>GeoParse Viewer</h1>
        </div>
        """,
        unsafe_allow_html=True,
    )


def reproject_if_needed(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reproject GeoDataFrame to target CRS if necessary.

    Args:
        gdf: Input GeoDataFrame.

    Returns:
        GeoDataFrame in target CRS (EPSG:4326).
    """
    if gdf.crs and gdf.crs.to_epsg() != 4326:
        source_epsg = gdf.crs.to_epsg()
        st.info(f"Reprojecting from EPSG:{source_epsg} to {TARGET_CRS}")
        return gdf.to_crs(TARGET_CRS)
    return gdf


def render_geodataframe(gdf: gpd.GeoDataFrame) -> Optional[str]:
    """
    Render a GeoDataFrame with preview and map visualization.

    Args:
        gdf: GeoDataFrame to visualize.

    Returns:
        HTML string of the map if successful, otherwise None.
    """
    with st.expander("Show data preview"):
        st.dataframe(gdf.head(PREVIEW_ROWS))

    try:
        map_obj = plp(gdf)
        map_html = map_obj._repr_html_()
        st.components.v1.html(map_html, height=MAP_HEIGHT)
        return map_html
    except Exception as e:
        logger.error(f"Visualization failed: {e}")
        st.error(f"Visualization failed: {e}")
        return None


def process_uploaded_file(uploaded_file: UploadedFile) -> Optional[gpd.GeoDataFrame]:
    """
    Read and process an uploaded geospatial file.

    Args:
        uploaded_file: Streamlit UploadedFile object.

    Returns:
        Processed GeoDataFrame, or None if reading fails.
    """
    try:
        ext = Path(uploaded_file).suffix.lower()
        if ext in [".csv", ".gz"]:
            gdf = pd.read_csv(uploaded_file)
            gdf = GeomUtils.data_to_geoms(gdf)
        elif ext == ".parquet":
            gdf = gpd.read_parquet(uploaded_file)
        else:
            gdf = gpd.read_file(uploaded_file)

        return reproject_if_needed(gdf)
    except Exception as e:
        logger.error(f"Failed to read file {uploaded_file.name}: {e}")
        st.error(f"Failed to read file: {e}")
        return None


def render_file_picker(key: str, label: str) -> Optional[str]:
    """
    Render a self-contained GUI directory browser and file picker.

    Args:
        key: Unique prefix for session state and widget keys.
        label: Display label shown above the picker.

    Returns:
        The confirmed file path string, or None if none selected yet.
    """
    nav_key = f"{key}_nav_path"
    selected_key = f"{key}_selected_file"

    if nav_key not in st.session_state:
        st.session_state[nav_key] = None
    if selected_key not in st.session_state:
        st.session_state[selected_key] = None

    st.markdown(f"**{label}**")

    with st.expander(
        "📂 Browse…" if not st.session_state[selected_key] else f"📄 `{Path(st.session_state[selected_key]).name}`",
        expanded=st.session_state[selected_key] is None,
    ):
        browse_data = fetch_directory(st.session_state[nav_key])
        if not browse_data:
            return st.session_state[selected_key]

        current_path = browse_data["current"]
        parent_path = browse_data["parent"]
        entries = browse_data["entries"]

        # Current path + Up button
        path_col, up_col = st.columns([0.82, 0.18])
        with path_col:
            st.code(current_path, language=None)
        with up_col:
            if parent_path and st.button("⬆ Up", key=f"{key}_up"):
                st.session_state[nav_key] = parent_path
                st.rerun()

        # Folders
        dirs = [e for e in entries if e["is_dir"]]
        if dirs:
            st.caption("Folders")
            for d in dirs:
                if st.button(f"📁  {d['name']}", key=f"{key}_dir_{d['path']}"):
                    st.session_state[nav_key] = d["path"]
                    st.rerun()

        # Supported files as selectable buttons
        files = [e for e in entries if e["is_supported"]]
        if files:
            st.caption("Files")
            for f in files:
                is_selected = st.session_state[selected_key] == f["path"]
                label_text = f"✅ `{f['name']}`" if is_selected else f"📄 `{f['name']}`"
                if st.button(label_text, key=f"{key}_file_{f['path']}"):
                    st.session_state[selected_key] = f["path"]
                    st.rerun()
        else:
            st.caption("No supported geospatial files in this folder.")

        skipped = [e for e in entries if not e["is_dir"] and not e["is_supported"]]
        if skipped:
            st.caption(f"{len(skipped)} unsupported file(s) hidden")

        # Clear selection
        if st.session_state[selected_key]:
            st.divider()
            if st.button("✖ Clear selection", key=f"{key}_clear"):
                st.session_state[selected_key] = None
                st.rerun()

    return st.session_state[selected_key]


def load_geo_file(path: str) -> Optional[gpd.GeoDataFrame]:
    """
    Read a geospatial file from a local path via the backend validation,
    then load it with geopandas.

    Args:
        path: Absolute path to the geospatial file.

    Returns:
        Reprojected GeoDataFrame, or None on failure.
    """
    try:
        # Validate via backend first
        r = requests.get(
            f"{BACKEND_URL}/read_geo_file/",
            params={"path": path},
            timeout=5,
        )
        r.raise_for_status()
    except requests.HTTPError as e:
        st.error(f"Backend rejected file: {e.response.json().get('detail', e)}")
        return None
    except Exception as e:
        st.error(f"Backend error: {e}")
        return None

    try:
        ext = Path(path).suffix.lower()
        if ext in [".csv", ".gz"]:
            gdf = pd.read_csv(path)
            gdf = GeomUtils.data_to_geoms(gdf)
        elif ext == ".parquet":
            gdf = gpd.read_parquet(path)
        else:
            gdf = gpd.read_file(path)

        return reproject_if_needed(gdf)
    except Exception as e:
        logger.error(f"Failed to read {path}: {e}")
        st.error(f"Failed to read file: {e}")
        return None


def render_compare_tab() -> None:
    """Render the file comparison tab with side-by-side GUI file pickers."""
    col1, col2 = st.columns(2)

    with col1:
        selected_1 = render_file_picker(key="compare_left", label="First file")
        if selected_1:
            gdf = load_geo_file(selected_1)
            if gdf is not None:
                render_geodataframe(gdf)

    with col2:
        selected_2 = render_file_picker(key="compare_right", label="Second file")
        if selected_2:
            gdf = load_geo_file(selected_2)
            if gdf is not None:
                render_geodataframe(gdf)


def filter_valid_files(uploaded_files: list[UploadedFile]) -> list[UploadedFile]:
    """
    Filter uploaded files to only include supported geospatial formats.

    Args:
        uploaded_files: List of uploaded files.

    Returns:
        Filtered list containing only supported file types.
    """
    return [f for f in uploaded_files if f.name.endswith(SUPPORTED_EXTENSIONS)]


def fetch_directory(path: Optional[str] = None) -> Optional[dict]:
    """Call the backend browse endpoint."""
    try:
        params = {"path": path} if path else {}
        r = requests.get(f"{BACKEND_URL}/browse/", params=params, timeout=5)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        st.error(f"Backend error: {e}")
        return None


def fetch_geo_files(path: str) -> list[dict]:
    """Return supported geospatial files in a directory via backend."""
    try:
        r = requests.get(f"{BACKEND_URL}/list_geo_files/", params={"path": path}, timeout=5)
        r.raise_for_status()
        return r.json().get("files", [])
    except Exception as e:
        st.error(f"Backend error: {e}")
        return []


def render_file_browser_tab() -> None:
    """Render the file browser tab with a GUI directory picker via backend."""
    st.subheader("File Browser")

    # --- Session state defaults ---
    if "browser_dir" not in st.session_state:
        st.session_state.browser_dir = None  # confirmed working directory
    if "nav_path" not in st.session_state:
        st.session_state.nav_path = None  # path being browsed in the picker

    # ── Directory Picker ──────────────────────────────────────────────────────
    with st.expander("📂 Choose folder", expanded=st.session_state.browser_dir is None):
        browse_data = fetch_directory(st.session_state.nav_path)
        if not browse_data:
            return

        current_path = browse_data["current"]
        parent_path = browse_data["parent"]
        entries = browse_data["entries"]

        # Current path + Up button on the same row
        path_col, up_col = st.columns([0.85, 0.15])
        with path_col:
            st.code(current_path, language=None)
        with up_col:
            if parent_path and st.button("⬆ Up"):
                st.session_state.nav_path = parent_path
                st.rerun()

        # Directory listing
        dirs = [e for e in entries if e["is_dir"]]
        files = [e for e in entries if e["is_supported"]]
        skipped = [e for e in entries if not e["is_dir"] and not e["is_supported"]]

        if dirs:
            st.caption("Folders")
            for d in dirs:
                if st.button(f"📁  {d['name']}", key=f"dir_{d['path']}"):
                    st.session_state.nav_path = d["path"]
                    st.rerun()

        if files:
            st.caption(f"{len(files)} supported file(s) in this folder")
            for f in files:
                st.markdown(f"&nbsp;&nbsp;📄 `{f['name']}`")
        elif not dirs:
            st.warning("No supported geospatial files here.")

        if skipped:
            st.caption(f"{len(skipped)} unsupported file(s) hidden")

        # Confirm selection
        st.divider()
        btn_col, _ = st.columns([0.4, 0.6])
        with btn_col:
            if st.button("✅ Select this folder", disabled=not files):
                st.session_state.browser_dir = current_path
                st.rerun()

    # ── File List + Visualisation ─────────────────────────────────────────────
    if not st.session_state.browser_dir:
        st.info("👆 Open **Choose folder** above and confirm a directory to begin.")
        return

    geo_files = fetch_geo_files(st.session_state.browser_dir)

    if not geo_files:
        st.warning("No supported geospatial files found in the selected directory.")
        return

    left_col, right_col = st.columns([0.3, 0.7])

    with left_col:
        st.markdown("**Available Files**")
        options = [f["path"] for f in geo_files]
        selected_path = st.radio(
            "Select a file to visualize",
            options=options,
            format_func=lambda p: Path(p).name,
            index=0,
            key="selected_file_radio",
            label_visibility="collapsed",
        )

    with right_col:
        if selected_path:
            st.markdown(f"**Visualizing:** `{Path(selected_path).name}`")
            try:
                ext = Path(selected_path).suffix.lower()
                if ext in [".csv", ".gz"]:
                    gdf = pd.read_csv(selected_path)
                    gdf = GeomUtils.data_to_geoms(gdf)
                elif ext == ".parquet":
                    gdf = gpd.read_parquet(selected_path)
                else:
                    gdf = gpd.read_file(selected_path)

                gdf = reproject_if_needed(gdf)
                map_html = render_geodataframe(gdf)

                # ── SAVE MAP BUTTON ──────────────────────────────────────────
                if map_html:
                    # Generate a filename for the exported HTML
                    base_name = Path(selected_path).stem
                    download_filename = f"{base_name}.html"

                    st.download_button(
                        label="💾 Save Map as HTML",
                        data=map_html,
                        file_name=download_filename,
                        mime="text/html",
                        key="save_map_button",
                    )
            except Exception as e:
                logger.error(f"Failed to read {selected_path}: {e}")
                st.error(f"Failed to read file: {e}")


def main() -> None:
    """Main application entry point."""
    configure_page()
    inject_custom_css()
    render_header(Path("graphics/geoparse_logo.png"))

    tab_browser, tab_compare = st.tabs(["📊 File Browser", "📁 Compare Files"])

    with tab_compare:
        render_compare_tab()

    with tab_browser:
        render_file_browser_tab()


if __name__ == "__main__":
    main()
