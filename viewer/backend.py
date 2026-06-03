import io
import os
from pathlib import Path
from typing import Annotated, Optional

import pyarrow.parquet as pq
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

app = FastAPI()

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory storage for the loaded dataframe
current_df = None

# Serve HTML files from a directory (e.g., "./html_files")
HTML_DIR = "html_files"
os.makedirs(HTML_DIR, exist_ok=True)
app.mount("/html_files", StaticFiles(directory=HTML_DIR), name="html_files")


def is_csv_gz(filename: str) -> bool:
    """Return True if filename ends with .csv.gz (case-insensitive)."""
    return filename.lower().endswith(".csv.gz")


@app.get("/browse/")
async def browse_directory(path: Optional[str] = None):
    """List directories and supported geospatial files at the given path."""
    root = Path(path) if path else Path.home()

    if not root.exists() or not root.is_dir():
        raise HTTPException(status_code=404, detail=f"Directory not found: {path}")

    supported = {".parquet", ".gpkg", ".geojson", ".shp", ".csv"}

    try:
        entries = sorted(root.iterdir(), key=lambda e: (e.is_file(), e.name.lower()))
        return {
            "current": str(root),
            "parent": str(root.parent) if root != root.parent else None,
            "entries": [
                {
                    "name": e.name,
                    "path": str(e),
                    "is_dir": e.is_dir(),
                    "is_supported": e.is_file() and (e.suffix in supported or is_csv_gz(e.name)),
                }
                for e in entries
                if not e.name.startswith(".")
            ],
        }
    except PermissionError as e:
        raise HTTPException(status_code=403, detail=f"Permission denied: {path}") from e


@app.get("/list_geo_files/")
async def list_geo_files(path: str):
    """Return only supported geospatial files in the given directory."""
    supported = {".parquet", ".gpkg", ".geojson", ".shp", ".csv"}
    directory = Path(path)

    if not directory.exists() or not directory.is_dir():
        raise HTTPException(status_code=404, detail=f"Directory not found: {path}")

    files = sorted(
        [
            {"name": f.name, "path": str(f)}
            for f in directory.iterdir()
            if f.is_file() and (f.suffix in supported or is_csv_gz(f.name))
        ],
        key=lambda x: x["name"],
    )
    return {"directory": str(directory), "files": files}


@app.get("/read_geo_file/")
async def read_geo_file(path: str):
    """Validate that a geo file exists and is readable."""
    supported = {".parquet", ".gpkg", ".geojson", ".shp", ".csv"}
    f = Path(path)

    if not f.exists() or not f.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {path}")
    if not (f.suffix in supported or is_csv_gz(f.name)):
        raise HTTPException(status_code=400, detail=f"Unsupported format: {f.suffix}")

    return {"path": str(f), "name": f.name}


@app.get("/list_html/")
async def list_html_files():
    try:
        html_files = [f for f in os.listdir(HTML_DIR) if f.endswith(".html")]
        return {"html_files": html_files}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.post("/upload/")
async def upload_parquet(file: Annotated[UploadFile, File()]):
    global current_df
    if not file.filename.endswith(".parquet"):
        raise HTTPException(status_code=400, detail="Only Parquet files are allowed")

    try:
        contents = await file.read()
        buffer = io.BytesIO(contents)
        current_df = pq.read_table(buffer).to_pandas()
        return {
            "message": "File uploaded successfully",
            "columns": list(current_df.columns),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) from e


@app.get("/columns/")
async def get_columns():
    global current_df
    if current_df is None:
        raise HTTPException(status_code=404, detail="No file uploaded yet")
    return {"columns": list(current_df.columns)}


@app.get("/value_counts/{column_name}")
async def get_value_counts(column_name: str):
    global current_df
    if current_df is None:
        raise HTTPException(status_code=404, detail="No file uploaded yet")

    if column_name not in current_df.columns:
        raise HTTPException(status_code=404, detail=f"Column '{column_name}' not found")

    value_counts = current_df[column_name].value_counts().to_dict()
    return {"column": column_name, "value_counts": value_counts}


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
