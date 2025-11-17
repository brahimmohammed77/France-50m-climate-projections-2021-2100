"""
Python Script for Mosaicking and Climatology Integration

Purpose
-------
- Read multiple tile rasters (e.g. anomalies for a given month/year).
- Mosaic them together.
- Align the mosaic to a reference climatology raster (8110 normal).
- Add the anomaly mosaic to the reference raster to obtain a final field.
- Save a single GeoTIFF per month/year.

Notes
-----
- VARIABLE can be "pr" (precipitation), "tas" (mean temperature),
  "tasmin" (min temperature), "tasmax" (max temperature), etc.
- REF_PREFIX typically refers to the climatology period, e.g. "8110" for 1981–2010.
"""

import os
import glob
import numpy as np
import rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling
from pathlib import Path
import time
from pyproj import datadir
from multiprocessing import Pool
import warnings
import contextlib


# ---------------------------------------------------------------------------
# CONFIGURATION (edit as needed)
# ---------------------------------------------------------------------------

class Config:
    # Input directories
    # INPUT_DIR contains subfolders like: anomaly_pr_YYYY_MM_...
    INPUT_DIR = Path("PATH/TO/INPUT/RESULTS/")

    # REF_DIR contains reference rasters (climatology) like:
    #   pr_8110_1.tif, pr_8110_2.tif, ..., pr_8110_12.tif
    REF_DIR = Path("PATH/TO/REFERENCE/RASTERS/")

    # Output directory for monthly mosaics
    OUTPUT_DIR = Path("PATH/TO/OUTPUT/MOSAICS/")

    # Variable name:
    #   pr     = precipitation
    #   tas    = mean temperature
    #   tasmin = min temperature
    #   tasmax = max temperature
    VARIABLE = "pr"

    # Prefix used in reference rasters (e.g. 8110 for 1981–2010 normals)
    REF_PREFIX = "8110"

    # Optional tags in output filename
    MODEL = "MODEL"          # climate model tag
    SCENARIO = "SCENARIO"    # scenario tag (e.g., ssp245, ssp585)

    # Parallel settings
    NUM_CORES = 6

    # Nodata value for output
    NODATA = -9999


# Initialization
os.environ["PROJ_LIB"] = datadir.get_data_dir()
Config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# PROCESSING OF ONE SUBDIRECTORY
# ---------------------------------------------------------------------------

def process_directory(dossier_path: str) -> str:
    """
    Process one folder containing tile anomalies for a given YYYY_MM.

    Steps:
    - Merge all tile .tif files into a global mosaic.
    - Align mosaic to reference raster (same CRS, grid, transform).
    - Add mosaic (anomaly) to reference climatology.
    - Save final GeoTIFF.
    """
    start_time = time.time()
    name = os.path.basename(dossier_path)

    try:
        # Expected folder pattern: anomaly_<var>_<year>_<month>...
        parts = name.split("_")
        if len(parts) < 4:
            return f"[IGNORED] Invalid folder name: {name}"

        year = parts[2]
        month = parts[3].zfill(2)

        out_file = Config.OUTPUT_DIR / f"{Config.VARIABLE}_{year}_{month}_{Config.MODEL}_{Config.SCENARIO}.tif"
        if out_file.exists():
            return f"[SKIP] {out_file.name} already exists"

        # Reference raster for this month (e.g., pr_8110_1.tif)
        ref_raster = Config.REF_DIR / f"{Config.VARIABLE}_{Config.REF_PREFIX}_{int(month)}.tif"
        if not ref_raster.exists():
            return f"[ERROR] Missing reference raster: {ref_raster.name}"

        # All tile TIF files in this folder
        tif_files = list(Path(dossier_path).glob("*.tif"))
        if not tif_files:
            return f"[EMPTY] No TIF files in {name}"

        # --- 1. Merge TIFs into a mosaic ---
        with contextlib.ExitStack() as stack:
            srcs = [stack.enter_context(rasterio.open(p)) for p in tif_files]
            mosaic, mosa_transform = merge(
                srcs,
                nodata=Config.NODATA,
                res=min(src.res[0] for src in srcs)
            )

        mosa_crs = srcs[0].crs
        mosaic = mosaic[0].astype(np.float32)
        mosaic[mosaic == Config.NODATA] = np.nan

        # --- 2. Load reference raster (climatology) ---
        with rasterio.open(ref_raster) as ref:
            ref_arr = ref.read(1).astype(np.float32)
            ref_arr[ref_arr == ref.nodata] = np.nan
            ref_meta = ref.meta.copy()

            # --- 3. Reproject/align mosaic if necessary ---
            need_reproject = (
                (mosa_crs != ref.crs) or
                (mosaic.shape != ref_arr.shape) or
                (not np.allclose(mosa_transform, ref.transform))
            )

            if need_reproject:
                aligned_mosaic = np.full_like(ref_arr, np.nan)
                reproject(
                    source=mosaic,
                    destination=aligned_mosaic,
                    src_transform=mosa_transform,
                    src_crs=mosa_crs,
                    dst_transform=ref.transform,
                    dst_crs=ref.crs,
                    resampling=Resampling.nearest
                )
                mosaic = aligned_mosaic
                mosa_transform = ref.transform
                mosa_crs = ref.crs

        # --- 4. Combine anomaly mosaic with reference (sum on overlap) ---
        valid_mask = ~np.isnan(mosaic) & ~np.isnan(ref_arr)
        result = np.full_like(mosaic, Config.NODATA, dtype=np.float32)
        result[valid_mask] = mosaic[valid_mask] + ref_arr[valid_mask]

        # --- 5. Compute basic stats for logging ---
        valid_vals = result[result != Config.NODATA]
        if valid_vals.size > 0:
            stats_msg = (
                f"min={valid_vals.min()}, "
                f"max={valid_vals.max()}, "
                f"mean={valid_vals.mean():.1f}"
            )
        else:
            stats_msg = "no valid values"

        # --- 6. Save result as GeoTIFF ---
        result = np.round(result).astype(np.int32)

        ref_meta.update({
            "driver": "GTiff",
            "height": result.shape[0],
            "width": result.shape[1],
            "transform": mosa_transform,
            "crs": mosa_crs,
            "nodata": Config.NODATA,
            "dtype": "int32",
            "compress": "DEFLATE",
            "predictor": 2,
        })

        with rasterio.open(out_file, "w", **ref_meta) as dst:
            dst.write(result, 1)

        elapsed = time.time() - start_time
        return f"[OK] {out_file.name} ({result.size / 1e6:.1f} Mpix, {stats_msg}) in {elapsed:.1f}s"

    except Exception as e:
        return f"[ERROR] {name}: {str(e)}"


# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

def main():
    print(f"\n🚀 Processing with {Config.NUM_CORES} core(s)")

    # Folders like: anomaly_pr_YYYY_MM*
    dossiers = sorted(
        glob.glob(str(Config.INPUT_DIR / f"anomaly_{Config.VARIABLE}_*"))
    )

    print(f"🎯 {len(dossiers)} folder(s) to process\n")

    start_time = time.time()

    with Pool(processes=Config.NUM_CORES) as pool:
        results = [pool.apply_async(process_directory, (d,)) for d in dossiers]
        for res in results:
            print(res.get())

    print(f"\n✅ Finished in {time.time() - start_time:.1f} seconds")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()
