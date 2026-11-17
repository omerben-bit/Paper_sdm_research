"""
Render currently available public bathymetry for Israeli Mediterranean waters.
Source: NOAA ETOPO 2022 (15 arc-sec) via ERDDAP — free, global, depth+land.
Also attempts Israel EEZ multibeam (Kanari et al. 2020) if zip is present.
"""
from __future__ import annotations

import io
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import urllib.request

OUT = Path(__file__).resolve().parent / "outputs_israel_sea_3d"
OUT.mkdir(exist_ok=True)

# Levantine / Israel EEZ-ish box
LAT_MIN, LAT_MAX = 29.5, 34.5
LON_MIN, LON_MAX = 32.0, 36.0

# Subsample stride so 15s grid stays plottable
STRIDE = 4

ERDDAP = (
    "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ETOPO_2022_v1_15s.nc"
    f"?z[({LAT_MIN}):{STRIDE}:({LAT_MAX})][({LON_MIN}):{STRIDE}:({LON_MAX})]"
)


def download_etopo() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    print("Downloading ETOPO 2022 subset (ERDDAP)...")
    print(ERDDAP)
    try:
        import xarray as xr
    
    except ImportError:
        raise SystemExit("Please: pip install xarray netCDF4 matplotlib numpy")

    # urllib may need a User-Agent on some networks
    req = urllib.request.Request(ERDDAP, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        data = resp.read()
    ds = xr.open_dataset(io.BytesIO(data))
    z = ds["z"].values.astype(float)
    lat = ds["latitude"].values.astype(float)
    lon = ds["longitude"].values.astype(float)
    ds.close()
    # ERDDAP often returns (lat, lon)
    if z.ndim != 2:
        raise RuntimeError(f"Unexpected z shape: {z.shape}")
    Lon, Lat = np.meshgrid(lon, lat)
    return Lon, Lat, z


def plot_map(Lon, Lat, z, title: str, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 7))
                    cv c d c cv cv cv    # OceanDepth-like: land positive, sea negative
    levels = np.linspace(-2500, 800, 60)
    cf = ax.contourf(Lon, Lat, z, levels=levels, cmap="terrain", extend="both")
    ax.contour(Lon, Lat, z, levels=[0], colors="k", linewidths=0.8)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title(title)
    ax.set_aspect("equal")
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85)
    cbar.set_label("Elevation / depth (m)")
    fig.tight_layout()
    path = OUT / fname
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("Wrote", path)


def plot_3d(Lon, Lat, z, title: str, fname: str, elev=35, azim=-60) -> None:
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    # Vertical exaggeration for visibility of shelf break
    z_plot = np.where(z < 0, z * 3.0, z)
    norm = plt.Normalize(vmin=-2500, vmax=800)
    surf = ax.plot_surface(
        Lon,
        Lat,
        z_plot,
        facecolors=plt.cm.terrain(norm(z)),
        linewidth=0,
        antialiased=False,
        rstride=1,
        cstride=1,
        shade=False,
    )
    ax.set_xlabel("Lon (°E)")
    ax.set_ylabel("Lat (°N)")
    ax.set_zlabel("Elev ×3 if seafloor (m)")
    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    m = plt.cm.ScalarMappable(cmap="terrain", norm=norm)
    m.set_array([])
    fig.colorbar(m, ax=ax, shrink=0.55, label="True elevation / depth (m)")
    fig.tight_layout()
    path = OUT / fname
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("Wrote", path)


def try_kanari_hillshade() -> None:
    """If user downloaded Kanari EEZ zip into this folder, show hillshade GeoTIFF."""
    zips = list(Path(__file__).resolve().parent.glob("*ptj2m36h5d*.zip")) + list(
        OUT.glob("*.zip")
    )
    # also search common download names
    zips += list(Path(__file__).resolve().parent.glob("*Israel*EEZ*.zip"))
    if not zips:
        print("Kanari Israel EEZ zip not found locally — skip high-res hillshade.")
        print("Download: https://data.mendeley.com/datasets/ptj2m36h5d/3")
        return
    zpath = zips[0]
    print("Opening", zpath)
    with zipfile.ZipFile(zpath) as zf:
        tifs = [n for n in zf.namelist() if n.lower().endswith((".tif", ".tiff"))]
        print("GeoTIFFs:", tifs[:10])
        # Prefer hillshade / relief
        prefer = [n for n in tifs if "hill" in n.lower() or "relief" in n.lower() or "shade" in n.lower()]
        target = prefer[0] if prefer else (tifs[0] if tifs else None)
        if not target:
            print("No GeoTIFF in zip.")
            return
        try:
            import rasterio
            from rasterio.io import MemoryFile
        except ImportError:
            print("pip install rasterio to plot Kanari GeoTIFF")
            return
        with zf.open(target) as src_bytes:
            data = src_bytes.read()
        with MemoryFile(data) as mem:
            with mem.open() as ds:
                arr = ds.read(1)
                fig, ax = plt.subplots(figsize=(8, 10))
                im = ax.imshow(arr, cmap="gray")
                ax.set_title(f"Israel EEZ — {Path(target).name}\n(Kanari et al. 2020 multibeam)")
                ax.set_xticks([])
                ax.set_yticks([])
                fig.colorbar(im, ax=ax, fraction=0.046, label="hillshade / value")
                fig.tight_layout()
                path = OUT / "israel_eez_kanari_hillshade.png"
                fig.savefig(path, dpi=160)
                plt.close(fig)
                print("Wrote", path)


def main() -> None:
    Lon, Lat, z = download_etopo()
    print(f"Grid shape {z.shape}; depth range sea: {np.nanmin(z):.0f} … land max {np.nanmax(z):.0f} m")

    plot_map(
        Lon,
        Lat,
        z,
        "Israel & Levantine Sea — ETOPO 2022 bathymetry/topography\n(public NOAA grid; currently available)",
        "01_israel_levantine_etopo_map.png",
    )
    plot_3d(
        Lon,
        Lat,
        z,
        "3D seafloor+land — Israel coast / EEZ region (ETOPO 2022)\nvertical exaggeration ×3 on negative depths",
        "02_israel_sea_etopo_3d.png",
        elev=28,
        azim=-55,
    )
    plot_3d(
        Lon,
        Lat,
        z,
        "Alternate view — looking along the Levantine margin",
        "03_israel_sea_etopo_3d_alt.png",
        elev=22,
        azim=-20,
    )
    try_kanari_hillshade()
    print("\nDone. Open folder:", OUT)


if __name__ == "__main__":
    main()
