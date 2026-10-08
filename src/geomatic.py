"""Geomatic - Tarut Bay baseline: helper functions.

Everything here works on the small clipped inputs shipped in ``data/``:
  * Sentinel-2 L2A clip (16 Aug 2026), 9 bands, 2 channels per file [reflectance, valid-mask]
  * NASA EMIT L2A reflectance clip (12 Aug 2026), 285 bands, WGS84, ~60 m x 54 m pixels
  * NASA EMIT L2A mask clip (same granule), 8 bands
  * 80 visually labelled sample points
No credentials or restricted imagery are used.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import tifffile
from pyproj import Transformer
from scipy import ndimage as ndi
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import cohen_kappa_score

S2_BANDS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11"]
CLASSES = ["dense_tree", "low_veg", "water_bare"]
# NDVI strata used when drawing the sample (area shares computed from the S2 clip, see notebook)
STRATA_RULE = {
    "A_dense_ndvi_gt0.6": "dense_tree",
    "B_mid_ndvi_0.3_0.6": "low_veg",
    "C_low_ndvi_0.1_0.3": "low_veg",
    "D_water_bare_ndvi_le0.1": "water_bare",
}
# EMIT bands used for NDVI. EMIT's wavelength table is not stored in the clip; we assume the
# standard EMIT L2A grid (~7.4 nm, starting ~381 nm): index 38 ~ 665 nm, index 65 ~ 865 nm.
EMIT_RED, EMIT_NIR = 38, 65

_to_utm = Transformer.from_crs(4326, 32639, always_xy=True)
_to_wgs = Transformer.from_crs(32639, 4326, always_xy=True)


# ----------------------------------------------------------------------------- loading
def load_s2(folder):
    """Return dict(band -> reflectance array), valid mask, and (x0, y0, dx, dy) in UTM 39N."""
    folder = Path(folder)
    refl, valid = {}, None
    for b in S2_BANDS:
        a = tifffile.imread(folder / f"{b}.tiff").astype("float64")
        refl[b] = a[..., 0]
        valid = (a[..., 1] > 0) if valid is None else valid & (a[..., 1] > 0)
    tags = tifffile.TiffFile(folder / "B04.tiff").pages[0].tags
    sc, tp = tags["ModelPixelScaleTag"].value, tags["ModelTiepointTag"].value
    return refl, valid, (tp[3], tp[4], sc[0], sc[1])


def load_emit(folder):
    """EMIT reflectance cube (rows, cols, 285), valid mask, mask cube, and geotransform (px, py, ox, oy)."""
    folder = Path(folder)
    cube = tifffile.imread(folder / "emit_rfl_clip.tif").astype("float64")
    mask = tifffile.imread(folder / "emit_mask_clip.tif").astype("float64")
    px, _, _, py, ox, oy = np.loadtxt(folder / "emit_rfl_clip.tfw")  # world file: A, D, B, E, C, F
    return cube, cube[..., 0] > -9000, mask, (px, py, ox, oy)


def ndvi(nir, red):
    return (nir - red) / (nir + red + 1e-9)


# ----------------------------------------------------------------------------- geometry
def emit_rowcol(x, y, gt, shift=(0.0, 0.0)):
    """EMIT (row, col) containing UTM point (x, y). ``shift`` moves the Sentinel-2 point (east, north) in metres."""
    px, py, ox, oy = gt
    lon, lat = _to_wgs.transform(np.asarray(x) + shift[0], np.asarray(y) + shift[1])
    col = np.floor((lon - (ox - px / 2)) / px).astype(int)
    row = np.floor(((oy + py / 2) - lat) / (-py)).astype(int)
    return row, col


def s2_xy(s2geo, shape):
    x0, y0, dx, dy = s2geo
    h, w = shape
    X, Y = np.meshgrid(x0 + (np.arange(w) + 0.5) * dx, y0 - (np.arange(h) + 0.5) * dy)
    return X, Y


def aggregate_ndvi_to_emit(s2, valid, s2geo, emit_valid, gt, shift=(0.0, 0.0)):
    """Mean Sentinel-2 NDVI inside each EMIT pixel (returns mean grid, pixel count grid)."""
    nd = ndvi(s2["B08"], s2["B04"])
    X, Y = s2_xy(s2geo, nd.shape)
    r, c = emit_rowcol(X, Y, gt, shift)
    nr, nc = emit_valid.shape
    ok = valid & (r >= 0) & (r < nr) & (c >= 0) & (c < nc)
    idx = (r * nc + c)[ok]
    n = np.bincount(idx, minlength=nr * nc)
    s = np.bincount(idx, weights=nd[ok], minlength=nr * nc)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (s / n).reshape(nr, nc)
    return mean, n.reshape(nr, nc)


def offset_search(s2, valid, s2geo, emit_ndvi, emit_valid, gt, rng_m=120, step=15, min_px=400, min_n=20):
    """Correlation between EMIT NDVI and aggregated S2 NDVI for a grid of east/north shifts of S2."""
    rows = []
    for e in range(-rng_m, rng_m + 1, step):
        for n in range(-rng_m, rng_m + 1, step):
            m, cnt = aggregate_ndvi_to_emit(s2, valid, s2geo, emit_valid, gt, (e, n))
            k = emit_valid & (cnt >= min_n) & np.isfinite(m)
            if k.sum() < min_px:
                continue
            a, b = emit_ndvi[k], m[k]
            rows.append(dict(shift_east_m=e, shift_north_m=n, n_px=int(k.sum()),
                             r=float(np.corrcoef(a, b)[0, 1]), rmse=float(np.sqrt(((a - b) ** 2).mean()))))
    return pd.DataFrame(rows).sort_values("r", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------------- samples
def point_features(points_geojson, labels_csv, s2, s2geo, cube, emit_valid, gt, shift):
    """One row per labelled point: 3x3-mean S2 bands, EMIT spectrum of the containing pixel (shifted S2 geometry)."""
    feats = {f["properties"]["pid"]: f["geometry"]["coordinates"]
             for f in json.load(open(points_geojson))["features"]}
    lab = pd.read_csv(labels_csv)
    x0, y0, dx, dy = s2geo
    rows = []
    for pid, (lon, lat) in feats.items():
        X, Y = _to_utm.transform(lon, lat)
        cx, cy = int((X - x0) / dx), int((y0 - Y) / dy)
        rec = {"pid": pid, "X": X, "Y": Y}
        for b in S2_BANDS:
            rec["s2_" + b] = s2[b][cy - 1:cy + 2, cx - 1:cx + 2].mean()
        r, c = emit_rowcol(X, Y, gt, shift)
        r, c = int(r), int(c)
        ok = 0 <= r < cube.shape[0] and 0 <= c < cube.shape[1] and emit_valid[r, c]
        rec["emit_ok"] = bool(ok)
        for j in range(cube.shape[2]):
            rec[f"e_{j}"] = cube[r, c, j] if ok else np.nan
        rows.append(rec)
    return pd.DataFrame(rows).merge(lab[["pid", "cls_corrected", "stratum"]], on="pid")


def threshold_class(v, hi=0.6, lo=0.1):
    return np.where(v > hi, "dense_tree", np.where(v > lo, "low_veg", "water_bare"))


def stratum_weights(s2, valid, erode=3):
    """Area share of each NDVI stratum in the S2 clip (valid pixels, eroded to stay inside the polygon)."""
    v = ndi.binary_erosion(valid, iterations=erode)
    nd = ndvi(s2["B08"], s2["B04"])
    masks = {"A_dense_ndvi_gt0.6": nd > 0.6, "B_mid_ndvi_0.3_0.6": (nd > 0.3) & (nd <= 0.6),
             "C_low_ndvi_0.1_0.3": (nd > 0.1) & (nd <= 0.3), "D_water_bare_ndvi_le0.1": nd <= 0.1}
    cnt = {k: int((v & m).sum()) for k, m in masks.items()}
    tot = sum(cnt.values())
    return {k: c / tot for k, c in cnt.items()}, cnt


def point_weights(df, shares):
    n = df.stratum.value_counts()
    w = df.stratum.map(lambda s: shares[s] / n[s]).values
    return w / w.sum()


def spatial_cv(X, y, groups, mask, weights, reps=8, folds=5, trees=100, seed=0):
    """Repeated spatial-block CV (random assignment of 250 m blocks to folds). Returns mean/sd of OA and weighted OA."""
    idx = np.where(mask)[0]
    ug = np.unique(groups)
    oa, woa = [], []
    for r in range(reps):
        rng = np.random.default_rng(seed + r)
        fold_of = {g: rng.integers(0, folds) for g in ug}
        fo = np.array([fold_of[g] for g in groups])
        pred = np.full(len(y), "", dtype=object)
        for k in range(folds):
            te, tr = idx[fo[idx] == k], idx[fo[idx] != k]
            if len(te) == 0 or len(set(y[tr])) < 2:
                continue
            m = RandomForestClassifier(trees, n_jobs=-1, random_state=seed, class_weight="balanced")
            m.fit(X[tr], y[tr])
            pred[te] = m.predict(X[te])
        ok = pred[idx] != ""
        hit = pred[idx][ok] == y[idx][ok]
        oa.append(hit.mean())
        woa.append((weights[idx][ok] * hit).sum() / weights[idx][ok].sum())
    return dict(oa=float(np.mean(oa)), oa_sd=float(np.std(oa)), oa_area_weighted=float(np.mean(woa)), n=int(len(idx)))


def score_fixed(pred, y, mask, weights):
    hit = pred[mask] == y[mask]
    return dict(oa=float(hit.mean()), oa_area_weighted=float((weights[mask] * hit).sum() / weights[mask].sum()),
                kappa=float(cohen_kappa_score(y[mask], pred[mask])), n=int(mask.sum()))
