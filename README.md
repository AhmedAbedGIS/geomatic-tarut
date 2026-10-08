# Geomatic – Mangrove and coastal vegetation baseline, Tarut Bay

**Arab Youth Space Hackathon 2026 · Challenge 813 – Ecosystem Health, Biodiversity & Blue Carbon · Proof of concept**

> **Working title:** *Geomatic: A satellite-first baseline for mangrove health and a blue-carbon scenario in Tarut Bay, using NASA EMIT hyperspectral and Sentinel-2 data.*
> (The largest dense patches were checked as mangrove on Esri high-resolution imagery by the analyst; the whole dense-tree class is not verified pixel by pixel – see Limitations.)

**الملخص بالعربية:** نقارن في خليج تاروت (≈233 هكتاراً) بين بيانات NASA EMIT الطيفية (285 حزمة) وSentinel-2 (9 حزم) في فصل الغطاء النباتي الكثيف عن النبات المنخفض والماء/الأرض الجرداء، اعتماداً على 80 نقطة مرجعية مسمّاة بصرياً. لم تتفوق EMIT على Sentinel-2 في هذه العينة الصغيرة (≈0.65 لكليهما، والفرق غير معنوي)، ويحتاج استخدام EMIT إلى معايرة مكانية دقيقة. كل النتائج تُعاد بتشغيل دفتر Jupyter. لا تُستخدم بيانات درون، ولا يُستنتج أي تدهور من مقارنة سنتين، وقد تحقق المحلل بصرياً من أن أكبر البقع الكثيفة مانغروف على صور Esri عالية الدقة، دون تحقق ميداني أو تأكيد لكل بكسل.

## What this repository does
1. Reads clipped **NASA EMIT L2A** reflectance (285 bands, 60 m) + mask and **Sentinel-2 L2A** (9 bands) for one study polygon.
2. Co-registers the two sensors (shift search, ±120 m) and reports EMIT quality flags (AOD ≈ 0.5, no cloud flagged).
3. Compares estimators of *dense tree / low vegetation / water-bare* against **80 visually labelled points** with repeated spatial-block cross-validation.
4. Derives the dense-tree area share (bootstrap interval) and a clearly labelled **soil-carbon scenario** from a published stock value.

Everything runs from `notebooks/01_tarut_baseline.ipynb` (executed outputs are stored in the notebook and in `outputs/`).

## Results (n = 80 points, approx. 95 % interval ±0.10)
| Estimator | Overall accuracy | Area-weighted |
|---|---|---|
| Sentinel-2 NDVI threshold (0.1 / 0.6) | 0.61 | 0.65 |
| Sentinel-2 random forest, red + NIR | 0.63 | 0.66 |
| Sentinel-2 random forest, 9 bands | 0.66 | 0.69 |
| EMIT random forest, 285 bands (co-registered) | 0.65 | 0.66 |
| EMIT random forest, 285 bands (no shift) | 0.51 | 0.53 |
| EMIT NDVI with Sentinel-2 thresholds | 0.30 | 0.36 |

* Hyperspectral EMIT did **not** beat multispectral Sentinel-2 on this reference; differences are within noise.
* Co-registration matters: EMIT↔S2 NDVI correlation rises from r = 0.73 to 0.85 with a ~(+30 m E, −60 m N) shift of the S2 grid, but RMSE stays ≈ 0.35.
* Dense-tree share ≈ 27 % (95 % bootstrap 21–33 %) of ≈ 237 ha ≈ 64 ha.
* Soil-carbon *scenario* (not a measurement): 64 ha × 43 ± 5 Mg C ha⁻¹ ≈ 2.8 kt C (range 1.9–3.7 kt), assuming the dense-tree area is mangrove with Red Sea soil stocks (verified for the largest patches only).

## Data (included in `data/`, small clips)
| Folder | Content | Source / identifier | Date |
|---|---|---|---|
| `emit_20260812/` | EMIT L2A reflectance clip (285 bands, float32, WGS84) + mask clip (8 bands) | NASA EMIT, LP DAAC – granule `EMIT_L2A_RFL_001_20260812T110954_2622407_050` (mask: `EMIT_L2A_MASK_001_20260812T110954_2622407_050`) | 2026-08-12 11:09:54 UTC |
| `sentinel2_20260816/` | Sentinel-2 L2A bands B02–B08, B8A, B11 (reflectance + valid mask) | Copernicus Sentinel-2 L2A via Copernicus Browser, tile expected **T39RVK** (MGRS square computed from the study-area centroid; full product name to be confirmed from Copernicus Browser search results), downloaded as a clip through the Copernicus Browser Analytical (Process API) download | 2026-08-16 |
| `reference/` | 80 sample points, raw and corrected visual labels, candidate dense-vegetation patches, the offline labelling tool | produced by the team | 2026-10-08/09 |
| `study_area.geojson` | Study polygon (WGS84, ≈233 ha) | drawn by the team | – |

Satellite data: NASA EMIT (public, LP DAAC) and Copernicus Sentinel-2 (free and open). The Arab satellite 813 data were not available for this area at submission time (confirmed by the organisers); EMIT and Sentinel-2 were accepted by the organisers as the satellite basis. **No drone data are used.**

## How to run
```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_tarut_baseline.ipynb
```
Tested with Python 3.13; runtime ≈ 1–2 minutes. Input → output example: `data/` → `outputs/` (`comparison_results.csv`, `comparison_accuracy.png`, `offset_search.png`, `emit_mask_summary.csv`, `carbon_scenario.csv`, `metrics.json`).

## Limitations (please read)
* **Reference:** 80 points labelled by **one analyst** from 10 m Sentinel-2 imagery (true colour, NIR-R-G, NDVI); 14 points first called *dense tree* were corrected to *water* after the NIR-R-G / NDVI check (raw labels kept). No inter-analyst agreement, no field or drone data. Accuracy intervals are wide (±0.10).
* **Mangrove identity:** the analyst visually verified the largest candidate dense patches as mangrove on Esri high-resolution imagery. This is a single-analyst visual check, with no field or drone confirmation, and it does not cover every pixel of the *dense tree* class; species cannot be confirmed at 10 m.
* **Dates differ** (EMIT 12 Aug 11:09 UTC, Sentinel-2 16 Aug ≈ 07:00 UTC) – tide, sun angle and the 4-day gap are uncontrolled. EMIT aerosol optical depth ≈ 0.5 (dust/haze).
* **EMIT geolocation** is imperfect; the shift was estimated from the same pair of scenes (a calibration, not an independent test). EMIT wavelengths are assumed to follow the standard EMIT grid (band 38 ≈ 665 nm, band 65 ≈ 865 nm).
* **No change/degradation claim**: only one date per sensor is used.
* **Spectral unmixing and EMIT-specific stress indicators were not validated** – no areal reference or stress reference exists. They remain future work.
* **Carbon** is a scenario using soil organic carbon stock of Red Sea mangroves (43 ± 5 Mg C ha⁻¹, top 1 m; Almahasheer et al. 2017, *Sci. Rep.* 7:9700, doi:10.1038/s41598-017-10424-9). Biomass carbon was not estimated (no height or field data).

## Credits and licences
Code: © the team (add a licence file before publishing, e.g. MIT). NASA EMIT data: NASA/JPL/LP DAAC. Contains modified Copernicus Sentinel data (2026). Challenge notebook (Tanager/Planet example) by Dr. Vincent Markiet / Space42 was consulted for the challenge framing; no Planet data are used here.
