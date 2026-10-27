# cell_edge_detection

Fiji/ImageJ macro that segments individual periderm cells in zebrafish
maximum-intensity projections (MIPs) and outputs one polygon ROI per complete
cell, shrunk to sit inside the junctional belt. The ROIs are meant as input to
a microridge quantification step (e.g. ridgeQuant), so every ROI is a plain
polygon (`selectionType == 2`).

Current version: **v7** (`periderm_cell_rois_v7.ijm`).

## Requirements

- [Fiji](https://fiji.sc/)
- [MorphoLibJ](https://imagej.net/plugins/morpholibj):
  *Help > Update > Manage update sites > tick **IJPB-plugins** > Apply*, then restart Fiji.
- Input image: single-plane (one channel, one slice, one frame) MIP,
  calibrated in microns with square pixels, saved on disk.
  The macro refuses to run otherwise, because every parameter is in µm.

## Usage

1. Open the MIP in Fiji and bring it to the front.
2. Run `periderm_cell_rois_v7.ijm` (*Plugins > Macros > Run…*). No dialogs.
3. Curate in the ROI Manager: delete bad ROIs, draw missed cells and press **T**
   to add them. ROIs named `auto_NN_check` have an irregular outline (often two
   merged cells) — look at those first.
4. Save the curated set (a companion `save_curated_rois.ijm` step was used for
   this in the original workflow; any ROI Manager save works).

### Outputs (written next to the image)

| File | Contents |
|---|---|
| `<image>_cells_auto.zip` | Detected ROIs, before curation |
| `<image>_cells_log.txt`  | Fate of every candidate cell: kept, edge, vanished, too small, oversized, fragment |
| `<image>_debug_seeds.tif`, `<image>_debug_labels.tif` | Only with `DEBUG = true` |

## How it works

1. **Blur and ratio.** Gaussian blur (0.8 µm) removes ridge texture but keeps
   the junctional belt. Dividing by a broad local background (5 µm blur) removes
   per-cell brightness differences, leaving junction-vs-interior contrast.
2. **Threshold seeds.** Li threshold on the ratio image gives cell interiors;
   these are shrunk 1.5 µm to cut thin bridges through faint junctions.
3. **Basins.** Each cell interior is also a dark basin walled by bright junction.
   Extended minima of the smoothed ratio image are found, with depth set
   *relative* to the image's own junction–interior contrast (0.8×), so dim and
   bright images are treated alike.
4. **Hybrid seeding.** A threshold seed containing ≥2 deep basins (bridged
   cells) is replaced by those basins; a deep basin no seed reached (e.g. a
   bright cell) gets its own seed.
5. **Watershed.** MorphoLibJ marker-controlled watershed on the blurred image
   grows seeds to the junction midline, so cell boundaries don't drift with
   staining intensity. A frame-ring marker absorbs fragments of cells cut by
   the image border.
6. **Edge rule.** A cell is cut by the frame if its interior (not its junction)
   runs along a 1.5 µm frame strip for more than 0.5 µm.
7. **ROIs.** Each label is shrunk 0.6 µm inward (clears the belt), smoothed and
   rebuilt as a polygon. If shrinking splits a cell, the largest piece is kept.
8. **Size filters.** ROIs > 2.5× or < 0.3× the image's median area are
   dropped. ROIs with solidity < 0.88 are flagged `_check`.

## Parameters

All set at the top of the macro; lengths in µm, converted per image.

| Parameter | Default | Meaning |
|---|---|---|
| `BLUR_UM` | 0.8 | Blur that removes ridges, keeps the belt |
| `BG_UM` | 5 | Local background scale for the ratio image |
| `SEED_ERODE_UM` | 1.5 | Threshold-seed shrink (cuts bridges) |
| `MIN_SEED_AREA` | 30 µm² | Minimum seed and final ROI area |
| `BASIN_SMOOTH_UM` | 0.5 | Extra smoothing before basin detection |
| `BASIN_REL` | 0.8 | Basin depth as a fraction of image contrast (0.7–0.9 all work) |
| `BASIN_DEPTH` | 0.4 | Absolute depth, used only if `BASIN_REL = 0` |
| `EDGE_BAND_UM` | 1.5 | Frame strip checked for cell interior |
| `EDGE_INTERIOR_UM` | 0.5 | Interior length along the frame that marks a cut cell |
| `MAX_AREA_FACTOR` | 2.5 | Drop ROIs larger than this × median |
| `MIN_AREA_FACTOR` | 0.3 | Drop ROIs smaller than this × median |
| `ERODE_UM` | 0.6 | Inward shrink from the junction midline |
| `CHECK_SOLIDITY` | 0.88 | Flag ROIs less convex than this |
| `VERTEX_SPACING` | 4 px | Polygon vertex spacing |
| `DEBUG` | false | Save seed and label images |

## Validation

Scored on a curated test set of 35 MIPs (501 hand-curated cells; zebrafish
periderm, 3–4 dpf, control and 0.25% / 0.5% methylcellulose). A curated cell
counts as found if one ROI has IoU ≥ 0.7 with it.

| Version | Found | Missed | False positives |
|---|---|---|---|
| v6.1 | 468 (93.4%) | 33 | 36 |
| **v7** | **491 (98.0%)** | **10** | **47** |

6 of the 10 cells v7 "misses" are detected but outlined differently from the
hand curation. Per-condition recall for v7 is 95–99%.

## Version history

- **v5** — seeds thresholded on blur / local-background ratio (bright cells no longer missed)
- **v6** — seeds shrunk 1.5 µm; frame marker
- **v6.1** — edge rule based on cell interior reaching the frame
- **v7** — hybrid threshold + basin seeding; per-cell log; fragment filter; `_check` flag; keeps largest piece when shrinking splits a cell
