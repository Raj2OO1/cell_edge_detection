

// ---------------- parameters ----------------
var BLUR_UM         = 0.8;   // um. Smooths away ridges, keeps the junction belt.
var BG_UM           = 5;     // um. Local background the blur is divided by.
var SEED_ERODE_UM   = 1.5;   // um. Threshold seeds shrunk by this (cuts thin bridges).
var MIN_SEED_AREA   = 30;    // um^2 floor for a threshold seed and for a final ROI.
var BASIN_SMOOTH_UM = 0.5;   // um. Extra smoothing of the ratio image before basins.
var BASIN_REL       = 0.8;   // a basin counts if it is at least this x the image's
                             // own junction-interior contrast deep. Relative, so
                             // dim and bright images are treated alike.
                             // 0.7-0.9 all beat v6.1 clearly; 0.8 is the centre.
var BASIN_DEPTH     = 0.4;   // used only if BASIN_REL = 0 (absolute ratio units)
var EDGE_BAND_UM     = 1.5;  // um. Strip along the frame checked for cell interior.
var EDGE_INTERIOR_UM = 0.5;  // um. Cell is "cut by the frame" above this.
var MAX_AREA_FACTOR = 2.5;   // drop ROIs > this x the image median area
var MIN_AREA_FACTOR = 0.3;   // drop ROIs < this x the image median area (fragments)
var ERODE_UM        = 0.6;   // um inward from the junction midline (belt excluded)
var CHECK_SOLIDITY  = 0.88;  // ROIs less convex than this get the _check suffix
var VERTEX_SPACING  = 4;     // px between polygon vertices
var THRESH_METHOD   = "Li";
var DEBUG           = false; // true -> also save seed and label images

var logLines;                // globals used inside functions must be declared with var
var saveDir;
var base;

if (nImages == 0) exit("Open an image first.");
orig = getTitle();
base = orig;
if (lastIndexOf(orig, ".") > 0) base = substring(orig, 0, lastIndexOf(orig, "."));
saveDir = getDirectory("image");
if (saveDir == "")
    exit("This image has no folder on disk (unsaved or duplicated).\nSave it first, then re-run.");

getDimensions(W, H, nC, nZ, nT);
if (nC * nZ * nT > 1)
    exit("Expected a single-plane MIP, got " + nC + " channel(s), " + nZ + " slice(s), " + nT + " frame(s).");
getPixelSize(unit, pw, ph);
if (unit == "pixel" || unit == "pixels" || unit == "")
    exit("Image is uncalibrated. Set the scale first (Analyze > Set Scale).");
if (unit == "nm" || unit == "mm" || unit == "cm")
    exit("Image is calibrated in " + unit + ", but every parameter here is in microns.");
if (abs(pw - ph) > 0.01 * pw)
    exit("Pixels are not square (" + pw + " x " + ph + ").");
List.setCommands;
if (List.get("Marker-controlled Watershed") == "" || List.get("Extended Min & Max") == "")
    exit("MorphoLibJ not found.\nHelp > Update > Manage update sites > tick IJPB-plugins > Apply, restart Fiji.");


selectWindow(orig);
run("Select None");
roiManager("reset");
run("ROI Manager...");
setBatchMode(true);
logLines = newArray(0);
erodePx = round(ERODE_UM / pw);
if (erodePx < 1) erodePx = 1;
addLog("image " + orig + "  pixel " + d2s(pw, 4) + " um  macro v7");

selectWindow(orig);
run("Set Measurements...", "area mean median shape redirect=None decimal=4");
run("Clear Results");
setOption("BlackBackground", true);

// ---- 1. blur (ridges gone, belt kept) and local background ----
selectWindow(orig); run("Duplicate...", "title=work");
run("Gaussian Blur...", "sigma=" + (BLUR_UM / pw));
selectWindow(orig); run("Duplicate...", "title=bg");
run("Gaussian Blur...", "sigma=" + (BG_UM / pw));
imageCalculator("Divide create 32-bit", "work", "bg");
rename("ratio");
selectWindow("bg"); close();

// ---- 2. threshold seeds (as v6.1) ----
selectWindow("ratio"); run("Duplicate...", "title=interior");
setAutoThreshold(THRESH_METHOD + " dark");
getThreshold(tLow, tHigh);
run("Convert to Mask");
run("Invert");                                   // 255 = cell interior
run("Duplicate...", "title=liseed");
run("Minimum...", "radius=" + (SEED_ERODE_UM / pw));
run("Analyze Particles...", "size=" + MIN_SEED_AREA + "-Infinity show=[Count Masks] clear");
rename("li");
selectWindow("liseed"); close();
selectWindow("li"); nLi = getValue("Max");

// junction-interior contrast of this image (for BASIN_REL)
selectWindow("ratio");
setThreshold(tLow, 1e30); run("Create Selection"); J = getValue("Median");
run("Make Inverse"); I = getValue("Median");
run("Select None"); resetThreshold();
depth = BASIN_DEPTH;
if (BASIN_REL > 0) depth = BASIN_REL * (J - I);
addLog("Li threshold " + d2s(tLow, 4) + "  contrast " + d2s(J - I, 3) + "  basin depth " + d2s(depth, 3));

// ---- 3. basins: deep regional minima of the smoothed ratio ----
selectWindow("ratio"); run("Duplicate...", "title=rsmooth");
if (BASIN_SMOOTH_UM > 0) run("Gaussian Blur...", "sigma=" + (BASIN_SMOOTH_UM / pw));
// MorphoLibJ reads 'dynamic' as a whole number, so a depth like 0.29 would
// become 0 (every tiny dip a basin -> hundreds of fragments). Scale the image
// x1000 and pass the depth x1000 instead: same result, integer dynamic.
run("Multiply...", "value=1000");
dyn = round(depth * 1000);
if (dyn < 1) dyn = 1;
run("Extended Min & Max", "operation=[Extended Minima] dynamic=" + dyn + " connectivity=8");
rename("basinmask");
run("Connected Components Labeling", "connectivity=8 type=[16 bits]");
rename("bas");
selectWindow("basinmask"); close();
selectWindow("rsmooth"); close();
selectWindow("bas"); nBas = getValue("Max");

// ---- 4. combine: split bridged threshold seeds, add orphan basins ----
newImage("seeds", "16-bit black", W, H, 1);
used = newArray(nBas + 1);
nextId = 1; nSplit = 0; nOrphan = 0;
for (c = 1; c <= nLi; c++) {
    selectWindow("li"); setThreshold(c, c); run("Create Selection");
    if (selectionType() < 0) continue;
    selectWindow("bas"); run("Restore Selection");
    getHistogram(hv, hc, 65536);
    run("Select None");
    inside = newArray(0);
    for (b = 1; b <= nBas; b++) if (hc[b] > 0) { inside = Array.concat(inside, b); used[b] = 1; }
    selectWindow("seeds");
    if (inside.length >= 2) {
        // two or more deep basins in one seed: it bridges cells -> use the basins
        for (k = 0; k < inside.length; k++) {
            selectWindow("bas"); setThreshold(inside[k], inside[k]); run("Create Selection");
            selectWindow("seeds"); run("Restore Selection"); setColor(nextId); fill(); run("Select None");
            nextId++;
        }
        nSplit++;
    } else {
        selectWindow("li"); run("Create Selection");    // threshold c..c still set
        selectWindow("seeds"); run("Restore Selection"); setColor(nextId); fill(); run("Select None");
        nextId++;
    }
}
selectWindow("li"); resetThreshold();
for (b = 1; b <= nBas; b++) {
    if (used[b] == 1) continue;
    // a basin no threshold seed reached (e.g. a bright cell) -> its own seed
    selectWindow("bas"); setThreshold(b, b); run("Create Selection");
    if (selectionType() < 0) continue;
    selectWindow("seeds"); run("Restore Selection"); setColor(nextId); fill(); run("Select None");
    nextId++; nOrphan++;
}
selectWindow("bas"); resetThreshold();
addLog("threshold seeds " + nLi + "  basins " + nBas + "  split seeds " + nSplit + "  orphan basins added " + nOrphan);
if (nBas > 4 * nLi + 20)
    addLog("WARNING: far more basins than seeds - basin depth may be wrong (dynamic used " + dyn + ")");

// frame ring = one more marker; MorphoLibJ relabels touching markers (binary)
selectWindow("seeds");
setColor(nextId);
drawLine(0, 0, W - 1, 0); drawLine(0, H - 1, W - 1, H - 1);
drawLine(0, 0, 0, H - 1); drawLine(W - 1, 0, W - 1, H - 1);
run("Select None");

// ---- 5. watershed: grow seeds to the junction midline ----
nBefore = nImages;
run("Marker-controlled Watershed", "input=work marker=seeds mask=None binary calculate use");
if (nImages == nBefore) { cleanup(); exit("Watershed produced no output."); }
rename("labels");
selectWindow("labels"); nLab = getValue("Max");
frameLab = getPixel(0, 0);
if (frameLab == 0) frameLab = getPixel(1, 1);

if (DEBUG) {
    selectWindow("seeds");  saveAs("Tiff", saveDir + base + "_debug_seeds.tif");  rename("seeds");
    selectWindow("labels"); saveAs("Tiff", saveDir + base + "_debug_labels.tif"); rename("labels");
}

// ---- 6. frame contact per label, from one histogram (no selections) ----
bandPx = round(EDGE_BAND_UM / pw); if (bandPx < 1) bandPx = 1;
newImage("band", "8-bit black", W, H, 1);
setColor(255);
fillRect(0, 0, W, bandPx); fillRect(0, H - bandPx, W, bandPx);
fillRect(0, 0, bandPx, H); fillRect(W - bandPx, 0, bandPx, H);
imageCalculator("AND", "band", "interior");      // frame strip AND cell interior
selectWindow("band");
run("Divide...", "value=255");                   // 0/1
imageCalculator("Multiply create 32-bit", "labels", "band");
rename("bandlab");
getHistogram(bvals, contactPx, nLab + 1, -0.5, nLab + 0.5);   // bin v = label v
selectWindow("bandlab"); close();

// ---- 7. labels -> eroded polygon ROIs, logging every decision ----
areas = newArray(0); sols = newArray(0); labOf = newArray(0);
nEdge = 0;
for (v = 1; v <= nLab; v++) {
    if (v == frameLab) continue;
    selectWindow("labels");
    setThreshold(v, v);
    run("Create Selection");
    if (selectionType() < 0) continue;
    labArea = getValue("Area");
    contact = contactPx[v] / bandPx * pw;
    if (contact > EDGE_INTERIOR_UM) {
        nEdge++;
        addLog("label " + v + "  area " + d2s(labArea, 1) + "  -> EDGE (interior along frame " + d2s(contact, 2) + " um)");
        continue;
    }
    selectWindow(orig);
    run("Restore Selection");
    run("Enlarge...", "enlarge=-" + erodePx + " pixel");
    if (selectionType() < 0) { addLog("label " + v + "  area " + d2s(labArea, 1) + "  -> VANISHED when shrunk"); continue; }
    note = "";
    if (selectionType() == 9) {
        // shrinking left several pieces: keep the largest instead of losing the cell
        n0 = roiManager("count");
        roiManager("Add"); roiManager("Select", n0); roiManager("Split");
        n1 = roiManager("count"); best = -1; bestA = -1;
        for (q = n0 + 1; q < n1; q++) {
            roiManager("Select", q); aq = getValue("Area");
            if (aq > bestA) { bestA = aq; best = q; }
        }
        roiManager("Select", best); getSelectionCoordinates(px, py);
        for (q = n1 - 1; q >= n0; q--) { roiManager("Select", q); roiManager("Delete"); }
        selectWindow(orig); makeSelection("freehand", px, py);
        note = "  (kept largest of " + (n1 - n0 - 1) + " pieces)";
    }
    // ridgeQuant rejects selectionType >= 4: rebuild as a true polygon
    run("Interpolate", "interval=" + VERTEX_SPACING + " smooth");
    getSelectionCoordinates(vx, vy);
    if (lengthOf(vx) < 3) { addLog("label " + v + "  -> DEGENERATE outline"); continue; }
    makeSelection("polygon", vx, vy);
    if (selectionType() != 2) { addLog("label " + v + "  -> NOT A POLYGON"); continue; }
    a0 = getValue("Area");
    if (a0 < MIN_SEED_AREA) { addLog("label " + v + "  area " + d2s(a0, 1) + "  -> TOO SMALL (< " + MIN_SEED_AREA + " um2)"); continue; }
    List.setMeasurements; sol = List.getValue("Solidity");
    Roi.setProperty("image", orig);
    roiManager("Add");
    areas = Array.concat(areas, a0); sols = Array.concat(sols, sol); labOf = Array.concat(labOf, v);
    addLog("label " + v + "  area " + d2s(a0, 1) + "  solidity " + d2s(sol, 3) + "  -> kept" + note);
}
selectWindow("labels"); resetThreshold();
cleanup();

n = roiManager("count");
addLog("ROIs in ROI Manager after the per-cell loop: " + n);
if (n == 0) { saveLog(); setBatchMode(false); exit("No ROIs survived. See " + base + "_cells_log.txt"); }

// ---- 8. size filters relative to this image's median cell ----
med = median(areas);
nBig = 0;
for (r = n - 1; r >= 0; r--) {
    if (areas[r] > MAX_AREA_FACTOR * med) {
        addLog("label " + labOf[r] + "  area " + d2s(areas[r], 1) + "  -> REMOVED, > " + MAX_AREA_FACTOR + " x median " + d2s(med, 1));
        roiManager("Select", r); roiManager("Delete");
        areas = Array.deleteIndex(areas, r); sols = Array.deleteIndex(sols, r); labOf = Array.deleteIndex(labOf, r);
        nBig++;
    }
}
n = roiManager("count");
med2 = median(areas);
nSmall = 0;
for (r = n - 1; r >= 0; r--) {
    if (areas[r] < MIN_AREA_FACTOR * med2) {
        addLog("label " + labOf[r] + "  area " + d2s(areas[r], 1) + "  -> REMOVED, < " + MIN_AREA_FACTOR + " x median " + d2s(med2, 1));
        roiManager("Select", r); roiManager("Delete");
        areas = Array.deleteIndex(areas, r); sols = Array.deleteIndex(sols, r); labOf = Array.deleteIndex(labOf, r);
        nSmall++;
    }
}

// ---- 9. name, save, hand over ----
n = roiManager("count");
nCheck = 0;
for (r = 0; r < n; r++) {
    roiManager("Select", r);
    nm = "auto_" + IJ.pad(r + 1, 2);
    if (sols[r] < CHECK_SOLIDITY) { nm = nm + "_check"; nCheck++; }
    roiManager("Rename", nm);
}
roiManager("Deselect");
if (n > 0) roiManager("Save", saveDir + base + "_cells_auto.zip");
addLog("SUMMARY  kept " + n + "  edge " + nEdge + "  oversized " + nBig + "  fragments " + nSmall + "  flagged _check " + nCheck);
saveLog();

setBatchMode(false);
// safety net: if the ROI Manager came back empty, reload what was saved
if (roiManager("count") == 0 && n > 0) roiManager("Open", saveDir + base + "_cells_auto.zip");
selectWindow(orig);
run("Select None");
roiManager("Show All with labels");
print(base + ": " + n + " ROIs (" + nCheck + " flagged _check). Curate, then run save_curated_rois.ijm.");

// ---------------- helpers ----------------
function addLog(s) { logLines = Array.concat(logLines, s); }
function saveLog() {
    f = File.open(saveDir + base + "_cells_log.txt");
    for (i = 0; i < logLines.length; i++) print(f, logLines[i]);
    File.close(f);
}
function median(a) {
    b = Array.copy(a); Array.sort(b); m = b.length;
    if (m == 0) return 0;
    if (m % 2 == 1) return b[(m - 1) / 2];
    return 0.5 * (b[m / 2 - 1] + b[m / 2]);
}
function cleanup() {
    names = newArray("work", "ratio", "interior", "li", "bas", "seeds", "labels", "band", "bandlab", "bg", "rsmooth", "basinmask", "liseed");
    for (k = 0; k < names.length; k++) if (isOpen(names[k])) { selectWindow(names[k]); close(); }
}
