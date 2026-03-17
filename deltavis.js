"use strict";

// ---------------------------------------------------------------------------
// Delta file parser
//
// Parses the .delta format produced by nucmer/promer (MUMmer).
// Returns a DeltaFile object:
//   { referencePath, queryPath, format, alignmentSections: [...] }
//
// Each AlignmentSection:
//   { refId, queryId, refLen, queryLen, alignments: [...] }
//
// Each Alignment:
//   { refStart, refEnd, queryStart, queryEnd, errors, simErrors, stopCodons }
// ---------------------------------------------------------------------------

function parseDelta(input) {
    const lines = input.split(/\r?\n/);
    let lineNum = 0;
    let referencePath = "";
    let queryPath = "";
    let format = "";
    const sections = [];
    let currentSection = null;
    let inDistances = false;

    for (const rawLine of lines) {
        lineNum++;
        const line = rawLine.trim();

        if (line === "") {
            continue;
        }

        if (lineNum === 1) {
            // First non-blank logical line: two file paths
            const parts = line.split(/\s+/);
            if (parts.length < 2) {
                throw new Error(`Line ${lineNum}: expected two file paths`);
            }
            referencePath = parts[0];
            queryPath = parts[parts.length - 1];

        } else if (format === "") {
            // Second non-blank logical line: format identifier
            if (line === "NUCMER" || line === "PROMER") {
                format = line;
            } else {
                throw new Error(`Line ${lineNum}: expected NUCMER or PROMER, got '${line}'`);
            }

        } else if (line.startsWith(">")) {
            // Section header: >refId queryId refLen queryLen
            if (currentSection) {
                currentSection.alignments.reverse();
                sections.push(currentSection);
            }
            const parts = line.slice(1).trim().split(/\s+/);
            if (parts.length !== 4) {
                throw new Error(
                    `Line ${lineNum}: invalid section header, expected >refId queryId refLen queryLen`
                );
            }
            const refLen = parseInt(parts[2], 10);
            const queryLen = parseInt(parts[3], 10);
            if (isNaN(refLen) || isNaN(queryLen)) {
                throw new Error(`Line ${lineNum}: invalid lengths in section header`);
            }
            currentSection = {
                refId: parts[0],
                queryId: parts[1],
                refLen,
                queryLen,
                alignments: [],
            };
            inDistances = false;

        } else if (inDistances) {
            // Distance value (single integer per line); 0 terminates the block
            const val = parseInt(line, 10);
            if (isNaN(val)) {
                throw new Error(
                    `Line ${lineNum}: expected integer distance value, got '${line}'`
                );
            }
            if (val === 0) {
                inDistances = false;
            }

        } else {
            // Either an alignment record (7 integers) or a distance value
            const parts = line.split(/\s+/);
            if (parts.length === 7) {
                const nums = parts.map((s) => parseInt(s, 10));
                if (nums.some(isNaN)) {
                    throw new Error(
                        `Line ${lineNum}: could not parse integers in alignment record`
                    );
                }
                if (!currentSection) {
                    throw new Error(
                        `Line ${lineNum}: alignment record without section header`
                    );
                }
                currentSection.alignments.push({
                    refStart: nums[0],
                    refEnd: nums[1],
                    queryStart: nums[2],
                    queryEnd: nums[3],
                    errors: nums[4],
                    simErrors: nums[5],
                    stopCodons: nums[6],
                });
                inDistances = true;
            } else if (parts.length === 1) {
                const val = parseInt(line, 10);
                if (isNaN(val)) {
                    throw new Error(
                        `Line ${lineNum}: expected integer distance value, got '${line}'`
                    );
                }
                if (val === 0) {
                    inDistances = false;
                }
            } else {
                throw new Error(
                    `Line ${lineNum}: unexpected format (got ${parts.length} fields)`
                );
            }
        }
    }

    // Flush last section
    if (currentSection) {
        currentSection.alignments.reverse();
        sections.push(currentSection);
    }

    if (!referencePath) {
        throw new Error("No file paths found (empty file?)");
    }
    if (!format) {
        throw new Error("No format line found");
    }

    return { referencePath, queryPath, format, alignmentSections: sections };
}


// ---------------------------------------------------------------------------
// Dot-plot SVG renderer
//
// Takes a DeltaFile and returns an SVG element showing the alignment dot plot.
// ---------------------------------------------------------------------------

const SVG_NS = "http://www.w3.org/2000/svg";

// Build a PLOT dimensions object.  dataSize controls the width and height
// of the data area (the square where alignments are drawn).
function makePlot(dataSize) {
    const p = {
        dataWidth: dataSize,
        dataHeight: dataSize,
        marginLeft: 100,
        marginBottom: 60,
        marginTop: 20,
        marginRight: 20,
    };
    p.totalWidth = p.marginLeft + p.dataWidth + p.marginRight;
    p.totalHeight = p.marginTop + p.dataHeight + p.marginBottom;
    return p;
}

const PLOT_DEFAULT_SIZE = 700;

// Collect unique sequences in order of first appearance.
function collectUnique(sections, idKey, lenKey) {
    const seen = new Set();
    const result = [];
    for (const sec of sections) {
        const id = sec[idKey];
        if (!seen.has(id)) {
            seen.add(id);
            result.push({ name: id, len: sec[lenKey] });
        }
    }
    return result;
}

// Compute cumulative offsets for a list of {name, len} entries.
function computeOffsets(seqs) {
    const offsets = new Map();
    let total = 0;
    for (const seq of seqs) {
        offsets.set(seq.name, total);
        total += seq.len;
    }
    return { offsets, total };
}

function computeLayout(delta, P) {
    const refs = collectUnique(delta.alignmentSections, "refId", "refLen");
    const queries = collectUnique(delta.alignmentSections, "queryId", "queryLen");
    const { offsets: refOffsets, total: totalRefLen } = computeOffsets(refs);
    const { offsets: queryOffsets, total: totalQueryLen } = computeOffsets(queries);
    return {
        refs,
        queries,
        refOffsets,
        queryOffsets,
        totalRefLen,
        totalQueryLen,
        scaleX: totalRefLen > 0 ? P.dataWidth / totalRefLen : 1,
        scaleY: totalQueryLen > 0 ? P.dataHeight / totalQueryLen : 1,
    };
}

// Coordinate transforms
function toSvgY(P, layout, seqName, pos) {
    const offset = layout.queryOffsets.get(seqName) || 0;
    return P.marginTop + P.dataHeight - (offset + pos) * layout.scaleY;
}

// SVG element helpers
function svgEl(tag, attrs, children) {
    const el = document.createElementNS(SVG_NS, tag);
    for (const [k, v] of Object.entries(attrs)) {
        el.setAttribute(k, String(v));
    }
    if (children) {
        for (const child of children) {
            el.appendChild(child);
        }
    }
    return el;
}

function svgText(content) {
    return document.createTextNode(content);
}

// Choose a "nice" tick interval for ~8 ticks across the given range.
function niceInterval(range) {
    const rough = range / 8;
    const magnitude = Math.pow(10, Math.floor(Math.log10(rough)));
    const normalized = rough / magnitude;
    if (normalized <= 1) return magnitude;
    if (normalized <= 2) return 2 * magnitude;
    if (normalized <= 5) return 5 * magnitude;
    return 10 * magnitude;
}

function generateTicks(start, end, interval) {
    if (interval <= 0) return [];
    const ticks = [];
    let t = Math.ceil(start / interval) * interval;
    while (t <= end + interval * 0.001) {
        ticks.push(t);
        t += interval;
    }
    return ticks;
}

function formatTickLabel(val) {
    if (val >= 1_000_000) {
        return (Math.round(val / 100_000) / 10) + "M";
    }
    if (val >= 1000) {
        return (Math.round(val / 100) / 10) + "k";
    }
    return String(Math.round(val));
}

// Render the dot plot SVG.
// rotation: fraction of the total reference length to shift the x-axis
//           origin (0–1).  Treats the reference as circular so coordinates
//           wrap around.  Tick labels always show original positions.
// plotSize: side length (in SVG units) of the square data area.
function renderDotPlot(delta, rotation, plotSize) {
    rotation = rotation || 0;
    plotSize = plotSize || PLOT_DEFAULT_SIZE;
    const P = makePlot(plotSize);
    const layout = computeLayout(delta, P);
    const L = layout.totalRefLen;
    const rotBp = rotation * L;

    // Shift a global reference position by the rotation offset.
    // When rotBp is 0 we skip the modulo so that position L maps to the
    // right edge rather than wrapping to 0.
    function shiftRef(globalPos) {
        if (rotBp === 0) return globalPos;
        return ((globalPos - rotBp) % L + L) % L;
    }

    // Global position of the wrap boundary (genome coordinate that maps
    // to display-x = 0).
    const boundaryGlobal = ((rotBp % L) + L) % L;

    const viewBox = `0 0 ${P.totalWidth} ${P.totalHeight}`;
    const svg = svgEl("svg", {
        viewBox,
        width: "100%",
        style: "max-width: " + (P.totalWidth + 20) + "px",
    });

    // Background
    svg.appendChild(
        svgEl("rect", {
            x: 0, y: 0,
            width: P.totalWidth, height: P.totalHeight,
            fill: "white",
        })
    );

    // Data area border
    svg.appendChild(
        svgEl("rect", {
            x: P.marginLeft, y: P.marginTop,
            width: P.dataWidth, height: P.dataHeight,
            fill: "white", stroke: "#ccc", "stroke-width": 1,
        })
    );

    // Alignment lines (with rotation-aware wrapping)
    const alignG = svgEl("g", {});
    for (const section of delta.alignmentSections) {
        const refOffset = layout.refOffsets.get(section.refId) || 0;
        for (const a of section.alignments) {
            const globalStart = refOffset + a.refStart;
            const globalEnd = refOffset + a.refEnd;
            const isForward = a.queryStart <= a.queryEnd;
            const color = isForward ? "#00BFFF" : "#9933FF";
            const lineAttrs = { stroke: color, "stroke-width": 2, "stroke-linecap": "round" };

            const y1 = toSvgY(P, layout, section.queryId, a.queryStart);
            const y2 = toSvgY(P, layout, section.queryId, a.queryEnd);

            const sStart = shiftRef(globalStart);
            const sEnd = shiftRef(globalEnd);

            // Tooltip with original genomic coordinates
            const tip = section.refId + ":" + a.refStart + "–" + a.refEnd
                + " vs " + section.queryId + ":" + a.queryStart + "–" + a.queryEnd
                + " (" + (isForward ? "fwd" : "rev") + ")";
            const titleEl = svgEl("title", {}, [svgText(tip)]);

            if (rotBp === 0 || sStart <= sEnd) {
                // No wrapping — single line
                const line = svgEl("line", {
                    x1: P.marginLeft + sStart * layout.scaleX, y1,
                    x2: P.marginLeft + sEnd * layout.scaleX, y2,
                    ...lineAttrs,
                });
                line.appendChild(titleEl);
                alignG.appendChild(line);
            } else {
                // Alignment crosses the wrap boundary — split into two segments.
                // Wrap both in a <g> so the tooltip covers both halves.
                const refSpan = globalEnd - globalStart;
                const frac = refSpan > 0 ? (boundaryGlobal - globalStart) / refSpan : 0;
                const qMid = a.queryStart + frac * (a.queryEnd - a.queryStart);
                const yMid = toSvgY(P, layout, section.queryId, qMid);

                const g = svgEl("g", {});
                g.appendChild(titleEl);
                // Right segment: sStart → right edge
                g.appendChild(svgEl("line", {
                    x1: P.marginLeft + sStart * layout.scaleX, y1,
                    x2: P.marginLeft + P.dataWidth, y2: yMid,
                    ...lineAttrs,
                }));
                // Left segment: left edge → sEnd
                g.appendChild(svgEl("line", {
                    x1: P.marginLeft, y1: yMid,
                    x2: P.marginLeft + sEnd * layout.scaleX, y2,
                    ...lineAttrs,
                }));
                alignG.appendChild(g);
            }
        }
    }
    svg.appendChild(alignG);

    // X-axis ticks — positions are in original genome coordinates,
    // placed at their rotated display positions.
    const tickInterval = niceInterval(L);
    const ticks = generateTicks(0, L, tickInterval);
    const axisY = P.marginTop + P.dataHeight;
    const xAxisG = svgEl("g", {});
    const seenX = new Set();
    for (const tick of ticks) {
        const shifted = shiftRef(tick);
        const x = P.marginLeft + shifted * layout.scaleX;
        // Deduplicate ticks that map to the same display position
        // (e.g. tick 0 and tick L in circular view).
        const xKey = Math.round(x * 10);
        if (seenX.has(xKey)) continue;
        seenX.add(xKey);
        xAxisG.appendChild(
            svgEl("line", {
                x1: x, y1: axisY, x2: x, y2: axisY + 6,
                stroke: "#333", "stroke-width": 1,
            })
        );
        xAxisG.appendChild(
            svgEl("text", {
                x, y: axisY + 20,
                "text-anchor": "middle", "font-size": 11, fill: "#333",
            }, [svgText(formatTickLabel(tick))])
        );
    }
    svg.appendChild(xAxisG);

    // Y-axis labels (one per query sequence, at its midpoint)
    const yAxisG = svgEl("g", {});
    for (const seq of layout.queries) {
        const offset = layout.queryOffsets.get(seq.name) || 0;
        const midpoint = offset + seq.len / 2;
        const y = P.marginTop + P.dataHeight - midpoint * layout.scaleY;
        yAxisG.appendChild(
            svgEl("text", {
                x: P.marginLeft - 8, y,
                "text-anchor": "end", "dominant-baseline": "central",
                "font-size": 11, fill: "#333",
            }, [svgText(seq.name)])
        );
    }
    svg.appendChild(yAxisG);

    // Horizontal dashed separators between query sequences
    const sepG = svgEl("g", {});
    for (const seq of layout.queries) {
        const offset = layout.queryOffsets.get(seq.name) || 0;
        if (offset > 0) {
            const y = P.marginTop + P.dataHeight - offset * layout.scaleY;
            sepG.appendChild(
                svgEl("line", {
                    x1: P.marginLeft, y1: y,
                    x2: P.marginLeft + P.dataWidth, y2: y,
                    stroke: "#ccc", "stroke-width": 1,
                    "stroke-dasharray": "4,4",
                })
            );
        }
    }
    svg.appendChild(sepG);

    // X-axis label (reference name)
    const refName =
        delta.alignmentSections.length > 0
            ? delta.alignmentSections[0].refId
            : "Reference";
    const xLabelY = P.marginTop + P.dataHeight + 45;
    svg.appendChild(
        svgEl("text", {
            x: P.marginLeft + P.dataWidth / 2, y: xLabelY,
            "text-anchor": "middle", "font-size": 13, fill: "#333",
        }, [svgText(refName)])
    );

    return svg;
}


// ---------------------------------------------------------------------------
// App UI – drag-and-drop zone, file reading, error display
// ---------------------------------------------------------------------------

function initDeltaVis(container) {
    let state = "waiting"; // "waiting" | "error" | "displaying"

    render();

    function render() {
        container.innerHTML = "";
        if (state === "waiting") {
            container.appendChild(buildDropZone());
        } else if (state === "error") {
            container.appendChild(buildError(state.message));
        } else if (state === "displaying") {
            container.appendChild(buildPlot(state.delta));
        }
    }

    function showPlot(delta) {
        state = "displaying";
        container.innerHTML = "";
        container.appendChild(buildPlot(delta));
    }

    function showError(message) {
        state = "error";
        container.innerHTML = "";
        container.appendChild(buildError(message));
    }

    function showDropZone() {
        state = "waiting";
        container.innerHTML = "";
        container.appendChild(buildDropZone());
    }

    function handleFileContent(content) {
        try {
            const delta = parseDelta(content);
            showPlot(delta);
        } catch (err) {
            showError(err.message);
        }
    }

    function readFile(file) {
        const reader = new FileReader();
        reader.onload = () => handleFileContent(reader.result);
        reader.onerror = () => showError("Failed to read file");
        reader.readAsText(file);
    }

    // --- Drop zone ---

    function buildDropZone() {
        const zone = document.createElement("div");
        Object.assign(zone.style, {
            border: "3px dashed #aaa",
            borderRadius: "12px",
            padding: "4rem 2rem",
            textAlign: "center",
            background: "#fafafa",
            cursor: "pointer",
            transition: "all 0.15s ease",
        });

        const prompt = document.createElement("p");
        Object.assign(prompt.style, { fontSize: "1.2rem", color: "#666", marginBottom: "0.5rem" });
        prompt.textContent = "Drop a .delta file here";

        const sub = document.createElement("p");
        Object.assign(sub.style, { fontSize: "0.9rem", color: "#999" });
        sub.textContent = "or click to browse";

        zone.append(prompt, sub);

        function setHover(hovering) {
            zone.style.border = hovering ? "3px solid #00BFFF" : "3px dashed #aaa";
            zone.style.background = hovering ? "#f0faff" : "#fafafa";
            prompt.textContent = hovering ? "Drop it!" : "Drop a .delta file here";
        }

        zone.addEventListener("dragenter", (e) => { e.preventDefault(); setHover(true); });
        zone.addEventListener("dragover", (e) => { e.preventDefault(); setHover(true); });
        zone.addEventListener("dragleave", () => setHover(false));
        zone.addEventListener("drop", (e) => {
            e.preventDefault();
            setHover(false);
            if (e.dataTransfer.files.length > 0) {
                readFile(e.dataTransfer.files[0]);
            }
        });
        zone.addEventListener("click", () => {
            const input = document.createElement("input");
            input.type = "file";
            input.addEventListener("change", () => {
                if (input.files.length > 0) {
                    readFile(input.files[0]);
                }
            });
            input.click();
        });

        return zone;
    }

    // --- Error view ---

    function buildError(message) {
        const wrapper = document.createElement("div");

        const box = document.createElement("div");
        Object.assign(box.style, {
            background: "#fff0f0",
            border: "1px solid #ffcccc",
            borderRadius: "8px",
            padding: "1.5rem",
            marginBottom: "1rem",
        });

        const title = document.createElement("p");
        Object.assign(title.style, { color: "#cc0000", fontWeight: "bold", marginBottom: "0.5rem" });
        title.textContent = "Parse Error";

        const detail = document.createElement("p");
        Object.assign(detail.style, { color: "#660000", fontFamily: "monospace", fontSize: "0.9rem", margin: "0" });
        detail.textContent = message;

        box.append(title, detail);
        wrapper.appendChild(box);
        wrapper.appendChild(buildLoadAnotherButton());

        return wrapper;
    }

    // --- Plot view ---

    function buildPlot(delta) {
        const wrapper = document.createElement("div");

        const totalAlignments = delta.alignmentSections.reduce(
            (sum, s) => sum + s.alignments.length, 0
        );
        const info = document.createElement("p");
        Object.assign(info.style, { color: "#666", fontSize: "0.85rem", marginBottom: "0.5rem" });
        info.textContent =
            `${delta.alignmentSections.length} section(s), ${totalAlignments} alignment(s)`;
        wrapper.appendChild(info);

        // Rotation control
        const rotRow = document.createElement("div");
        rotRow.className = "form-group mb-2";
        Object.assign(rotRow.style, { display: "flex", alignItems: "center", gap: "0.75rem" });

        const rotLabel = document.createElement("label");
        rotLabel.textContent = "Rotate reference";
        Object.assign(rotLabel.style, { margin: "0", whiteSpace: "nowrap", fontSize: "0.85rem" });

        const slider = document.createElement("input");
        slider.type = "range";
        slider.className = "custom-range";
        slider.min = "0";
        slider.max = "100";
        slider.value = "0";
        slider.step = "1";
        slider.style.flex = "1";

        const pctLabel = document.createElement("span");
        pctLabel.textContent = "0%";
        Object.assign(pctLabel.style, { minWidth: "3rem", textAlign: "right", fontSize: "0.85rem" });

        rotRow.append(rotLabel, slider, pctLabel);
        wrapper.appendChild(rotRow);

        // Size control
        const sizeRow = document.createElement("div");
        sizeRow.className = "form-group mb-2";
        Object.assign(sizeRow.style, { display: "flex", alignItems: "center", gap: "0.75rem" });

        const sizeLabel = document.createElement("label");
        sizeLabel.textContent = "Plot size";
        Object.assign(sizeLabel.style, { margin: "0", whiteSpace: "nowrap", fontSize: "0.85rem" });

        const sizeSlider = document.createElement("input");
        sizeSlider.type = "range";
        sizeSlider.className = "custom-range";
        sizeSlider.min = "300";
        sizeSlider.max = "1500";
        sizeSlider.value = String(PLOT_DEFAULT_SIZE);
        sizeSlider.step = "50";
        sizeSlider.style.flex = "1";

        const sizeValueLabel = document.createElement("span");
        sizeValueLabel.textContent = PLOT_DEFAULT_SIZE + "px";
        Object.assign(sizeValueLabel.style, { minWidth: "4rem", textAlign: "right", fontSize: "0.85rem" });

        sizeRow.append(sizeLabel, sizeSlider, sizeValueLabel);
        wrapper.appendChild(sizeRow);

        // SVG container (replaced on rotation or size change)
        const svgContainer = document.createElement("div");
        svgContainer.appendChild(renderDotPlot(delta, 0, PLOT_DEFAULT_SIZE));
        wrapper.appendChild(svgContainer);

        function rerender() {
            const rot = parseInt(slider.value, 10) / 100;
            const size = parseInt(sizeSlider.value, 10);
            pctLabel.textContent = parseInt(slider.value, 10) + "%";
            sizeValueLabel.textContent = size + "px";
            svgContainer.innerHTML = "";
            svgContainer.appendChild(renderDotPlot(delta, rot, size));
        }

        slider.addEventListener("input", rerender);
        sizeSlider.addEventListener("input", rerender);

        const btnRow = document.createElement("div");
        btnRow.style.marginTop = "1rem";
        Object.assign(btnRow.style, { display: "flex", gap: "0.5rem" });
        btnRow.appendChild(buildLoadAnotherButton());

        const dlBtn = document.createElement("button");
        dlBtn.type = "button";
        dlBtn.className = "btn btn-outline-secondary";
        dlBtn.textContent = "Download SVG";
        dlBtn.addEventListener("click", () => {
            const svgEl = svgContainer.querySelector("svg");
            if (!svgEl) return;
            const serializer = new XMLSerializer();
            const svgStr = '<?xml version="1.0" encoding="UTF-8"?>\n' +
                serializer.serializeToString(svgEl);
            const blob = new Blob([svgStr], { type: "image/svg+xml;charset=utf-8" });
            const url = URL.createObjectURL(blob);
            const a = document.createElement("a");
            a.href = url;
            a.download = "dotplot.svg";
            a.click();
            URL.revokeObjectURL(url);
        });
        btnRow.appendChild(dlBtn);

        wrapper.appendChild(btnRow);

        return wrapper;
    }

    function buildLoadAnotherButton() {
        const btn = document.createElement("button");
        btn.type = "button";
        btn.className = "btn btn-outline-secondary";
        btn.textContent = "Load another file";
        btn.addEventListener("click", showDropZone);
        return btn;
    }

    // Public API for loading content programmatically (used by "load example")
    return { handleFileContent };
}
