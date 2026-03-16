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

const PLOT = {
    dataWidth: 700,
    dataHeight: 700,
    marginLeft: 100,
    marginBottom: 60,
    marginTop: 20,
    marginRight: 20,
};
PLOT.totalWidth = PLOT.marginLeft + PLOT.dataWidth + PLOT.marginRight;
PLOT.totalHeight = PLOT.marginTop + PLOT.dataHeight + PLOT.marginBottom;

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

function computeLayout(delta) {
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
        scaleX: totalRefLen > 0 ? PLOT.dataWidth / totalRefLen : 1,
        scaleY: totalQueryLen > 0 ? PLOT.dataHeight / totalQueryLen : 1,
    };
}

// Coordinate transforms
function toSvgX(layout, seqName, pos) {
    const offset = layout.refOffsets.get(seqName) || 0;
    return PLOT.marginLeft + (offset + pos) * layout.scaleX;
}

function toSvgY(layout, seqName, pos) {
    const offset = layout.queryOffsets.get(seqName) || 0;
    return PLOT.marginTop + PLOT.dataHeight - (offset + pos) * layout.scaleY;
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

function renderDotPlot(delta) {
    const layout = computeLayout(delta);
    const viewBox = `0 0 ${PLOT.totalWidth} ${PLOT.totalHeight}`;

    const svg = svgEl("svg", {
        viewBox,
        width: "100%",
        style: "max-width: 900px",
    });

    // Background
    svg.appendChild(
        svgEl("rect", {
            x: 0, y: 0,
            width: PLOT.totalWidth, height: PLOT.totalHeight,
            fill: "white",
        })
    );

    // Data area border
    svg.appendChild(
        svgEl("rect", {
            x: PLOT.marginLeft, y: PLOT.marginTop,
            width: PLOT.dataWidth, height: PLOT.dataHeight,
            fill: "white", stroke: "#ccc", "stroke-width": 1,
        })
    );

    // Alignment lines
    const alignG = svgEl("g", {});
    for (const section of delta.alignmentSections) {
        for (const a of section.alignments) {
            const isForward = a.queryStart <= a.queryEnd;
            alignG.appendChild(
                svgEl("line", {
                    x1: toSvgX(layout, section.refId, a.refStart),
                    y1: toSvgY(layout, section.queryId, a.queryStart),
                    x2: toSvgX(layout, section.refId, a.refEnd),
                    y2: toSvgY(layout, section.queryId, a.queryEnd),
                    stroke: isForward ? "#00BFFF" : "#9933FF",
                    "stroke-width": 2,
                    "stroke-linecap": "round",
                })
            );
        }
    }
    svg.appendChild(alignG);

    // X-axis ticks
    const tickInterval = niceInterval(layout.totalRefLen);
    const ticks = generateTicks(0, layout.totalRefLen, tickInterval);
    const axisY = PLOT.marginTop + PLOT.dataHeight;
    const xAxisG = svgEl("g", {});
    for (const tick of ticks) {
        const x = PLOT.marginLeft + tick * layout.scaleX;
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
        const y = PLOT.marginTop + PLOT.dataHeight - midpoint * layout.scaleY;
        yAxisG.appendChild(
            svgEl("text", {
                x: PLOT.marginLeft - 8, y,
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
            const y = PLOT.marginTop + PLOT.dataHeight - offset * layout.scaleY;
            sepG.appendChild(
                svgEl("line", {
                    x1: PLOT.marginLeft, y1: y,
                    x2: PLOT.marginLeft + PLOT.dataWidth, y2: y,
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
    const xLabelY = PLOT.marginTop + PLOT.dataHeight + 45;
    svg.appendChild(
        svgEl("text", {
            x: PLOT.marginLeft + PLOT.dataWidth / 2, y: xLabelY,
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
        wrapper.appendChild(renderDotPlot(delta));

        const btnRow = document.createElement("div");
        btnRow.style.marginTop = "1rem";
        btnRow.appendChild(buildLoadAnotherButton());
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
