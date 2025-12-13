// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

// --- D3 Setup ---
const svg = d3.select("#graph");
const g = svg.append("g");
const tooltip = d3.select("#tooltip");

svg.call(d3.zoom().scaleExtent([0.2, 8])
    .on("zoom", (e) => g.attr("transform", e.transform)));

let originalData = {}; 
let currentGraphKey = "json/graph1.json"; 

const shapeMap = { circle: d3.symbolCircle, box: d3.symbolCircle, triangle: d3.symbolTriangle, diamond: d3.symbolDiamond, hexagon: d3.symbolCross, octagon: d3.symbolStar, parallelogram: d3.symbolWye, trapezium: d3.symbolSquare,  };

// --- DATA LOADING & INITIALIZATION ---
function loadAndRenderGraph(fileKey) {
    currentGraphKey = fileKey;
    resetFilters(false);

    const isColoView = fileKey.endsWith("graph2.json");

    if (isColoView) {
        disableFiltersForColoView();
    } else {
        enableAllFilters();
    }
    if (originalData[fileKey]) {
        populateFilters(originalData[fileKey]);
        applyFiltersAndDraw();
    } else {
        d3.json(fileKey).then(data => {
            originalData[fileKey] = data;
            populateFilters(data);
            applyFiltersAndDraw();
        }).catch(error => console.error("Error loading JSON:", error));
    }
}

// --- POPULATE FILTER DROPDOWNS ---
function populateFilters(data) {
    const menu = document.querySelector("#mgeGroupMenu");
    menu.innerHTML = `
        <li><a class="dropdown-item active" data-value="all">All Groups</a></li>
    `;
    const nodeSource = currentGraphKey.includes("graph1")
        ? data.nodes.filter(n => !n.isARG)
        : data.nodes;
    const groups = [...new Set(nodeSource.map(n => n.mgeGroup).filter(Boolean))].sort();
    groups.forEach(g => {
        menu.innerHTML += `<li><a class="dropdown-item" data-value="${g}">${g}</a></li>`;
    });
    bindCustomDropdownHandlers();

}

function bindCustomDropdownHandlers() {
  document
    .querySelectorAll(".dropdown-select .dropdown-item")
    .forEach(item => {

      item.onclick = function () {
        const value = this.dataset.value;

        const dropdown = this.closest(".dropdown");
        const hiddenSelector = dropdown.dataset.targetInput;
        const hiddenInput = hiddenSelector ? document.querySelector(hiddenSelector) : null;
        const button = dropdown.querySelector("button.dropdown-btn");

        hiddenInput.value = value;
        button.textContent = this.textContent.trim();

        this.closest(".dropdown-menu")
            .querySelectorAll(".dropdown-item")
            .forEach(i => i.classList.remove("active"));
        this.classList.add("active");

        // trigger D3 logic
        hiddenInput.dispatchEvent(new Event("change", { bubbles: true }));
      };
    });
}



// --- RESET FILTERS ---
function resetFilters(redraw = true) {
    // hidden values used by D3
    d3.select("#diseaseFilter").property("value", "all");
    d3.select("#mgeGroupFilter").property("value", "all");
    d3.select("#argSearch").property("value", "");
    d3.select("#mgeSearch").property("value", "");
    d3.selectAll(".timepoint-checkbox").property("checked", true);

    resetSingleDropdown("#diseaseFilter", "All Diseases");
    resetSingleDropdown("#mgeGroupFilter", "All Groups");
    updateTimepointButtonText();

    if (redraw) {
        applyFiltersAndDraw();
    }
}

function resetSingleDropdown(hiddenSelector, labelText) {
    const dropdown = document.querySelector(
        `.dropdown[data-target-input='${hiddenSelector}']`
    );
    if (!dropdown) return;

    const button = dropdown.querySelector(".dropdown-btn");
    if (button) button.textContent = labelText;

    const menuItems = dropdown.querySelectorAll(".dropdown-item");
    menuItems.forEach(item => {
        const isAll = item.dataset.value === "all";
        item.classList.toggle("active", isAll);
    });
}


// --- CORE FILTERING LOGIC ---
function applyFiltersAndDraw() {
    if (!originalData[currentGraphKey]) return;

    let data = JSON.parse(JSON.stringify(originalData[currentGraphKey])); 

    const filters = {
        disease: d3.select("#diseaseFilter").property("value"),
        mgeGroup: d3.select("#mgeGroupFilter").property("value"),
        timepoints: Array.from(
            d3.selectAll(".timepoint-checkbox").nodes()
        )
        .filter(cb => cb.checked)
        .map(cb => cb.value),
        argSearchTerm: d3.select("#argSearch").property("value").trim().toLowerCase(),
        mgeSearchTerm: d3.select("#mgeSearch").property("value").trim().toLowerCase(),
    };

    let { nodes, links } = data;
    
    // --- FILTERING PIPELINE ---

    let strictlyFilteredNodeIds = getStrictlyFilteredNodeIds(nodes, links, filters);

    let seedNodeIds = getSeedNodeIds(nodes, filters, strictlyFilteredNodeIds);

    let finalVisibleNodeIds;
    // Add patientCount to colocalization links (graph1 only)
    if (currentGraphKey.includes("graph1")) {
        const activeDisease = d3.select("#diseaseFilter").property("value");

        data.links.forEach(link => {
            if (!link.isColo) return;

            const dc = link.diseaseCounts || {};

            if (activeDisease !== "all") {
                link.patientCount = Number(dc[activeDisease]) || 0;
            } else {
                link.patientCount = Object.values(dc)
                    .reduce((sum, v) => sum + (Number(v) || 0), 0);
            }
        });
    }

    if (filters.mgeGroup !== 'all' || filters.argSearchTerm || filters.mgeSearchTerm) {
        if (seedNodeIds.size === 0) {
            finalVisibleNodeIds = new Set();
        } else {
            const neighborIds = getNeighborIds(links, seedNodeIds);
            const allowedNeighbors = [...neighborIds].filter(id => strictlyFilteredNodeIds.has(id));
            const allowedSeeds = [...seedNodeIds].filter(id => strictlyFilteredNodeIds.has(id));
            finalVisibleNodeIds = new Set([...allowedSeeds, ...allowedNeighbors]);
        }
    } else {
        finalVisibleNodeIds = strictlyFilteredNodeIds;
    }

    const finalNodes = nodes.filter(n => finalVisibleNodeIds.has(n.id));
    let finalLinks = links.filter(l => {
        const sourceId = typeof l.source === 'object' ? l.source.id : l.source;
        const targetId = typeof l.target === 'object' ? l.target.id : l.target;
        return finalVisibleNodeIds.has(sourceId) && finalVisibleNodeIds.has(targetId);
    });
        
    updateVisualization({ nodes: finalNodes, links: finalLinks });
}
// -- Disable filters for Colocalization View ---
function disableFiltersForColoView() {
    d3.select("#argSearch").property("value", "").attr("disabled", true);
    d3.select("#mgeSearch").property("value", "").attr("disabled", true);
    d3.select("#toggleColo").property("checked", true).attr("disabled", true);
}

// -- Enable filters ---
function enableAllFilters() {
    d3.select("#argSearch").attr("disabled", null);
    d3.select("#mgeSearch").attr("disabled", null);
    d3.select("#toggleColo").attr("disabled", null);
}


// --- FILTERING HELPERS ---
function getStrictlyFilteredNodeIds(nodes, links, filters) {
    let visibleNodeIds = new Set(nodes.map(n => n.id));

    // Timepoint Filter
    if (filters.timepoints && filters.timepoints.length > 0) {
        const allowedIds = new Set(
            nodes
                .filter(node => filters.timepoints.includes(node.timepointCategory))
                .map(n => n.id)
        );
        visibleNodeIds = new Set([...visibleNodeIds].filter(id => allowedIds.has(id)));
    } else {
        visibleNodeIds.clear();
    }
    
    // Disease Filter
    if (filters.disease !== 'all') {
        let diseaseFilteredIds = new Set();
        if (currentGraphKey.includes('graph1')) {
            links.forEach(link => {
                if (link.isColo && link.diseases && link.diseases.includes(filters.disease)) {
                    diseaseFilteredIds.add(typeof link.source === 'object' ? link.source.id : link.source);
                    diseaseFilteredIds.add(typeof link.target === 'object' ? link.target.id : link.target);
                }
            });
        } else {
            nodes.forEach(node => {
                if (node.diseases && node.diseases.includes(filters.disease)) {
                    diseaseFilteredIds.add(node.id);
                }
            });
        }
        visibleNodeIds = new Set([...visibleNodeIds].filter(id => diseaseFilteredIds.has(id)));
    }
    
    return visibleNodeIds;
}

function getSeedNodeIds(nodes, filters, availableNodeIds) {
    let mgeGroupSeedIds = new Set();
    let searchSeedIds = new Set();

    const availableNodes = nodes.filter(n => availableNodeIds.has(n.id));

    if (filters.mgeGroup !== 'all') {
        availableNodes.forEach(n => {
            if (n.mgeGroup === filters.mgeGroup) {
                mgeGroupSeedIds.add(n.id);
            }
        });
    }

    const hasArgSearch = !!filters.argSearchTerm;
    const hasMgeSearch = !!filters.mgeSearchTerm;

    if (hasArgSearch || hasMgeSearch) {
        availableNodes.forEach(n => {
            const label = n.label ? n.label.toLowerCase() : "";
            if (hasArgSearch && n.isARG && label.includes(filters.argSearchTerm)) {
                searchSeedIds.add(n.id);
            }
            if (hasMgeSearch && !n.isARG && label.includes(filters.mgeSearchTerm)) {
                searchSeedIds.add(n.id);
            }
        });
    }

    const isMgeGroupFiltered = filters.mgeGroup !== 'all';
    const isSearchFiltered = hasArgSearch || hasMgeSearch;

    if (isMgeGroupFiltered && isSearchFiltered) {
        return new Set([...mgeGroupSeedIds].filter(id => searchSeedIds.has(id)));
    }
    if (isMgeGroupFiltered) {
        return mgeGroupSeedIds;
    }
    if (isSearchFiltered) {
        return searchSeedIds;
    }
    
    return new Set();
}

function getNeighborIds(links, seedNodeIds) {
    const neighborIds = new Set();
    links.forEach(link => {
        const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
        const targetId = typeof link.target === 'object' ? link.target.id : link.target;
        if (seedNodeIds.has(sourceId)) neighborIds.add(targetId);
        if (seedNodeIds.has(targetId)) neighborIds.add(sourceId);
    });
    return neighborIds;
}


// --- D3 RENDERING ---
function updateVisualization(data) {
    g.selectAll("*").remove();
    if (!data.nodes.length) return;

    const defs = svg.append("defs");
    const colors = [...new Set(data.links.map(d => d.color || "#999"))];
    colors.forEach(c => {
        defs.append("marker")
            .attr("id", `arrow-${c.replace("#", "")}`)
            .attr("viewBox", "0 -5 10 10")
            .attr("refX", 12)
            .attr("refY", 0)
            .attr("markerWidth", 3)
            .attr("markerHeight", 3)
            .attr("orient", "auto")
            .attr("markerUnits", "strokeWidth")
            .append("path")
            .attr("d", "M0,-5L10,0L0,5")
            .attr("fill", c);
    });

    const simNodes = data.nodes.map(d => ({ ...d }));
    const simLinks = data.links.map(d => ({ ...d }));

    // --- compute patient count only if graph2.json - colocalization view ---
    if (currentGraphKey.includes("graph2")) {
        const activeDisease = d3.select("#diseaseFilter").property("value");
        simNodes.forEach(n => {
            const dc = (n && typeof n.diseaseCounts === "object" && n.diseaseCounts) ? n.diseaseCounts : {};
            console.log(`Node: ${n.id}, raw diseaseCounts:`, dc);
            if (activeDisease && activeDisease !== "all") {
                n.patientCount = Number(dc[activeDisease]) || 0;
            } else {
                n.patientCount = Object.values(dc).reduce((sum, v) => sum + (Number(v) || 0), 0);
            }
            console.log(`Node: ${n.id}, PatientCount: ${n.patientCount}`);
        });
        const uniqCounts = new Set(simNodes.map(n => n.patientCount));
        if (uniqCounts.size === 1) {
            console.warn("[graph2] All nodes have identical patientCount =", [...uniqCounts][0], "— scaling disabled.");
            simNodes.forEach(n => { n.__disableScaling = true; });
        }
    }

    // --- links ---
    const linkSelection = g.selectAll("path.link")
        .data(simLinks, d => `${d.source.id}-${d.target.id}-${d.type}`)
        .join("path")
        .attr("class", "link")
        .attr("stroke", d => d.color || "#999")
        .attr("marker-end", d => d.isColo ? null : `url(#arrow-${(d.color || "#999").replace("#", "")})`)
        .attr("stroke-width", d => Math.max(1, d.penwidth || 1))
        .attr("stroke-dasharray", d => d.isColo ? null : "4 2")

        .on("mouseover", function(event, d) {
            if (!currentGraphKey.includes("graph1")) return;  // Only graph1.json

            const count = d.individualCount ?? d.patientCount ?? 0;

            tooltip
                .style("opacity", 1)
                .html(`<strong>Patients:</strong> ${count}`);

            d3.select(this).attr("stroke-width", (d.penwidth || 2) + 2);
        })
        .on("mousemove", function(event) {
            tooltip
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 20) + "px");
        })
        .on("mouseout", function(event, d) {
            if (!currentGraphKey.includes("graph1")) return;

            tooltip.style("opacity", 0);
            d3.select(this).attr("stroke-width", d.penwidth || 2);
        });


    // --- nodes ---
    const nodeSelection = g.selectAll("path.node")
        .data(simNodes, d => d.id)
        .join("path")
        .attr("class", "node")
        .attr("d", d3.symbol()
            .type(d => shapeMap[d.shape] || d3.symbolCircle)
            .size(d => {
                if (!currentGraphKey.includes("graph2") || d.__disableScaling) {
                    return 300;
                }
                const count = (typeof d.patientCount === "number") ? d.patientCount : 0;

                // Smooth scaling between 150 and 4000 for counts 1–30
                const MIN_COUNT = 1;
                const MAX_COUNT = 30;
                const MIN_SIZE = 150;
                const MAX_SIZE = 4000;

                // Use square-root scaling for gradual increase
                const scale = d3.scaleSqrt()
                    .domain([MIN_COUNT, MAX_COUNT])
                    .range([MIN_SIZE, MAX_SIZE]);

                return scale(Math.max(MIN_COUNT, Math.min(MAX_COUNT, count)));
            })
        )
        .attr("fill", d => d.color)
        .attr("stroke", "grey")
        .attr("stroke-width", 1.5)
        .call(d3.drag().on("start", dragstart).on("drag", dragged).on("end", dragend));

    nodeSelection.append("title").text(d => {
        if (currentGraphKey.includes("graph2")) {
            return `${d.label}\nPatients: ${d.patientCount ?? 0}`;
        }
        return d.label;
    });

    // --- labels ---
    const labelSelection = g.selectAll("text.label")
        .data(simNodes, d => d.id)
        .join("text")
        .attr("class", "label")
        .attr("dy", -12)
        .text(d => d.label);

    g.selectAll("text.label").style("display", d3.select("#toggleLabels").property("checked") ? "block" : "none");

    // --- simulation ---
    const sim = d3.forceSimulation(simNodes)
        .force("link", d3.forceLink(simLinks).id(d => d.id).distance(d => d.isColo ? 40 : 60))
        .force("charge", d3.forceManyBody().strength(d => -(40 + (d.degree || 0) * 15)))
        .force("center", d3.forceCenter(500, 350))
        .force("collision", d3.forceCollide().radius(20))
        .force("x", d3.forceX(500).strength(0.05))
        .force("y", d3.forceY(350).strength(0.05))
        .on("tick", ticked);

    function ticked() {
        linkSelection.attr("d", d => linkArc(d));
        nodeSelection.attr("transform", d => `translate(${d.x},${d.y})`);
        labelSelection.attr("x", d => d.x).attr("y", d => d.y);
    }

    function dragstart(event, d) { if (!event.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; }
    function dragged(event, d) { d.fx = event.x; d.fy = event.y; }
    function dragend(event, d) { if (!event.active) sim.alphaTarget(0); d.fx = null; d.fy = null; }

    updateLinkVisibility();
}



function linkArc(d) {
    const r = Math.hypot(d.target.x - d.source.x, d.target.y - d.source.y);
    return `M${d.source.x},${d.source.y}A${r},${r} 0 0,1 ${d.target.x},${d.target.y}`;
}

// ---- EDGE VISIBILITY TOGGLING ---
function updateLinkVisibility() {
    const showColo = d3.select("#toggleColo").property("checked");
    const showTemporal = d3.select("#toggleTemporal").property("checked");

    g.selectAll("path.link")
        .style("display", d => {
            if (d.isColo) {
                return showColo ? "inline" : "none";
            } else {
                return showTemporal ? "inline" : "none";
            }
        });
}

// --- Generic handler for custom single-select dropdowns (e.g., Disease) ---
document.querySelectorAll(".dropdown-select .dropdown-item").forEach(item => {
  item.addEventListener("click", function () {
    const value = this.dataset.value;

    // Find the wrapper <div class="dropdown" ...>
    const dropdown = this.closest(".dropdown");
    if (!dropdown) return;
    const hiddenSelector = dropdown.dataset.targetInput;
    const hiddenInput = hiddenSelector ? document.querySelector(hiddenSelector) : null;
    if (!hiddenInput) return;

    const button = dropdown.querySelector("button.dropdown-btn");
    if (!button) return;
    hiddenInput.value = value;
    button.textContent = this.textContent.trim();

    // Mark active item
    this.closest(".dropdown-menu")
      .querySelectorAll(".dropdown-item")
      .forEach(i => i.classList.remove("active"));
    this.classList.add("active");

    // Fire a real "change" event on the hidden input so D3 listeners run
    hiddenInput.dispatchEvent(new Event("change", { bubbles: true }));
  });
});



// --- PATIENT STAGES FILTER DROPDOWN HANDLING ---
function bindTimepointListeners() {
  const dropdownMenu = document.querySelector(".dropdown-menu[aria-labelledby='timepointFilter']");
  if (!dropdownMenu) return;

  // Keep dropdown open when clicking inside
  dropdownMenu.addEventListener("click", (event) => {
    if (event.target.classList.contains("timepoint-checkbox")) {
      event.stopPropagation();
    }
  });

  // Apply filtering when checkbox state changes
  dropdownMenu.addEventListener("change", (event) => {
    if (event.target.classList.contains("timepoint-checkbox")) {
      console.log("Checkbox changed:", event.target.value, event.target.checked);
      applyFiltersAndDraw();
      updateTimepointButtonText();
    }
  });
}

// Update the button text based on selected timepoints in patient stages filter
function updateTimepointButtonText() {
  const checked = Array.from(document.querySelectorAll(".timepoint-checkbox:checked"));
  const button = document.getElementById("timepointFilter");
  if (!button) return;

  if (checked.length === 5) { button.textContent = "All Stages"; return; }
  if (checked.length === 0) { button.textContent = "No Stages"; return; }

  if (checked.length === 5) {
    button.textContent = "All Stages";
    button.title = "All Stages";
    return;
  }
  if (checked.length === 0) {
    button.textContent = "No Stages";
    button.title = "No Stages";
    return;
  }

  const pretty = {
    donor: "Donor",
    pre: "Pre-FMT",
    post1: "Post-FMT (1–30 d)",
    post2: "Post-FMT (31–60 d)",
    post3: "Post-FMT (61+ d)"
  };
  button.textContent = checked.map(cb => pretty[cb.value]).join(", ");
}

// --- SVG DOWNLOAD ---
function downloadCurrentGraphAsSVG() {
  const svgEl = document.querySelector("#graph");
  if (!svgEl) return;
  const clone = svgEl.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

  // Inline any external CSS 
  const cssStyles = Array.from(document.styleSheets)
    .map(ss => {
      try { return Array.from(ss.cssRules).map(r => r.cssText).join("\n"); }
      catch (e) { return ""; } 
    })
    .join("\n");

  const styleEl = document.createElement("style");
  styleEl.textContent = cssStyles;
  clone.insertBefore(styleEl, clone.firstChild);

  // Serialize SVG to string
  const serializer = new XMLSerializer();
  const svgString = serializer.serializeToString(clone);

  // Create downloadable Blob
  const blob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);

  // Create a temporary <a> to trigger the download
  const a = document.createElement("a");
  const fileName = currentGraphKey.replace(/^.*[\\/]/, '').replace('.json', '') + ".svg";
  a.href = url;
  a.download = fileName;
  document.body.appendChild(a);
  a.click();

  // Clean up
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}


// --- EVENT LISTENERS ---
d3.select("#dataset").on("change", function() { loadAndRenderGraph(this.value); });
d3.select("#diseaseFilter").on("change", applyFiltersAndDraw);
d3.select("#mgeGroupFilter").on("change", applyFiltersAndDraw);
d3.select("#timepointFilter").on("change", applyFiltersAndDraw);
d3.select("#searchBtn").on("click", applyFiltersAndDraw);
d3.select("#resetBtn").on("click", resetFilters);
d3.select("#argSearch").on("keydown", event => { if (event.key === 'Enter') { applyFiltersAndDraw(); } });
d3.select("#mgeSearch").on("keydown", event => { if (event.key === 'Enter') { applyFiltersAndDraw(); } });
d3.select("#toggleLabels").on("change", () => g.selectAll("text.label").style("display", d3.select("#toggleLabels").property("checked") ? "block" : "none"));
d3.select("#toggleColo").on("change", updateLinkVisibility);
d3.select("#toggleTemporal").on("change", updateLinkVisibility);
d3.selectAll(".timepoint-checkbox").on("change", applyFiltersAndDraw);
document.getElementById("downloadSvgBtn").addEventListener("click", downloadCurrentGraphAsSVG);


const legendOverlay = document.getElementById("legendOverlay");
document.getElementById("toggleLegend").addEventListener("change", function() {
	if (this.checked) {
		legendOverlay.classList.add("visible");
	} else {
		legendOverlay.classList.remove("visible");
	}
});

// --- INITIAL LOAD ---
bindTimepointListeners();
loadAndRenderGraph(currentGraphKey);
updateTimepointButtonText();

