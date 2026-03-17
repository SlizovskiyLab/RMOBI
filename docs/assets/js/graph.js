// Released under the GNU GPLv3; see LICENSE for details.
// Developed by Boucher Lab and Slizovskiy Lab.

/* 
ENRICHMENT FUNCTIONS
*/

function getTimepointCategory(timepoint) {
  const t = Number(timepoint);
  if (t === 1000) return "donor";
  if (t === 0) return "pre";
  if (t > 0 && t < 31) return "post1";
  if (t > 30 && t < 61) return "post2";
  if (t > 60) return "post3";
  return "unknown";
}

function getTimepointColor(timepoint) {
  const t = Number(timepoint);
  if (t === 1000) return "yellow";
  if (t === 0) return "red";
  if (t > 0 && t < 31) return "#99D2FF";
  if (t > 30 && t < 61) return "#4D9DFF";
  if (t > 60) return "#3A6EFF";
  return "green"; // fallback
}

function getMGEGroupShape(groupName) {
  if (
    groupName === "plasmid" ||
    groupName === "Colicin_plasmid" ||
    groupName === "Inc_plasmid"
  ) {
    return "diamond";
  }
  if (groupName === "prophage") {
    return "hexagon";
  }
  if (groupName === "virus") {
    return "triangle";
  }
  if (groupName === "ICE" || groupName === "ICEberg") {
    return "octagon";
  }
  if (groupName === "replicon") {
    return "parallelogram";
  }
  if (groupName === "likely IS/TE") {
    return "trapezium";
  }
  return "box"; // default
}

function isPostCategory(tpCategory) {
  return tpCategory !== "donor" && tpCategory !== "pre";
}

function getTemporalLinkColor(srcCategory, tgtCategory) {
  const tgtIsPost = isPostCategory(tgtCategory);

  if (srcCategory === "donor" && tgtCategory === "pre") return "#006400";  // donor -> pre
  if (srcCategory === "donor" && tgtIsPost) return "#4B0082";             // donor -> post
  if (srcCategory === "pre" && tgtIsPost) return "orange";                // pre -> post
  return "black";                                                          // post -> post (and any other)
}

function getPenWidth(weight) {
  const w = Number(weight ?? 1);
  let penwidth = 4.0;
  if (w > 1) penwidth = 4.0 + (w - 1) * 2.0;
  return Math.min(10.0, penwidth);
}

function getDiseasesFromCounts(diseaseCounts) {
  if (!diseaseCounts) return [];
  return Object.keys(diseaseCounts);
}


function enrichNodes(nodes) {
  for (const n of nodes) {
    n.timepointCategory = getTimepointCategory(n.timepoint);
    n.color = getTimepointColor(n.timepoint);
    if (!n.isARG) {
        n.shape = getMGEGroupShape(n.mgeGroup);
    } else {
        n.shape = "circle";
    }
    n.diseases = getDiseasesFromCounts(n.diseaseCounts);
  }
  return nodes;
}

function enrichLinks(links, nodeById /* optional */) {
  for (const e of links) {
    const isColo = !!e.isColo;

    // Colocalization styling
    if (isColo) {
      e.type = "colocalization";
      e.style = "solid";
      e.color = "#696969";
      e.penwidth = getPenWidth(e.weight);
      continue;
    }

    // Temporal styling
    e.type = "temporal";
    e.style = "dashed";
    e.penwidth = getPenWidth(e.weight);

    // pre-attached timepoints (fast)
    let srcTp = e.sourceTimepoint;
    let tgtTp = e.targetTimepoint;

    // Fallback
    if ((srcTp == null || tgtTp == null) && nodeById) {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      srcTp = nodeById.get(sid)?.timepoint ?? null;
      tgtTp = nodeById.get(tid)?.timepoint ?? null;
      e.sourceTimepoint = srcTp;
      e.targetTimepoint = tgtTp;
    }

    const srcCat = getTimepointCategory(srcTp);
    const tgtCat = getTimepointCategory(tgtTp);

    e.color = getTemporalLinkColor(srcCat, tgtCat);
  }

  return links;
}

/********************************************************************/ 

// --- D3 Setup ---
const svg = d3.select("#graph");
const g = svg.select("g.main")
  .empty()
  ? svg.append("g").attr("class", "main")
  : svg.select("g.main");
const tooltip = d3.select("#tooltip");
let currentRenderedData = { nodes: [], links: [] };

resizeGraphViewBox();
window.addEventListener("resize", resizeGraphViewBox);

const zoom = d3.zoom()
  .scaleExtent([0.1, 5])
  .on("zoom", (e) => {
    g.attr("transform", e.transform);
  });

svg.call(zoom);
setInitialZoom(0.6); 

function setInitialZoom(k = 0.25) {
  const vb = svg.node().viewBox.baseVal; // 0 0 1000 1000
  const cx = vb.x + vb.width / 2;
  const cy = vb.y + vb.height / 2;

  // zoom around the viewBox center
  svg.call(
    zoom.transform,
    d3.zoomIdentity.translate(cx, cy).scale(k).translate(-cx, -cy)
  );
}

function resizeGraphViewBox() {
  const el = document.getElementById("graph");
  if (!el) return;

  const r = el.getBoundingClientRect();
  const w = Math.max(1, Math.round(r.width));
  const h = Math.max(1, Math.round(r.height));

  d3.select(el).attr("viewBox", `0 0 ${w} ${h}`);
}



// Optional: disable dblclick zoom
// svg.on("dblclick.zoom", null);

const zoomStep = 1.2;

// Zoom around the SVG center, for buttons
function zoomAtCenter(k) {
  const node = svg.node();
  if (!node) return;

  const { width, height } = node.getBoundingClientRect();
  const center = [width / 2, height / 2];

  svg.transition().duration(180).call(zoom.scaleBy, k, center);
}

document.getElementById("zoom-in")?.addEventListener("click", () => {
  zoomAtCenter(zoomStep);
});

document.getElementById("zoom-out")?.addEventListener("click", () => {
  zoomAtCenter(1 / zoomStep);
});

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
    updateGraphStatsVisibility();
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

    enrichNodes(data.nodes);
    const nodeById = new Map(data.nodes.map(n => [n.id, n]));
    enrichLinks(data.links, nodeById);

    const filters = {
        disease: d3.select("#diseaseFilter").property("value"),
        mgeGroup: d3.select("#mgeGroupFilter").property("value"),
        timepoints: Array.from(d3.selectAll(".timepoint-checkbox").nodes())
            .filter(cb => cb.checked)
            .map(cb => cb.value),
        argSearchTerm: d3.select("#argSearch").property("value").trim().toLowerCase(),
        mgeSearchTerm: d3.select("#mgeSearch").property("value").trim().toLowerCase(),
    };

    let { nodes, links } = data;

    // -------------------------------------------------
    // STEP 1: Build working link set for graph1 disease
    // -------------------------------------------------
    let workingLinks = links.map(l => ({ ...l }));

    if (currentGraphKey.includes("graph1")) {
        workingLinks = workingLinks
            .filter(l => {
                if (!l.isColo) return true; // keep temporal links
                if (filters.disease === "all") return true;

                const dc = l.diseaseCounts || {};
                return (Number(dc[filters.disease]) || 0) > 0;
            })
            .map(l => {
                if (!l.isColo) return l;

                const dc = l.diseaseCounts || {};
                l.patientCount =
                    filters.disease === "all"
                        ? Object.values(dc).reduce((sum, v) => sum + (Number(v) || 0), 0)
                        : Number(dc[filters.disease]) || 0;

                return l;
            });
    }

    // -------------------------------------------------
    // STEP 2: For graph2, assign node patientCount 
    // -------------------------------------------------
    let workingNodes = nodes.map(n => ({ ...n }));

    if (!currentGraphKey.includes("graph1")) {
        workingNodes = workingNodes.map(n => {
            const dc = n.diseaseCounts || {};
            n.patientCount =
                filters.disease === "all"
                    ? Object.values(dc).reduce((sum, v) => sum + (Number(v) || 0), 0)
                    : Number(dc[filters.disease]) || 0;
            return n;
        });
    }

    // -------------------------------------------------
    // STEP 3: Strict filtering should use workingLinks
    // -------------------------------------------------
    let strictlyFilteredNodeIds = getStrictlyFilteredNodeIds(
        workingNodes,
        workingLinks,
        filters
    );

    let seedNodeIds = getSeedNodeIds(workingNodes, filters, strictlyFilteredNodeIds);

    let finalVisibleNodeIds;

    if (filters.mgeGroup !== "all" || filters.argSearchTerm || filters.mgeSearchTerm) {
        if (seedNodeIds.size === 0) {
            finalVisibleNodeIds = new Set();
        } else {
            // IMPORTANT: use workingLinks, not original links
            const neighborIds = getNeighborIds(workingLinks, seedNodeIds);
            const allowedNeighbors = [...neighborIds].filter(id => strictlyFilteredNodeIds.has(id));
            const allowedSeeds = [...seedNodeIds].filter(id => strictlyFilteredNodeIds.has(id));
            finalVisibleNodeIds = new Set([...allowedSeeds, ...allowedNeighbors]);
        }
    } else {
        finalVisibleNodeIds = strictlyFilteredNodeIds;
    }

    // -------------------------------------------------
    // STEP 4: Build final nodes/links from working sets
    // -------------------------------------------------
    let finalNodes = workingNodes.filter(n => finalVisibleNodeIds.has(n.id));

    let finalLinks = workingLinks.filter(l => {
        const sourceId = typeof l.source === "object" ? l.source.id : l.source;
        const targetId = typeof l.target === "object" ? l.target.id : l.target;
        return finalVisibleNodeIds.has(sourceId) && finalVisibleNodeIds.has(targetId);
    });

    // remove isolated nodes after link pruning
    if (currentGraphKey.includes("graph1")) {
        const linkedNodeIds = new Set();
        finalLinks.forEach(l => {
            const sourceId = typeof l.source === "object" ? l.source.id : l.source;
            const targetId = typeof l.target === "object" ? l.target.id : l.target;
            linkedNodeIds.add(sourceId);
            linkedNodeIds.add(targetId);
        });
        finalNodes = finalNodes.filter(n => linkedNodeIds.has(n.id));
    }

    populateSearchSuggestionsFromNodes(finalNodes);

    currentRenderedData = { nodes: finalNodes, links: finalLinks };

    updateVisualization(currentRenderedData);
    const stats = computeGraphStats(finalNodes, finalLinks);
    renderGraphStats(stats);
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
    if (filters.disease !== "all") {
        let diseaseFilteredIds = new Set();

        if (currentGraphKey.includes("graph1")) {
            links.forEach(link => {
                if (!link.isColo) return;

                const dc = link.diseaseCounts || {};
                if ((Number(dc[filters.disease]) || 0) > 0) {
                    diseaseFilteredIds.add(typeof link.source === "object" ? link.source.id : link.source);
                    diseaseFilteredIds.add(typeof link.target === "object" ? link.target.id : link.target);
                }
            });
        } else {
            nodes.forEach(node => {
                const dc = node.diseaseCounts || {};
                if ((Number(dc[filters.disease]) || 0) > 0) {
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

// --- CENTERING FUNCTION ---
function centerOnNodes(svg, zoom, nodes, width, height) {
  if (!nodes || nodes.length === 0) return;

  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;

  for (const n of nodes) {
    if (!isFinite(n.x) || !isFinite(n.y)) continue;
    if (n.x < minX) minX = n.x;
    if (n.y < minY) minY = n.y;
    if (n.x > maxX) maxX = n.x;
    if (n.y > maxY) maxY = n.y;
  }

  if (!isFinite(minX)) return;

  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;

  // get current zoom scale, keep it
  const t = d3.zoomTransform(svg.node());
  const k = t.k;

  const tx = width / 2 - cx * k;
  const ty = height / 2 - cy * k;

  svg.transition().duration(350).call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(k));
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
            .attr("fill", c)
            .attr("stroke", c)
            .attr("stroke-width", 0.5);
    });

    const simNodes = data.nodes.map(d => ({ ...d }));
    const simLinks = data.links.map(d => ({ ...d }));

    const neighborMap = new Map(simNodes.map(n => [n.id, new Set()]));

    simLinks.forEach(l => {
      const s = typeof l.source === "object" ? l.source.id : l.source;
      const t = typeof l.target === "object" ? l.target.id : l.target;
      neighborMap.get(s)?.add(t);
      neighborMap.get(t)?.add(s);
    });

    function getDraggedGroup(node) {
      const ids = new Set([node.id, ...(neighborMap.get(node.id) || [])]);
      return simNodes.filter(n => ids.has(n.id));
    }

    function dragstart(event, d) {
      d._dragGroup = getDraggedGroup(d);
      d._lastX = event.x;
      d._lastY = event.y;
    }

    function dragged(event, d) {
      if (!d._dragGroup) return;

      const dx = event.x - d._lastX;
      const dy = event.y - d._lastY;

      d._dragGroup.forEach(n => {
        n.x += dx;
        n.y += dy;
      });

      d._lastX = event.x;
      d._lastY = event.y;

      ticked();
    }

    function dragend(event, d) {
      d._dragGroup = null;
    }

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
            if (!d.isColo) return;  // Only colocalization links

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
                    return 200;
                }
                const count = (typeof d.patientCount === "number") ? d.patientCount : 0;

                // Smooth scaling between 150 and 4000 for counts 1–30
                const MIN_COUNT = 1;
                const MAX_COUNT = 60;
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
        .attr("stroke", d => "#999")
        .attr("stroke-width", 0.5)
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

    
    function getCenter() {
      const vb = svg.node().viewBox.baseVal;   // {x,y,width,height}
      return [vb.x + vb.width / 2, vb.y + vb.height / 2];
    }
    const [cx, cy] = getCenter();

    // --- simulation ---
    const sim = d3.forceSimulation(simNodes)
      .force("link", d3.forceLink(simLinks).id(d => d.id).distance(d => d.isColo ? 40 : 60))
      .force("charge", d3.forceManyBody().strength(d => -(40 + (d.degree || 0) * 15)))
      .force("collision", d3.forceCollide().radius(20))
      .force("center", d3.forceCenter(cx, cy))
      .force("x", d3.forceX(cx).strength(0.05))
      .force("y", d3.forceY(cy).strength(0.05))
      .alpha(1)
      .alphaDecay(0.03)
      .on("tick", ticked)
      .on("end", () => {
        console.log("simulation finished");
      });;

    function ticked() {
      linkSelection.attr("d", d => linkArc(d));
      nodeSelection.attr("transform", d => `translate(${d.x},${d.y})`);
      labelSelection.attr("x", d => d.x).attr("y", d => d.y);
    }

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

// ----ARG, MGE Search suggestions population based on current visible nodes after filtering---
function populateSearchSuggestionsFromNodes(nodes) {
  const argList = document.getElementById("argSuggestions");
  const mgeList = document.getElementById("mgeSuggestions");
  const argInput = document.getElementById("argSearch");
  const mgeInput = document.getElementById("mgeSearch");

  if (!argList || !mgeList || !argInput || !mgeInput) return;

  argList.innerHTML = "";
  mgeList.innerHTML = "";

  const argQuery = argInput.value.trim().toLowerCase();
  const mgeQuery = mgeInput.value.trim().toLowerCase();

  function rankAndLimit(labels, query) {
    const uniqueLabels = [...new Set(labels.map(s => s.trim()).filter(Boolean))];

    if (!query) {
      return uniqueLabels
        .sort((a, b) => a.localeCompare(b))
        .slice(0, 5);
    }

    const prefixMatches = [];
    const containsMatches = [];

    uniqueLabels.forEach(label => {
      const lower = label.toLowerCase();

      if (lower.startsWith(query)) {
        prefixMatches.push(label);
      } else if (lower.includes(query)) {
        containsMatches.push(label);
      }
    });

    prefixMatches.sort((a, b) => a.localeCompare(b));
    containsMatches.sort((a, b) => a.localeCompare(b));

    return [...prefixMatches, ...containsMatches].slice(0, 5);
  }

  const argLabels = rankAndLimit(
    nodes.filter(n => n.isARG && n.label).map(n => n.label),
    argQuery
  );

  const mgeLabels = rankAndLimit(
    nodes.filter(n => !n.isARG && n.label).map(n => n.label),
    mgeQuery
  );

  argLabels.forEach(label => {
    const option = document.createElement("option");
    option.value = label;
    argList.appendChild(option);
  });

  mgeLabels.forEach(label => {
    const option = document.createElement("option");
    option.value = label;
    mgeList.appendChild(option);
  });
}


// --- SVG DOWNLOAD ---
function downloadCurrentGraphAsSVG() {
  console.log("SVG clicked");
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

// --- PDF DOWNLOAD ---
async function downloadCurrentGraphAsPDF() {
  const svgEl = document.querySelector("#graph");        // your <svg id="graph">
  if (!svgEl) return;

  // ---- Clone SVG ----
  const clone = svgEl.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");

  // ---- Inline CSS ----
  const cssStyles = Array.from(document.styleSheets)
    .map(ss => {
      try { return Array.from(ss.cssRules).map(r => r.cssText).join("\n"); }
      catch { return ""; }
    })
    .join("\n");
  const styleEl = document.createElement("style");
  styleEl.textContent = cssStyles;
  clone.insertBefore(styleEl, clone.firstChild);

  // Use the on-screen size of the SVG element
  const rect = svgEl.getBoundingClientRect();
  const outW = Math.max(1, Math.round(rect.width));
  const outH = Math.max(1, Math.round(rect.height));

  // Match exported pixel size to what user sees
  clone.setAttribute("width", outW);
  clone.setAttribute("height", outH);

  // Ensure viewBox corresponds to current viewBox (your SVG already has it)
  // If you don't have one, set it:
  if (!clone.getAttribute("viewBox")) {
    clone.setAttribute("viewBox", `0 0 ${outW} ${outH}`);
  }

  // ---- Serialize ----
  const svgString = new XMLSerializer().serializeToString(clone);

  // ---- SVG -> Canvas ----
  const svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
  const svgUrl = URL.createObjectURL(svgBlob);

  const img = new Image();
  await new Promise((resolve, reject) => {
    img.onload = resolve;
    img.onerror = () => reject(new Error("SVG->IMG load failed (CORS/fonts/images?)"));
    img.src = svgUrl;
  });

  const scale = 2; // sharper
  const canvas = document.createElement("canvas");
  canvas.width = outW * scale;
  canvas.height = outH * scale;

  const ctx = canvas.getContext("2d");
  ctx.setTransform(scale, 0, 0, scale, 0, 0);
  ctx.drawImage(img, 0, 0, outW, outH);

  URL.revokeObjectURL(svgUrl);

  // ---- Canvas -> PDF ----
  const { jsPDF } = window.jspdf;
  const pdf = new jsPDF({
    orientation: outW >= outH ? "landscape" : "portrait",
    unit: "px",
    format: [outW, outH]
  });

  pdf.addImage(canvas.toDataURL("image/png"), "PNG", 0, 0, outW, outH);

  const fileName =
    (typeof currentGraphKey === "string"
      ? currentGraphKey.replace(/^.*[\\/]/, "").replace(".json", "")
      : "graph") + ".pdf";

  pdf.save(fileName);
}

// --- Text DOWNLOAD ---
function downloadCurrentGraphAsText() {
  console.log("CSV clicked");

  if (!currentRenderedData || !currentRenderedData.nodes.length) {
    console.warn("No rendered data available for export");
    return;
  }

  const activeDisease = d3.select("#diseaseFilter").property("value");

  const nodes = currentRenderedData.nodes;
  const links = currentRenderedData.links || [];

  const rows = [];
  rows.push([
    "ARG_Label",
    "MGE_Label",
    "MGE_Group",
    "Disease",
    "Donor",
    "Pre",
    "Post",
    "PatientCount"
  ]);

  // -----------------------------
  // graph2: nodes already represent ARG-MGE colocalizations
  // -----------------------------
  if (currentGraphKey.includes("graph2.json")) {
    for (const node of nodes) {
      const tp = Number(node.timepoint);
      const category = getTimepointCategory(tp);

      const diseaseCounts = node.diseaseCounts || {};

      let diseases = [];
      let patientCount = 0;

      if (activeDisease !== "all") {
        patientCount = Number(diseaseCounts[activeDisease]) || 0;
        if (patientCount <= 0) continue; // only export selected disease
        diseases = [activeDisease];
      } else {
        diseases = Object.keys(diseaseCounts).filter(
          d => Number(diseaseCounts[d]) > 0
        );
        patientCount = Object.values(diseaseCounts).reduce(
          (sum, v) => sum + (Number(v) || 0),
          0
        );
      }

      const donor = category === "donor" ? patientCount : 0;
      const pre = category === "pre" ? patientCount : 0;
      const post =
        category === "post1" || category === "post2" || category === "post3"
          ? patientCount
          : 0;

      let argLabel = "";
      let mgeLabel = node.label || "";

      if (typeof node.label === "string" && node.label.includes("+")) {
        const parts = node.label.split("+");
        argLabel = parts[0]?.trim() || "";
        mgeLabel = parts.slice(1).join("+").trim() || "";
      }

      rows.push([
        argLabel,
        mgeLabel,
        node.mgeGroup || "",
        diseases.join(";"),
        donor,
        pre,
        post,
        patientCount
      ]);
    }
  }

  // -----------------------------
  // graph1: export visible colocalization links
  // -----------------------------
  else {
    const nodeById = new Map(nodes.map(n => [n.id, n]));
    const coloLinks = links.filter(l => l.isColo);

    for (const link of coloLinks) {
      const sourceId = typeof link.source === "object" ? link.source.id : link.source;
      const targetId = typeof link.target === "object" ? link.target.id : link.target;

      const sourceNode = nodeById.get(sourceId);
      const targetNode = nodeById.get(targetId);
      if (!sourceNode || !targetNode) continue;

      let argNode, mgeNode;

      if (sourceNode.isARG && !targetNode.isARG) {
        argNode = sourceNode;
        mgeNode = targetNode;
      } else if (!sourceNode.isARG && targetNode.isARG) {
        argNode = targetNode;
        mgeNode = sourceNode;
      } else {
        continue;
      }

      const dc = link.diseaseCounts || {};

      let disease = "";
      let patientCount = 0;

      if (activeDisease !== "all") {
        patientCount = Number(dc[activeDisease]) || 0;
        if (patientCount <= 0) continue; // only export selected disease
        disease = activeDisease;
      } else {
        disease = link.diseases ? link.diseases.join(";") : "";
        patientCount = Object.values(dc).reduce(
          (sum, v) => sum + (Number(v) || 0),
          0
        );

        // fallback if diseaseCounts missing
        if (patientCount === 0) {
          patientCount = Number(link.patientCount ?? link.individualCount ?? 0);
        }
      }

      const tp = Number(link.sourceTimepoint ?? link.targetTimepoint);
      const category = getTimepointCategory(tp);

      const donor = category === "donor" ? patientCount : 0;
      const pre = category === "pre" ? patientCount : 0;
      const post =
        category === "post1" || category === "post2" || category === "post3"
          ? patientCount
          : 0;

      rows.push([
        argNode.label || "",
        mgeNode.label || "",
        mgeNode.mgeGroup || "",
        disease,
        donor,
        pre,
        post,
        patientCount
      ]);
    }
  }

  const csvContent = rows
    .map(row =>
      row
        .map(value => {
          const s = String(value ?? "");
          return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
        })
        .join(",")
    )
    .join("\n");

  const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);

  const a = document.createElement("a");
  const fileName = currentGraphKey.replace(/^.*[\\/]/, "").replace(".json", "") + ".csv";
  a.href = url;
  a.download = fileName;

  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

// --- GRAPH STATISTICS ---
function computeGraphStats(nodes, links) {
  const nodeCount = nodes.length;

  const argCount = nodes.reduce((acc, n) => acc + (n.isARG ? 1 : 0), 0);
  const mgeCount = nodeCount - argCount;

  const edgeCount = links.length;

  // Based on link objects: type: "temporal" | "colocalization"
  const colocCount = links.filter(e => e.type === "colocalization" || e.isColo).length;
  const temporalCount = links.filter(e => e.type === "temporal").length;
//   const colocCount = links.reduce((acc, e) => acc + (e.type === "colocalization" ? 1 : 0), 0);
//   const temporalCount = links.reduce((acc, e) => acc + (e.type === "temporal" ? 1 : 0), 0);

  return { nodeCount, argCount, mgeCount, edgeCount, colocCount, temporalCount };
}



function renderGraphStats(stats) {
  document.getElementById("st-nodes").textContent = stats.nodeCount;
  document.getElementById("st-arg").textContent = stats.argCount;
  document.getElementById("st-mge").textContent = stats.mgeCount;

  document.getElementById("st-edges").textContent = stats.edgeCount;
  document.getElementById("st-coloc").textContent = stats.colocCount;
  document.getElementById("st-temp").textContent = stats.temporalCount;
}

function updateGraphStatsVisibility() {
  const dataset = document.getElementById("dataset")?.value || "";
  const isColoc = dataset.includes("graph2");

  // elements
  const show = (id, v) =>
    document.getElementById(id)?.classList.toggle("d-none", !v);

  if (isColoc) {
    // Only Nodes + Edges
    show("stat-nodes", true);
    show("stat-edges", true);

    show("stat-arg", false);
    show("stat-mge", false);
    show("stat-coloc", false);
    show("stat-temp", false);
    show("stat-divider", false);
  } else {
    // Full stats
    show("stat-nodes", true);
    show("stat-edges", true);
    show("stat-arg", true);
    show("stat-mge", true);
    show("stat-coloc", true);
    show("stat-temp", true);
    show("stat-divider", true);
  }
}

function restrictLinksToVisibleNodesRobust(visibleNodes, allLinks) {
  const visible = new Set(visibleNodes.map(n => n.id));

  return allLinks.filter(l => {
    const s = typeof l.source === "object" ? l.source.id : l.source;
    const t = typeof l.target === "object" ? l.target.id : l.target;
    return visible.has(s) && visible.has(t);
  });
}

const fab = document.getElementById("zoom-fab");
const graphBox = document.getElementById("graph-container");


// --- FLOATING ACTION BUTTON POSITIONING TO BOTTOM---
function positionFab() {
  if (!fab || !graphBox) return;

  const r = graphBox.getBoundingClientRect();

  const padding = 16; // inside the graph box
  const fabWidth = fab.offsetWidth;
  const fabHeight = fab.offsetHeight;

  // Target: bottom-right INSIDE the graph box
  let left = r.right - fabWidth - padding;
  let top  = r.bottom - fabHeight - padding;

  // Clamp so it never goes off the viewport
  const minLeft = 8;
  const minTop = 8;
  const maxLeft = window.innerWidth - fabWidth - 8;
  const maxTop = window.innerHeight - fabHeight - 8;

  left = Math.max(minLeft, Math.min(left, maxLeft));
  top  = Math.max(minTop, Math.min(top, maxTop));

  fab.style.left = `${left}px`;
  fab.style.top = `${top}px`;

  // Optional: hide if graph is completely off-screen
  const visible =
    r.bottom > 0 &&
    r.top < window.innerHeight &&
    r.right > 0 &&
    r.left < window.innerWidth;

  fab.style.display = visible ? "flex" : "none";
}



document.addEventListener("DOMContentLoaded", function () {
  const popoverTriggerList = document.querySelectorAll('[data-bs-toggle="popover"]');
  const popoverList = [...popoverTriggerList].map(el => {
    return new bootstrap.Popover(el, {
      container: 'body'
    });
  });

  // Optional: keep only one popover open at a time
  document.querySelectorAll('.help-icon').forEach(btn => {
    btn.addEventListener('click', function () {
      popoverList.forEach(pop => {
        const triggerEl = pop._element;
        if (triggerEl !== btn) {
          pop.hide();
        }
      });
    });
  });
});



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
document.getElementById("downloadPdfBtn").addEventListener("click", downloadCurrentGraphAsPDF);
document.getElementById("downloadDataBtn").addEventListener("click", downloadCurrentGraphAsText);
document.getElementById("argSearch")?.addEventListener("input", () => {populateSearchSuggestionsFromNodes(currentRenderedData?.nodes || []);});
document.getElementById("mgeSearch")?.addEventListener("input", () => {populateSearchSuggestionsFromNodes(currentRenderedData?.nodes || []);});

// Update on scroll + resize + initial load
window.addEventListener("scroll", positionFab, { passive: true });
window.addEventListener("resize", positionFab);
positionFab();

const legendOverlay = document.getElementById("legendOverlay");
// document.getElementById("toggleLegend").addEventListener("change", function() {
// 	if (this.checked) {
// 		legendOverlay.classList.add("visible");
// 	} else {
// 		legendOverlay.classList.remove("visible");
// 	}
// });

// --- INITIAL LOAD ---
bindTimepointListeners();
loadAndRenderGraph(currentGraphKey);
updateTimepointButtonText();
updateGraphStatsVisibility();

