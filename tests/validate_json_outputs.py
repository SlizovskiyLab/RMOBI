#!/usr/bin/env python3
"""Validate JSON outputs consumed by the static web app."""

# High-level CI guard: confirm generated backend JSON keeps the graph and
# temporal-dynamics fields that the RMOBI web interface expects to render.

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
JSON_DIR = ROOT / "docs" / "json"
REQUIRED_DISEASES = {"MDRB", "Melanoma", "rCDI"}


def fail(message: str) -> None:
    raise AssertionError(message)


def require_keys(obj: dict[str, Any], keys: set[str], context: str) -> None:
    missing = sorted(keys - obj.keys())
    if missing:
        fail(f"{context} missing required field(s): {', '.join(missing)}")


def validate_graph(path: Path, data: dict[str, Any]) -> None:
    require_keys(data, {"nodes"}, path.name)
    edge_key = "edges" if "edges" in data else "links" if "links" in data else None
    if edge_key is None:
        fail(f"{path.name} missing required edge collection: expected 'edges' or 'links'")

    nodes = data["nodes"]
    edges = data[edge_key]
    if not isinstance(nodes, list) or not nodes:
        fail(f"{path.name} field 'nodes' must be a non-empty list")
    if not isinstance(edges, list) or not edges:
        fail(f"{path.name} field '{edge_key}' must be a non-empty list")

    node_ids: set[str] = set()
    for index, node in enumerate(nodes):
        if not isinstance(node, dict):
            fail(f"{path.name} nodes[{index}] must be an object")
        require_keys(node, {"id", "timepoint", "mgeGroup"}, f"{path.name} nodes[{index}]")

        node_id = node["id"]
        if not isinstance(node_id, str) or not node_id:
            fail(f"{path.name} nodes[{index}].id must be a non-empty string")
        if node_id in node_ids:
            fail(f"{path.name} has duplicate node id: {node_id}")
        node_ids.add(node_id)

        if not isinstance(node["timepoint"], int):
            fail(f"{path.name} node {node_id} timepoint must be an integer")
        if "MGE" in node_id and not node["mgeGroup"]:
            fail(f"{path.name} node {node_id} is an MGE but has no mgeGroup")

    seen_diseases: set[str] = set()
    for index, edge in enumerate(edges):
        if not isinstance(edge, dict):
            fail(f"{path.name} {edge_key}[{index}] must be an object")
        require_keys(
            edge,
            {"source", "target", "sourceTimepoint", "targetTimepoint", "diseases"},
            f"{path.name} {edge_key}[{index}]",
        )

        source = edge["source"]
        target = edge["target"]
        if source not in node_ids:
            fail(f"{path.name} {edge_key}[{index}] source not found in nodes: {source}")
        if target not in node_ids:
            fail(f"{path.name} {edge_key}[{index}] target not found in nodes: {target}")
        if not isinstance(edge["sourceTimepoint"], int):
            fail(f"{path.name} {edge_key}[{index}].sourceTimepoint must be an integer")
        if not isinstance(edge["targetTimepoint"], int):
            fail(f"{path.name} {edge_key}[{index}].targetTimepoint must be an integer")
        if not isinstance(edge["diseases"], list) or not edge["diseases"]:
            fail(f"{path.name} {edge_key}[{index}].diseases must be a non-empty list")
        seen_diseases.update(str(disease) for disease in edge["diseases"])

    missing_diseases = sorted(REQUIRED_DISEASES - seen_diseases)
    if missing_diseases:
        fail(f"{path.name} does not include disease(s): {', '.join(missing_diseases)}")


def validate_temporal_dynamics(path: Path, data: dict[str, Any]) -> None:
    missing_diseases = sorted(REQUIRED_DISEASES - data.keys())
    if missing_diseases:
        fail(f"{path.name} missing disease group(s): {', '.join(missing_diseases)}")

    required_fields = {"colocalization", "mgeGroup", "status", "patients"}
    for disease, records in data.items():
        if not isinstance(records, list) or not records:
            fail(f"{path.name} disease {disease} must contain a non-empty list")
        for index, record in enumerate(records):
            if not isinstance(record, dict):
                fail(f"{path.name} {disease}[{index}] must be an object")
            require_keys(record, required_fields, f"{path.name} {disease}[{index}]")
            if not record["mgeGroup"]:
                fail(f"{path.name} {disease}[{index}].mgeGroup must be non-empty")
            if not isinstance(record["patients"], int) or record["patients"] < 0:
                fail(f"{path.name} {disease}[{index}].patients must be a non-negative integer")


def validate_json_file(path: Path) -> None:
    with path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        fail(f"{path.name} must contain a JSON object at the top level")

    if "nodes" in data or "edges" in data or "links" in data:
        validate_graph(path, data)
    else:
        validate_temporal_dynamics(path, data)


def main() -> int:
    json_files = sorted(JSON_DIR.glob("*.json"))
    if not json_files:
        fail(f"No JSON outputs found in {JSON_DIR}")

    for path in json_files:
        validate_json_file(path)
        print(f"validated {path.relative_to(ROOT)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"JSON validation failed: {error}", file=sys.stderr)
        raise SystemExit(1)
