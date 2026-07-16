#!/usr/bin/env python3
"""Serve docs locally and verify the static site plus JSON assets load."""

# CI guard: Apply the GitHub Pages-style static app over HTTP so broken page paths, missing data files, and invalid JSON responses fail fast.

from __future__ import annotations

import contextlib
import functools
import http.server
import json
import socket
import socketserver
import sys
import threading
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return


def get_free_port() -> int:
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def fetch(base_url: str, path: str) -> bytes:
    url = f"{base_url}/{path.lstrip('/')}"
    with urllib.request.urlopen(url, timeout=10) as response:
        if response.status != 200:
            raise AssertionError(f"{url} returned HTTP {response.status}")
        return response.read()


def assert_page_loads(base_url: str, page: str, required_snippets: list[str]) -> None:
    body = fetch(base_url, page).decode("utf-8")
    for snippet in required_snippets:
        if snippet not in body:
            raise AssertionError(f"{page} missing expected content: {snippet}")
    print(f"loaded {page}")


def assert_json_loads(base_url: str, path: str) -> None:
    body = fetch(base_url, path).decode("utf-8")
    data = json.loads(body)
    if not isinstance(data, dict) or not data:
        raise AssertionError(f"{path} did not load as a non-empty JSON object")
    print(f"loaded {path}")


def main() -> int:
    if not (DOCS / "index.html").exists():
        raise AssertionError(f"Missing docs/index.html at {DOCS}")

    port = get_free_port()
    handler = functools.partial(QuietHandler, directory=str(DOCS))
    server = socketserver.TCPServer(("127.0.0.1", port), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        base_url = f"http://127.0.0.1:{port}"
        assert_page_loads(base_url, "index.html", ["RMOBI"])
        assert_page_loads(base_url, "tool.html", ["json/graph1.json", "json/graph2.json"])
        assert_page_loads(base_url, "analysis.html", ["json/temporal_dynamics_disease.json"])
        assert_json_loads(base_url, "json/graph1.json")
        assert_json_loads(base_url, "json/graph2.json")
        assert_json_loads(base_url, "json/temporal_dynamics_disease.json")
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Web smoke test failed: {error}", file=sys.stderr)
        raise SystemExit(1)
