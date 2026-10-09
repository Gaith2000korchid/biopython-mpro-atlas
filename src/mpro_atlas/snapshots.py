"""Record successful UTF-8 API payloads and replay them without network access."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path

import requests


class RecordingSession:
    retrieval_mode = "live_recorded"

    def __init__(self, session, directory):
        self.session = session
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        if (self.directory / "manifest.json").exists():
            raise ValueError("Snapshot already exists; choose a new output directory")
        self.manifest = {"format_version": 1, "encoding": "utf-8",
                         "sources": {}}

    def get(self, url, timeout):
        response = self.session.get(url, timeout=timeout)
        response.raise_for_status()
        # Store the decoded text used by parsers; these are not original wire bytes.
        payload = response.text.encode("utf-8")
        filename = sha256(url.encode("utf-8")).hexdigest() + ".txt"
        (self.directory / filename).write_bytes(payload)
        self.manifest["sources"][url] = {
            "file": filename, "sha256": sha256(payload).hexdigest(),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        (self.directory / "manifest.json").write_text(
            json.dumps(self.manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return response


class ReplaySession:
    retrieval_mode = "offline_snapshot"

    def __init__(self, directory):
        self.directory = Path(directory)
        self.manifest = json.loads(
            (self.directory / "manifest.json").read_text(encoding="utf-8")
        )
        if self.manifest.get("format_version") != 1:
            raise ValueError("Unsupported snapshot format")

    def get(self, url, timeout):
        del timeout
        source = self.manifest["sources"].get(url)
        if source is None:
            raise ValueError(f"URL is absent from offline snapshot: {url}")
        filename = sha256(url.encode("utf-8")).hexdigest() + ".txt"
        if source["file"] != filename:
            raise ValueError("Invalid snapshot filename")
        payload = (self.directory / filename).read_bytes()
        if sha256(payload).hexdigest() != source["sha256"]:
            raise ValueError("Snapshot payload checksum mismatch")
        response = requests.Response()
        response.status_code = 200
        response.url = url
        response.encoding = "utf-8"
        response._content = payload
        return response
