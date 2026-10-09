import json
import pytest
import requests

from mpro_atlas.core import build_report
from mpro_atlas.snapshots import RecordingSession, ReplaySession
from test_core import FakeSession


def test_recorded_report_can_be_replayed_without_a_network_session(tmp_path):
    sources = tmp_path / "sources"
    live, fasta, table = build_report(
        RecordingSession(FakeSession(), sources), include_structure=False
    )
    replay, replay_fasta, replay_table = build_report(
        ReplaySession(sources), include_structure=False
    )
    assert replay["retrieval_mode"] == "offline_snapshot"
    for key in ["profile", "domain", "coverage_tiers", "selected_sifts",
                "domain_structures", "provenance"]:
        assert replay[key] == live[key]
    assert replay_table == table and str(replay_fasta.seq) == str(fasta.seq)


def test_replay_rejects_modified_payload_and_missing_url(tmp_path):
    recorder = RecordingSession(FakeSession(), tmp_path)
    url = "https://example.test/PS51442.txt"
    recorder.get(url, timeout=1)
    replay = ReplaySession(tmp_path)
    with pytest.raises(ValueError, match="absent"):
        replay.get("https://example.test/unknown.txt", timeout=1)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    (tmp_path / manifest["sources"][url]["file"]).write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        replay.get(url, timeout=1)
    with pytest.raises(ValueError, match="already exists"):
        RecordingSession(FakeSession(), tmp_path)


def test_failed_http_response_is_not_archived(tmp_path):
    class FailedSession:
        def get(self, url, timeout):
            response = requests.Response()
            response.status_code = 503
            response.url = url
            return response
    with pytest.raises(requests.HTTPError):
        RecordingSession(FailedSession(), tmp_path).get("https://example.test", 1)
    assert not (tmp_path / "manifest.json").exists()
