import pytest
from mpro_atlas.core import parse_aux, parse_uniprot, build_report


def test_aux_separates_profile_structures_from_protein_references():
    aux = parse_aux("ID   M_PRO; MATRIX.\nAC   PS51442;\nNR   /RELEASE=2026_03,575748;\n"
                    "DR   P0DTD1    , R1AB_SARS2 , T; P0C6X7    , R1AB_SARS  , T;\n"
                    "3D   6LU7; 1LVO;\n3D   6LU7; 7NH7;\n")
    assert aux == {"accession": "PS51442", "release": "2026_03",
                   "pdb_ids": ["1LVO", "6LU7", "7NH7"],
                   "uniprot_accessions": ["P0C6X7", "P0DTD1"]}


def test_uniprot_properties_by_key_not_order():
    entry = {"primaryAccession": "P0DTD1", "sequence": {"value": "ACD"},
             "uniProtKBCrossReferences": [
                 {"database": "PDB", "id": "6LU7", "properties": [
                     {"key": "Chains", "value": "A=3264-3569"}, {"key": "Method", "value": "X-ray"}]},
                 {"database": "PROSITE", "id": "PS51442"}]}
    result = parse_uniprot(entry)
    assert result["pdb"]["6LU7"]["Method"] == "X-ray"
    assert result["pdb"]["6LU7"]["Chains"] == "A=3264-3569"
    assert result["prosite"] == ["PS51442"]


class FakeResponse:
    text = ""
    def __init__(self, text="", data=None):
        self.text, self.data = text, data
    def raise_for_status(self):
        pass
    def json(self):
        return self.data


class FakeSession:
    def get(self, url, timeout):
        if url.endswith(".txt"):
            return FakeResponse("ID   M_PRO; MATRIX.\nAC   PS51442;\nDE   Coronavirus main protease domain.\n//\n")
        if url.endswith(".aux"):
            return FakeResponse("AC   PS51442;\nDR   P0DTD1    , R1AB_SARS2 , T;\n3D   6LU7; 1LVO;\n")
        if url.endswith(".json"):
            return FakeResponse(data={"primaryAccession": "P0DTD1", "sequence": {"value": "ACDE"},
                "uniProtKBCrossReferences": [{"database": "PDB", "id": "6LU7"},
                                               {"database": "PDB", "id": "9ZZZ"}]})
        raise AssertionError(url)


def test_build_report_offline_and_validate_selection():
    report, fasta = build_report(FakeSession(), include_structure=False)
    assert report["shared_pdb_ids"] == ["6LU7"]
    assert report["selected_pdb"] == "6LU7"
    assert report["profile"]["primary_record_pdb_count"] == 0
    assert str(fasta.seq) == "ACDE"
    with pytest.raises(ValueError, match="not cross-referenced"):
        build_report(FakeSession(), pdb_id="1LVO", include_structure=False)
