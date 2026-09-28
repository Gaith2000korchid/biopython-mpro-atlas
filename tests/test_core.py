import pytest
from mpro_atlas.core import (
    domain_range_from_uniprot,
    domain_structures_csv,
    filter_domain_structures,
    inclusive_overlap,
    map_pdb_to_domain,
    parse_aux,
    parse_chain_ranges,
    parse_cif,
    parse_resolution_angstrom,
    parse_uniprot,
    build_report,
)


def test_aux_separates_profile_structures_from_protein_references():
    aux = parse_aux(
        "ID   M_PRO; MATRIX.\nAC   PS51442;\nNR   /RELEASE=2026_03,575748;\n"
        "DR   P0DTD1    , R1AB_SARS2 , T; P0C6X7    , R1AB_SARS  , T;\n"
        "3D   6LU7; 1LVO;\n3D   6LU7; 7NH7;\n"
    )
    assert aux == {
        "accession": "PS51442",
        "release": "2026_03",
        "pdb_ids": ["1LVO", "6LU7", "7NH7"],
        "uniprot_accessions": ["P0C6X7", "P0DTD1"],
    }


def test_aux_rejects_tokens_that_are_not_pdb_ids():
    aux = parse_aux("AC   PS51442;\n3D   NOTE; 6LU7; ABCD;\nDR   P0DTD1    , R1AB_SARS2 , T;\n")
    assert aux["pdb_ids"] == ["6LU7"]


def test_aux_requires_accession_and_annotation_lines():
    with pytest.raises(ValueError, match="Missing PROSITE accession"):
        parse_aux("ID   M_PRO; MATRIX.\n3D   6LU7;\n")
    with pytest.raises(ValueError, match="no 3D or DR"):
        parse_aux("AC   PS51442;\nNR   /RELEASE=2026_03,1;\n")


def test_uniprot_properties_by_key_not_order():
    entry = {
        "primaryAccession": "P0DTD1",
        "sequence": {"value": "ACD"},
        "features": [{
            "type": "Domain",
            "description": "3C-like proteinase nsp5",
            "location": {"start": {"value": 3264}, "end": {"value": 3569}},
            "featureCrossReferences": [{"database": "PROSITE", "id": "PS51442"}],
        }],
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "6LU7", "properties": [
                {"key": "Chains", "value": "A=3264-3569"}, {"key": "Method", "value": "X-ray"}]},
            {"database": "PROSITE", "id": "PS51442"},
        ],
    }
    result = parse_uniprot(entry)
    assert result["pdb"]["6LU7"]["Method"] == "X-ray"
    assert result["pdb"]["6LU7"]["Chains"] == "A=3264-3569"
    assert result["prosite"] == ["PS51442"]
    assert result["domain"]["start"] == 3264 and result["domain"]["source"] == "uniprot_feature_prosite"


def test_domain_falls_back_to_description_keywords():
    domain = domain_range_from_uniprot({"features": [{
        "type": "Chain", "description": "3C-like proteinase",
        "location": {"start": {"value": 3264}, "end": {"value": 3569}},
    }]})
    assert domain["source"] == "uniprot_feature_description" and domain["start"] == 3264


def test_chain_range_parsing_and_overlap():
    assert parse_chain_ranges("A=3264-3569") == [{"chain": "A", "start": 3264, "end": 3569}]
    assert parse_chain_ranges("A/B=3264-3569") == [
        {"chain": "A", "start": 3264, "end": 3569}, {"chain": "B", "start": 3264, "end": 3569}]
    assert inclusive_overlap(3264, 3569, 1, 932) == 0
    assert inclusive_overlap(3264, 3569, 3264, 3569) == 306


def test_map_keeps_mpro_and_drops_nonoverlapping_polyprotein_chains():
    mpro = map_pdb_to_domain("6LU7", {"Chains": "A=3264-3569", "Resolution": "2.16 A"}, 3264, 3569)
    polymerase = map_pdb_to_domain("7LNN", {"Chains": "A=1-932"}, 3264, 3569)
    assert mpro["covers_domain"] and mpro["coverage_pct"] == 100.0
    assert not polymerase["covers_domain"]


def test_filter_sorts_by_coverage_then_resolution():
    annotations = {
        "9ZZZ": {"Chains": "A=3264-3400", "Resolution": "1.50 A"},
        "6LU7": {"Chains": "A=3264-3569", "Resolution": "2.16 A"},
        "1ABC": {"Chains": "A=3264-3569", "Resolution": "1.80 A"},
        "7LNN": {"Chains": "A=1-932", "Resolution": "2.90 A"},
    }
    assert [row["pdb_id"] for row in filter_domain_structures(annotations, 3264, 3569, min_coverage=50)] == ["1ABC", "6LU7"]


def test_csv_escapes_commas_in_ligands():
    csv = domain_structures_csv([map_pdb_to_domain("6LU7", {"Chains": "A=3264-3569", "Ligands": "N3, inhibitor"}, 3264, 3569)])
    assert '"N3, inhibitor"' in csv


def test_parse_resolution_and_cif_standard_residues():
    assert parse_resolution_angstrom("2.16 A") == 2.16
    cif = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 20.00 1 ALA A CA 1
ATOM   2 N N  . ALA A 1 1 ? 0.500 0.000 0.000 1.00 20.00 1 ALA A N  1
ATOM   3 C CA . GLY A 1 2 ? 3.800 0.000 0.000 1.00 20.00 2 GLY A CA 1
HETATM 4 C C1 . N3  C 2 1 ? 6.000 0.000 0.000 1.00 20.00 1 N3  C C1 1
"""
    parsed = parse_cif(cif, "TEST")
    by_id = {chain["id"]: chain["residue_count"] for chain in parsed["chains"]}
    assert by_id["A"] == 2 and by_id.get("C", 0) == 0


class FakeResponse:
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
            return FakeResponse("AC   PS51442;\nDR   P0DTD1    , R1AB_SARS2 , T;\n3D   6LU7; 1LVO; 7LNN;\n")
        if url.endswith(".json"):
            return FakeResponse(data={
                "primaryAccession": "P0DTD1", "sequence": {"value": "ACDE"},
                "features": [{"type": "Domain", "description": "3C-like proteinase",
                              "location": {"start": {"value": 3264}, "end": {"value": 3569}},
                              "featureCrossReferences": [{"database": "PROSITE", "id": "PS51442"}]}],
                "uniProtKBCrossReferences": [
                    {"database": "PDB", "id": "6LU7", "properties": [{"key": "Chains", "value": "A=3264-3569"}]},
                    {"database": "PDB", "id": "7LNN", "properties": [{"key": "Chains", "value": "A=1-932"}]},
                    {"database": "PDB", "id": "9ZZZ", "properties": [{"key": "Chains", "value": "A=3264-3569"}]},
                ],
            })
        raise AssertionError(url)


def test_build_report_offline_filters_to_domain_and_rejects_nonoverlapping_choice():
    report, fasta, table = build_report(FakeSession(), include_structure=False)
    assert report["shared_pdb_ids"] == ["6LU7", "7LNN"]
    assert report["selected_pdb"] == "6LU7"
    assert report["domain_structure_count"] == 1
    assert "7LNN" not in table and "6LU7" in table
    with pytest.raises(ValueError, match="does not overlap"):
        build_report(FakeSession(), pdb_id="7LNN", include_structure=False)
    with pytest.raises(ValueError, match="not cross-referenced"):
        build_report(FakeSession(), pdb_id="1LVO", include_structure=False)
