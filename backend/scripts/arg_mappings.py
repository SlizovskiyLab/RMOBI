import csv
import json
from collections import OrderedDict
from typing import Dict, List, Optional, Iterable, Tuple
import re


csv_file_path = "../../data/patientwise_colocalization_by_timepoint.csv"
megares_db_path = "../../data/megares_database_v3.00.fasta"
output_map_path = "../../data/arg_id_map.cpp"
output_map_json_path = "../../data/arg_group_id_map.json"
arg_group_cpp_path = "../../data/arg_group_id_map.cpp"
arg_resistance_group_cpp_path = "../../data/arg_resistance_group_id_map.cpp"
arg_snpconfirm_cpp_path = "../../data/arg_snpconfirm_id_map.cpp"



def create_unique_group_id_map(
    csv_file_path: str,
    output_map_txt_path: str,
    output_map_json_path: str,
    group_col: str = "MEGARes group",
    start_id: int = 1,
    sort_groups: bool = True,
    include_header_comment: bool = False,
) -> Dict[str, int]:
    """
    1) Reads the CSV and assigns a unique integer ID to each distinct value in `group_col`.
    2) Writes:
       - a txt file with {id, "Name"} lines
       - a json file with {"Name": id, ...} to be reused
    Returns: dict mapping group_name -> id
    """
    groups_seen = set()
    groups: List[str] = []

    with open(csv_file_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or group_col not in reader.fieldnames:
            raise ValueError(
                f"Column '{group_col}' not found. Available columns: {reader.fieldnames}"
            )

        for row in reader:
            val = (row.get(group_col) or "").strip()
            if not val or val.lower() in {"na", "nan", "none", "null"}:
                continue
            if val not in groups_seen:
                groups_seen.add(val)
                groups.append(val)

    if sort_groups:
        groups.sort(key=lambda s: s.upper())

    mapping: Dict[str, int] = OrderedDict((g, start_id + i) for i, g in enumerate(groups))

    # Write TXT (C-style lines)
    with open(output_map_txt_path, "w", encoding="utf-8") as out:
        out.write('const std::unordered_map<int, std::string> argIdMap = {\n')
        for g, i in mapping.items():
            out.write(f'{{{i}, "{g}"}},\n')
        out.write('};\n')

    # Write JSON for reuse
    with open(output_map_json_path, "w", encoding="utf-8") as jf:
        json.dump(mapping, jf, indent=2, sort_keys=False)

    return dict(mapping)



def load_arg_id_map(arg_id_map_path: str) -> Dict[str, int]:

    name_to_id: Dict[str, int] = {}

    # matches {123, "SOMETHING"} with optional spaces and trailing comma
    pattern = re.compile(r'\{\s*(\d+)\s*,\s*"([^"]+)"\s*\}\s*,?\s*$')

    with open(arg_id_map_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("//"):
                continue
            m = pattern.match(line)
            if not m:
                continue
            idx = int(m.group(1))
            arg_name = m.group(2).strip()
            if arg_name:
                name_to_id[arg_name] = idx


    return name_to_id




def create_arg_id_to_group_map_from_megares_fasta(
    megares_fasta_path: str,
    arg_id_map_path: str,
    output_cpp_path: str,
    group_field_index: int = 2,     # 3rd column in header: group (e.g., Aminoglycosides)
    arg_name_field_index: int = 4,  # 5th column in header: ARG name (e.g., A16S)
    drop_missing_arg_names: bool = True,
) -> Dict[int, str]:
    """
    MEGARes FASTA header example:
      >MEG_1|Drugs|Aminoglycosides|...|A16S|RequiresSNPConfirmation

    Writes:
      const std::unordered_map<int, std::string> argGroupMap = {
        {1, "Aminoglycosides"},
        {2, "Aminoglycosides"},
        ...
      };
    Returns: {arg_id: group}
    """
    arg_name_to_id = load_arg_id_map(arg_id_map_path)
    id_to_group: Dict[int, str] = {}

    def iter_fasta_headers(path: str) -> Iterable[str]:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    yield line[1:].strip()

    for hdr in iter_fasta_headers(megares_fasta_path):
        parts = [p.strip() for p in hdr.split("|")]
        if len(parts) <= max(group_field_index, arg_name_field_index):
            continue

        group = parts[group_field_index]
        arg_name = parts[arg_name_field_index]

        if not group or not arg_name:
            continue

        if arg_name not in arg_name_to_id:
            if drop_missing_arg_names:
                continue
            raise ValueError(
                f'ARG name "{arg_name}" found in FASTA (5th column) but not in arg_id_map.'
            )

        arg_id = arg_name_to_id[arg_name]
        # If duplicates appear, they should map to the same group; keep first unless conflict
        if arg_id in id_to_group and id_to_group[arg_id] != group:
            raise ValueError(
                f"Conflict for arg_id={arg_id} ({arg_name}): "
                f'"{id_to_group[arg_id]}" vs "{group}"'
            )
        id_to_group[arg_id] = group

    with open(output_cpp_path, "w", encoding="utf-8") as out:
        out.write('const std::unordered_map<int, std::string> argGroupMap = {\n')
        for arg_id in sorted(id_to_group.keys()):
            out.write(f'  {{{arg_id}, "{id_to_group[arg_id]}"}},\n')
        out.write('};\n')

    return id_to_group



def create_arg_id_to_resistance_group_map_from_megares_fasta(
    megares_fasta_path: str,
    arg_id_map_path: str,
    output_cpp_path: str,
    group_field_index: int = 1,     # 2nd column in header: group (e.g., Drugs)
    arg_name_field_index: int = 4,  # 5th column in header: ARG name (e.g., A16S)
    drop_missing_arg_names: bool = True,
) -> Dict[int, str]:
    """
    MEGARes FASTA header example:
      >MEG_1|Drugs|Aminoglycosides|...|A16S|RequiresSNPConfirmation

    Writes:
      const std::unordered_map<int, std::string> argGroupMap = {
        {1, "Drugs"},
        {2, "Drugs"},
        ...
      };
    Returns: {arg_id: group}
    """
    arg_name_to_id = load_arg_id_map(arg_id_map_path)
    id_to_group: Dict[int, str] = {}

    def iter_fasta_headers(path: str) -> Iterable[str]:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    yield line[1:].strip()

    for hdr in iter_fasta_headers(megares_fasta_path):
        parts = [p.strip() for p in hdr.split("|")]
        if len(parts) <= max(group_field_index, arg_name_field_index):
            continue

        group = parts[group_field_index]
        arg_name = parts[arg_name_field_index]

        if not group or not arg_name:
            continue

        if arg_name not in arg_name_to_id:
            if drop_missing_arg_names:
                continue
            raise ValueError(
                f'ARG name "{arg_name}" found in FASTA (5th column) but not in arg_id_map.'
            )

        arg_id = arg_name_to_id[arg_name]
        # If duplicates appear, they should map to the same group; keep first unless conflict
        if arg_id in id_to_group and id_to_group[arg_id] != group:
            raise ValueError(
                f"Conflict for arg_id={arg_id} ({arg_name}): "
                f'"{id_to_group[arg_id]}" vs "{group}"'
            )
        id_to_group[arg_id] = group

    with open(output_cpp_path, "w", encoding="utf-8") as out:
        out.write('const std::unordered_map<int, std::string> argGroupMap = {\n')
        for arg_id in sorted(id_to_group.keys()):
            out.write(f'  {{{arg_id}, "{id_to_group[arg_id]}"}},\n')
        out.write('};\n')

    return id_to_group



def create_arg_id_to_snp_conf_map_from_megares_fasta(
    megares_fasta_path: str,
    arg_id_map_path: str,
    output_cpp_path: str,
    snp_confirmation_index: int = 5,   # 6th column
    arg_name_field_index: int = 4,     # 5th column
    snp_token: str = "RequiresSNPConfirmation",
    drop_missing_arg_names: bool = True,
) -> Dict[int, bool]:
    """
    MEGARes FASTA header example:
      >MEG_1|Drugs|Aminoglycosides|...|A16S|RequiresSNPConfirmation

    Output: map ARG_ID -> bool
      true  if the 6th column exists AND equals snp_token (case-insensitive optional)
      false otherwise (including if the column is missing)
    """

    arg_name_to_id = load_arg_id_map(arg_id_map_path)

    # Default ALL known ARG IDs to False (so "if not exists false")
    id_to_snp_conf: Dict[int, bool] = {arg_id: False for arg_id in arg_name_to_id.values()}

    def iter_fasta_headers(path: str) -> Iterable[str]:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    yield line[1:].strip()

    for hdr in iter_fasta_headers(megares_fasta_path):
        parts = [p.strip() for p in hdr.split("|")]

        # make sure arg name exists
        if len(parts) <= arg_name_field_index:
            continue
        arg_name = parts[arg_name_field_index]
        if not arg_name:
            continue

        if arg_name not in arg_name_to_id:
            if drop_missing_arg_names:
                continue
            raise ValueError(
                f'ARG name "{arg_name}" found in FASTA (5th column) but not in arg_id_map.'
            )

        arg_id = arg_name_to_id[arg_name]

        # SNP confirmation is true ONLY if the field exists and matches token
        has_snp = (
            len(parts) > snp_confirmation_index
            and parts[snp_confirmation_index].strip().lower() == snp_token.lower())

        if has_snp:
            id_to_snp_conf[arg_id] = True

    # Write C++ unordered_map<int, bool>
    with open(output_cpp_path, "w", encoding="utf-8") as out:
        out.write("const std::unordered_map<int, bool> argIDSNPConfirmation = {\n")
        for arg_id in sorted(id_to_snp_conf.keys()):
            out.write(f"  {{{arg_id}, {'true' if id_to_snp_conf[arg_id] else 'false'}}},\n")
        out.write("};\n")

    return id_to_snp_conf





# -------------------------
# 1) Build and save group IDs from CSV
create_unique_group_id_map(csv_file_path, output_map_path, output_map_json_path, group_col="MEGARes group", start_id=1,sort_groups=True)

# 2) Use the saved IDs to build ID->GROUP map from MEGARes FASTA
create_arg_id_to_group_map_from_megares_fasta(megares_db_path,output_map_path,arg_group_cpp_path)
create_arg_id_to_resistance_group_map_from_megares_fasta(megares_db_path,output_map_path,arg_resistance_group_cpp_path)
create_arg_id_to_snp_conf_map_from_megares_fasta(megares_db_path,output_map_path,arg_snpconfirm_cpp_path)
