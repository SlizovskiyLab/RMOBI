import csv
from collections import OrderedDict
from typing import Dict, List

csv_file_path = "../../data/patientwise_colocalization_by_timepoint.csv"
output_file_name = '../../scripts/mge_id.txt'

def create_unique_group_id_map(
    csv_file_path: str,
    output_map_txt_path: str,
    group_col: str = "MGE gene",
    start_id: int = 1001,
    sort_groups: bool = True,
    include_header_comment: bool = False,
) -> Dict[str, int]:
    """
    1) Reads the CSV and assigns a unique integer ID to each distinct value in `group_col`.
    2) Writes:
       - a txt file with {id, "Name"} lines
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

    with open(output_map_txt_path, "w", encoding="utf-8") as out:
        for g, i in mapping.items():
            out.write(f'{{{i}, "{g}"}},\n')


    return dict(mapping)

create_unique_group_id_map(csv_file_path,output_file_name)
