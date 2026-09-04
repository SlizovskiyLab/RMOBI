#!/usr/bin/env python3
"""Prepare SNP-filtered emergence/disappearance data for a separate R GLMM.

This script exports patient- and event-level CSV files for the R analysis.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_DATA = ROOT / "data" / "patientwise_colocalization_by_timepoint.csv"
DEFAULT_SAMPLE_METADATA = ROOT / "data" / "meta_data.csv"
DEFAULT_CLINICAL_ARGS = ROOT / "data" / "MedImpAMR_taxa.csv"
DISAPPEAR_DETAILED = ROOT / "docs" / "output" / "temporal_dynamics" / "disappear_detailed.csv"
EMERGE_DETAILED = ROOT / "docs" / "output" / "temporal_dynamics" / "emerge_detailed.csv"
PERSIST_DETAILED = ROOT / "docs" / "output" / "temporal_dynamics" / "persist_detailed.csv"
DEFAULT_OUTPUT = HERE / "output" / "model_1"
MODEL_2_OUTPUT = HERE / "output" / "model_2"
MODEL_3_OUTPUT = HERE / "output" / "model_3"
CLINICAL_MODEL_1_OUTPUT = HERE / "output" / "model_1_clinically_relevant"
TRANSITION_SCRIPT = ROOT / "analysis_scripts" / "2_transitionRate.py"
COHORTS = ["rCDI", "MDRB", "Melanoma"]
INPUT = HERE / "output" / "model_1" / "emergence_disappearance_patients.csv"
OUTPUT_DIR = HERE / "output" / "leave_one_study_out"
OUTPUT = OUTPUT_DIR / "loso_glmm_patients.csv"


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--sample-metadata", type=Path, default=DEFAULT_SAMPLE_METADATA)
    parser.add_argument("--clinical-arg-list", type=Path, default=DEFAULT_CLINICAL_ARGS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def load_transition_module():
    spec = importlib.util.spec_from_file_location("rmobi_transition_rate", TRANSITION_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {TRANSITION_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalize_patient(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def build_patient_metadata(sample_metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"Patient", "Sample", "disease_type", "study_data", "donor_pre_post"}
    if not required.issubset(sample_metadata.columns):
        raise ValueError(f"meta_data.csv is missing: {sorted(required - set(sample_metadata.columns))}")
    data = sample_metadata.copy()
    data["Patient"] = normalize_patient(data["Patient"])
    if data["Patient"].isna().any():
        raise ValueError("meta_data.csv contains an invalid Patient value.")
    consistency = data.groupby("Patient").agg(
        disease_values=("disease_type", "nunique"), study_values=("study_data", "nunique")
    )
    if (consistency[["disease_values", "study_values"]] != 1).any().any():
        raise ValueError("Each patient must map to exactly one cohort and study in meta_data.csv.")
    phase = data["donor_pre_post"].astype(str).str.lower().str.replace(r"[^a-z]", "", regex=True)
    post_counts = data.loc[phase.eq("postfmt")].groupby("Patient")["Sample"].nunique()
    metadata = data.groupby("Patient", as_index=False).agg(
        Disease_type=("disease_type", "first"), Study_ID=("study_data", "first")
    )
    metadata["Post-fmt_samples"] = metadata["Patient"].map(post_counts).fillna(0).astype(int)
    if (metadata["Post-fmt_samples"] <= 0).any():
        bad = metadata.loc[metadata["Post-fmt_samples"] <= 0, "Patient"].tolist()
        raise ValueError(f"Patients without post-FMT samples in meta_data.csv: {bad}")
    return metadata


# Count post-FMT columns with observed signal; return the inference flag too.
def infer_k_post(source: pd.DataFrame) -> pd.DataFrame:
    post_columns = [column for column in source if re.fullmatch(r"PostFMT_\d{3}", column)]
    numeric = source[post_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    observed = numeric.gt(0).groupby(source["Patient"]).any().sum(axis=1).astype(int)
    return observed.rename("k_post_detected").reset_index()


def build_contig_summaries(sample_metadata: pd.DataFrame) -> pd.DataFrame:
    """Return median assembly metrics across each patient's recipient samples."""
    required = {"Patient", "donor_receipient", "num_contigs", "total_length", "non_host_length"}
    missing = required.difference(sample_metadata.columns)
    if missing:
        raise ValueError(f"meta_data.csv is missing columns: {sorted(missing)}")
    assembly = sample_metadata.copy()
    assembly["Patient"] = normalize_patient(assembly["Patient"])
    for column in ["num_contigs", "total_length", "non_host_length"]:
        assembly[column] = pd.to_numeric(assembly[column], errors="coerce")
        if (assembly[column].dropna() <= 0).any():
            raise ValueError(f"{column} must be positive when present.")
    recipient = assembly["donor_receipient"].astype(str).str.strip().str.lower().eq("recipient")
    result = assembly.loc[recipient].groupby("Patient", as_index=False).agg(
        median_recipient_num_contigs=("num_contigs", "median"),
        median_recipient_total_length=("total_length", "median"),
        median_recipient_non_host_length=("non_host_length", "median"),
    )
    model_columns = [
        "median_recipient_num_contigs", "median_recipient_total_length",
        "median_recipient_non_host_length",
    ]
    if result[model_columns].isna().any().any():
        raise ValueError("A patient has no valid recipient assembly-metric value.")
    return result


def build_events(source: pd.DataFrame, transition) -> pd.DataFrame:
    filtered = transition.filter_to_snp_confirmed_args(source)
    phased = transition.add_phase_presence(filtered)
    events = phased.drop_duplicates(["Patient", "MEGARes group", "MGE gene"]).copy()
    donor = events["Donor_p"].astype(bool)
    pre = events["PreFMT_p"].astype(bool)
    post = events["PostAny_p"].astype(bool)
    events["event"] = np.select(
        [(~donor) & (~pre) & post, pre & (~post)],
        ["Emergence", "Disappearance"],
        default="Other",
    )
    events = events.loc[events["event"] != "Other", [
        "Patient", "Disease_type", "MEGARes group", "MGE gene", "event"
    ]].copy()
    events["Patient"] = normalize_patient(events["Patient"])
    return events.sort_values(["Patient", "event", "MEGARes group", "MGE gene"])


# def export_clinically_relevant_patients(
#     source: pd.DataFrame,
#     metadata: pd.DataFrame,
#     sample_metadata: pd.DataFrame,
#     clinical_arg_path: Path,
#     transition,
# ) -> None:
#     """Export Model 1 input restricted to clinically relevant ARG groups."""
#     clinical = pd.read_csv(clinical_arg_path)
#     if "Group" not in clinical.columns:
#         raise ValueError(f"{clinical_arg_path} must contain a 'Group' column.")
#     if "MEGARes group" not in source.columns:
#         raise ValueError("Model source must contain a 'MEGARes group' column.")
#     clinical_groups = set(clinical["Group"].dropna().astype(str).str.strip().str.upper())
#     clinical_source = source.loc[
#         source["MEGARes group"].astype(str).str.strip().str.upper().isin(clinical_groups)
#     ].copy()
#     clinical_events = build_events(clinical_source, transition)
#     clinical_patients = build_patients(
#         clinical_source, metadata, clinical_events, sample_metadata, transition
#     )
#     CLINICAL_MODEL_1_OUTPUT.mkdir(parents=True, exist_ok=True)
#     clinical_patients.to_csv(
#         CLINICAL_MODEL_1_OUTPUT / "emergence_disappearance_patients.csv", index=False
#     )
#     print(
#         f"Wrote {len(clinical_patients)} clinically relevant ARG patients to "
#         f"{CLINICAL_MODEL_1_OUTPUT / 'emergence_disappearance_patients.csv'}"
#     )


def build_patients(source: pd.DataFrame, metadata: pd.DataFrame, events: pd.DataFrame,
                   sample_metadata: pd.DataFrame, transition) -> pd.DataFrame:
    if "Post-fmt_samples" not in metadata:
        raise ValueError("Derived patient metadata must contain 'Post-fmt_samples'.")
    metadata = metadata.rename(columns={"Disease_type": "cohort", "Study_ID": "study"}).copy()
    metadata["Patient"] = normalize_patient(metadata["Patient"])
    metadata = metadata.dropna(subset=["Patient", "cohort", "study"])
    if metadata["Patient"].duplicated().any():
        raise ValueError("Derived metadata must have exactly one row per Patient.")

    counts = events.pivot_table(index="Patient", columns="event", aggfunc="size", fill_value=0)
    counts = counts.rename(columns={"Emergence": "n_Em", "Disappearance": "n_Dis"})
    for column in ["n_Em", "n_Dis"]:
        if column not in counts:
            counts[column] = 0
    eligible = transition.filter_to_snp_confirmed_args(source)[["Patient"]].copy()
    eligible["Patient"] = normalize_patient(eligible["Patient"])
    eligible = eligible.dropna().drop_duplicates()
    patients = metadata.merge(eligible, on="Patient", how="inner")
    patients = patients.merge(counts[["n_Em", "n_Dis"]], on="Patient", how="left")
    patients[["n_Em", "n_Dis"]] = patients[["n_Em", "n_Dis"]].fillna(0)

    source = source.copy()
    source["Patient"] = normalize_patient(source["Patient"])
    detected = infer_k_post(source)
    detected["Patient"] = normalize_patient(detected["Patient"])
    patients = patients.merge(detected, on="Patient", how="left")
    patients["k_post_detected"] = patients["k_post_detected"].fillna(0).astype(int)

    patients["k_post"] = pd.to_numeric(patients["Post-fmt_samples"], errors="coerce")
    if patients["k_post"].isna().any():
        missing = patients.loc[patients["k_post"].isna(), "Patient"].tolist()
        raise ValueError(f"Missing or nonnumeric Post-fmt_samples for patients: {missing}")
    if (patients["k_post"] <= 0).any():
        invalid = patients.loc[patients["k_post"] <= 0, "Patient"].tolist()
        raise ValueError(f"Post-fmt_samples must be positive for patients: {invalid}")
    if not np.allclose(patients["k_post"], np.round(patients["k_post"])):
        raise ValueError("Post-fmt_samples must contain whole-number sample counts.")
    patients["k_post_source"] = "meta_data:count_unique_PostFMT_samples"
    patients = patients.merge(build_contig_summaries(sample_metadata), on="Patient", how="left")
    contig_columns = [
        "median_recipient_num_contigs", "median_recipient_total_length",
        "median_recipient_non_host_length",
    ]
    if patients[contig_columns].isna().any().any():
        missing = patients.loc[patients[contig_columns].isna().any(axis=1), "Patient"].tolist()
        raise ValueError(f"No valid recipient assembly metadata for patients: {missing}")

    patients[["n_Em", "n_Dis", "k_post"]] = patients[["n_Em", "n_Dis", "k_post"]].astype(int)
    patients["cohort"] = pd.Categorical(patients["cohort"], COHORTS)
    if patients["cohort"].isna().any():
        raise ValueError("Unexpected cohort label in metadata.")
    return patients.sort_values(["cohort", "study", "Patient"])[[
        "Patient", "cohort", "study", "n_Em", "n_Dis", "k_post",
        "median_recipient_num_contigs", "median_recipient_total_length",
        "median_recipient_non_host_length",
        "k_post_detected", "k_post_source"
    ]]

# Write the four minimal patient-level inputs for Model 2.
def export_model_2(source: pd.DataFrame, metadata: pd.DataFrame,
                   sample_metadata: pd.DataFrame, transition) -> None:
    
    counts = transition.compute_patient_rates_from_definitions(source).rename(columns={
        "Patient": "patient", "Disease_type": "cohort",
        "n_persist": "persistence", "n_disappear": "disappearance",
        "n_emerge": "emergence", "n_transfer": "transfer",
    })
    meta = metadata.rename(columns={
        "Patient": "patient", "Disease_type": "cohort", "Study_ID": "study",
        "Post-fmt_samples": "k_post",
    })[["patient", "cohort", "study", "k_post"]].copy()
    counts["patient"] = normalize_patient(counts["patient"])
    meta["patient"] = normalize_patient(meta["patient"])
    pat = counts.merge(meta, on=["patient", "cohort"], how="left", validate="one_to_one")
    transition_columns = ["persistence", "disappearance", "emergence", "transfer"]
    pat["n_total"] = pat[transition_columns].sum(axis=1).astype(int)
    pat[transition_columns + ["k_post"]] = pat[transition_columns + ["k_post"]].astype(int)
    if pat[["study", "k_post"]].isna().any().any() or (pat["k_post"] <= 0).any():
        raise ValueError("Model 2 requires study and positive k_post for every patient.")

    donor_rows = sample_metadata.loc[
        sample_metadata["donor_receipient"].astype(str).str.strip().str.lower().eq("donor"),
        ["Patient"],
    ].copy()
    donor_rows["patient"] = normalize_patient(donor_rows["Patient"])
    donor_patients = set(donor_rows["patient"].dropna())

    MODEL_2_OUTPUT.mkdir(parents=True, exist_ok=True)
    for transition_type in transition_columns:
        model_data = pat.copy()
        if transition_type == "transfer":
            model_data = model_data.loc[model_data["patient"].isin(donor_patients)].copy()
        model_data = model_data.loc[model_data["n_total"] > 0].copy()
        model_data["n_type"] = model_data[transition_type]
        model_data[["patient", "cohort", "study", "n_type", "n_total", "k_post"]].to_csv(
            MODEL_2_OUTPUT / f"{transition_type}.csv", index=False
        )

def read_detailed_transition(path: Path) -> pd.DataFrame:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle))
        group_index = header.index("MGE_Group_Name")
        patient_index = header.index("Patient_ID")
        for fields in csv.reader(handle):
            if len(fields) <= max(group_index, patient_index):
                raise ValueError(f"Malformed row in {path}: {fields}")
            rows.append({"MGE_Group_Name": fields[group_index], "patient": fields[patient_index]})
    return pd.DataFrame(rows)


def export_model_3(metadata: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
    emerged = read_detailed_transition(EMERGE_DETAILED)
    disappeared = read_detailed_transition(DISAPPEAR_DETAILED)
    for frame in (emerged, disappeared):
        frame["patient"] = normalize_patient(frame["patient"])
        frame["MGE_class"] = frame["MGE_Group_Name"].astype(str).str.strip()

    em_counts = emerged.groupby(["patient", "MGE_class"]).size().rename("n_Em")
    dis_counts = disappeared.groupby(["patient", "MGE_class"]).size().rename("n_Dis")
    model_data = pd.concat([em_counts, dis_counts], axis=1).fillna(0).astype(int).reset_index()
    meta = metadata.rename(columns={
        "Patient": "patient", "Disease_type": "cohort", "Study_ID": "study",
        "Post-fmt_samples": "k_post",
    })[["patient", "cohort", "study", "k_post"]].copy()
    meta["patient"] = normalize_patient(meta["patient"])
    model_data = model_data.merge(meta, on="patient", how="left", validate="many_to_one")
    assembly = build_contig_summaries(sample_metadata).rename(columns={"Patient": "patient"})
    model_data = model_data.merge(assembly, on="patient", how="left", validate="many_to_one")
    required_metadata = [
        "cohort", "study", "k_post", "median_recipient_num_contigs",
        "median_recipient_total_length", "median_recipient_non_host_length",
    ]
    if model_data[required_metadata].isna().any().any():
        raise ValueError("Model 3 metadata are incomplete.")

    represented = model_data.groupby("MGE_class")["cohort"].nunique()
    retained = represented[represented == model_data["cohort"].nunique()].index
    model_data = model_data.loc[model_data["MGE_class"].isin(retained)].copy()
    MODEL_3_OUTPUT.mkdir(parents=True, exist_ok=True)
    model_data[[
        "patient", "study", "cohort", "MGE_class", "n_Em", "n_Dis", "k_post",
        "median_recipient_num_contigs", "median_recipient_total_length",
        "median_recipient_non_host_length",
    ]].sort_values(["patient", "MGE_class"]).to_csv(
        MODEL_3_OUTPUT / "patient_mge_emergence_disappearance.csv", index=False
    )


# "Export patient-level CSVs for leave-one-study-out analysis."""
def extract_data_for_loso(INPUT: Path, OUTPUT_DIR: Path, OUTPUT: Path) -> None:
    pat = pd.read_csv(INPUT)
    columns = ["Patient", "cohort", "study", "n_Em", "n_Dis", "k_post"]
    missing = set(columns) - set(pat.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    pat = pat.loc[pat["n_Em"] + pat["n_Dis"] > 0, columns].copy()
    if pat[columns].isna().any().any() or (pat["k_post"] <= 0).any():
        raise ValueError("LOSO input contains missing values or nonpositive k_post.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pat.to_csv(OUTPUT, index=False)
    print(f"Wrote {len(pat)} patients across {pat['study'].nunique()} studies to {OUTPUT}")

def main() -> None:
    args = arguments()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(args.data)
    sample_metadata = pd.read_csv(args.sample_metadata)
    metadata = build_patient_metadata(sample_metadata)
    transition = load_transition_module()
    events = build_events(source, transition)
    patients = build_patients(source, metadata, events, sample_metadata, transition)
    events.to_csv(args.output_dir / "emergence_disappearance_events.csv", index=False)
    patients.to_csv(args.output_dir / "emergence_disappearance_patients.csv", index=False)
    # export_clinically_relevant_patients(
    #     source, metadata, sample_metadata, args.clinical_arg_list, transition
    # )
    export_model_2(source, metadata, sample_metadata, transition)
    export_model_3(metadata, sample_metadata)
    extract_data_for_loso(INPUT, OUTPUT_DIR, OUTPUT)


if __name__ == "__main__":
    main()
