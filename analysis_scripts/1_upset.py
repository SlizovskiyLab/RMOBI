import re
from pathlib import Path
from typing import Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ANALYSIS_DIR = Path(__file__).resolve().parent
DATA_PATH = ANALYSIS_DIR.parent / "data" / "patientwise_colocalization_by_timepoint.csv"
OUTPUT_DIR = Path("images/upset")
DISEASES = ["MDRB", "Melanoma", "rCDI"]

PHASES = ["Donor", "PreFMT", "PostFMT1", "PostFMT2", "PostFMT3"]
PHASE_LABELS = {"Donor": "Donor", "PreFMT": "Pre-FMT", "PostFMT1": "Post-FMT-30", "PostFMT2": "Post-FMT-60", "PostFMT3": "Post-FMT-60+"}

def postfmt_bin(colname: str):
    m = re.fullmatch(r"PostFMT_(\d{3})", colname)
    if not m:
        return None
    day = int(m.group(1))
    if 1 <= day <= 30:
        return "PostFMT1"
    if 31 <= day <= 60:
        return "PostFMT2"
    return "PostFMT3"

def build_coloc_pattern_counts_and_set_sizes(df: pd.DataFrame, disease="MDRB",
                                            coloc_id_cols=("Patient","MEGARes group","MGE gene")):
    df = df.copy()
    df["Disease_type"] = df["Disease_type"].astype(str).str.strip()
    df = df[df["Disease_type"].str.lower() == disease.lower()].copy()

    # group PostFMT columns
    post_cols = [c for c in df.columns if postfmt_bin(c) is not None]
    groups = {"PostFMT1": [], "PostFMT2": [], "PostFMT3": []}
    for c in post_cols:
        groups[postfmt_bin(c)].append(c)

    # row-level phase presence (0/1)
    df["Donor_p"]  = (pd.to_numeric(df["Donor"],  errors="coerce").fillna(0) > 0).astype(int)
    df["PreFMT_p"] = (pd.to_numeric(df["PreFMT"], errors="coerce").fillna(0) > 0).astype(int)

    for g in ["PostFMT1", "PostFMT2", "PostFMT3"]:
        if groups[g]:
            block = df[groups[g]].apply(pd.to_numeric, errors="coerce").fillna(0)
            df[f"{g}_p"] = (block.max(axis=1) > 0).astype(int)
        else:
            df[f"{g}_p"] = 0

    # pattern per row (exact combination)
    def row_pattern(r):
        present = []
        if r["Donor_p"]:    present.append("Donor")
        if r["PreFMT_p"]:   present.append("PreFMT")
        if r["PostFMT1_p"]: present.append("PostFMT1")
        if r["PostFMT2_p"]: present.append("PostFMT2")
        if r["PostFMT3_p"]: present.append("PostFMT3")
        return "+".join(present) if present else None

    df["pattern"] = df.apply(row_pattern, axis=1)
    df = df[df["pattern"].notna()].copy()

    # Ensure "colocalization" counted once (in case of duplicates)
    df_unique = df.drop_duplicates(list(coloc_id_cols) + ["pattern"])

    # Intersection size = # colocalizations per pattern
    counts = df_unique.groupby("pattern").size().sort_values(ascending=False)

    # Left set sizes = # colocalizations that appear in each phase (any time within that phase bin)
    set_sizes = {
        "Donor":    int(df_unique[df_unique["Donor_p"]    == 1].shape[0]),
        "PreFMT":   int(df_unique[df_unique["PreFMT_p"]   == 1].shape[0]),
        "PostFMT1": int(df_unique[df_unique["PostFMT1_p"] == 1].shape[0]),
        "PostFMT2": int(df_unique[df_unique["PostFMT2_p"] == 1].shape[0]),
        "PostFMT3": int(df_unique[df_unique["PostFMT3_p"] == 1].shape[0]),
    }

    return counts, set_sizes


def plot_upset_like_clear(counts: pd.Series,
                          set_sizes: dict,
                          top_n=31,
                          figsize=(26, 8.0),
                          dpi=400,
                          dot_grey=80, dot_black=100, line_w=0.6,
                          bar_label_font=20,
                          label_only_if_ge=1,   # don’t label tiny bars (reduces clutter)
                          out_prefix="MDRB_upset_coloc"):
    out_prefix = Path(out_prefix)
    counts = counts.sort_values(ascending=False).head(top_n)
    patterns = counts.index.tolist()
    values = counts.values.astype(int)

    # membership matrix (rows=phases, cols=patterns)
    M = np.zeros((len(PHASES), len(patterns)), dtype=bool)
    for j, pat in enumerate(patterns):
        s = set(pat.split("+"))
        for i, ph in enumerate(PHASES):
            M[i, j] = (ph in s)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = GridSpec(
        nrows=2, ncols=2,
        height_ratios=[3.0, 1.55],
        width_ratios=[1.8, 6.0],
        hspace=0.15, wspace=0.34,
        left=0.16, right=0.95, top=0.95, bottom=0.12
    )

    ax_top  = fig.add_subplot(gs[0, 1])
    ax_mat  = fig.add_subplot(gs[1, 1], sharex=ax_top)
    ax_left = fig.add_subplot(gs[1, 0], sharey=ax_mat)

    x = np.arange(len(patterns))

    # --- Top bars: intersection size (# colocalizations) ---
    colors = plt.cm.YlGnBu(np.linspace(0.45, 0.95, len(patterns)))
    bars = ax_top.bar(x, values, width=0.58, color=colors, edgecolor="black", linewidth=0.9)
    ax_top.set_ylabel("No. of colocalizations", fontsize=24, fontweight="bold", labelpad=10)
    ax_top.set_xticks([])

    ymax = values.max() if len(values) else 0
    ax_top.set_ylim(0, ymax * 1.12 + 1)

    # Bar labels: smaller + only label if big enough
    for b, v in zip(bars, values):
        if v >= label_only_if_ge:
            ax_top.text(
                b.get_x() + b.get_width()/2,
                v + max(1, int(0.01 * ymax)),
                str(v),
                ha="center", va="bottom",
                fontsize=bar_label_font
            )

    ax_top.tick_params(axis="y", labelsize=22)
    for label in ax_top.get_yticklabels():
        label.set_fontweight("bold")

    # --- Matrix ---
    ax_mat.set_yticks(np.arange(len(PHASES)))
    ax_mat.set_yticklabels([PHASE_LABELS[p] for p in PHASES], fontsize=24, fontweight="bold",)
    ax_mat.tick_params(axis="y", bottom=False, labelbottom=False)
    ax_mat.set_xlim(-0.5, len(patterns)-0.5)

    # grey background dots
    for i in range(len(PHASES)):
        ax_mat.scatter(x, np.full_like(x, i), s=dot_grey, color="lightgrey", zorder=1)

    # black dots + connecting lines
    for j in range(len(patterns)):
        rows = np.where(M[:, j])[0]
        ax_mat.scatter([j]*len(rows), rows, s=dot_black, color="black", zorder=3)
        if len(rows) > 1:
            ax_mat.plot([j, j], [rows.min(), rows.max()], color="black", linewidth=line_w, zorder=2)

    for spine in ["top", "right"]:
        ax_mat.spines[spine].set_visible(False)

    # --- Left set sizes (# colocalizations per phase) ---
    phase_sizes = np.array([set_sizes.get(p, 0) for p in PHASES], dtype=int)
    left_colors = plt.cm.YlGnBu(np.linspace(0.45, 0.95, len(PHASES)))
    ax_left.barh(
        np.arange(len(PHASES)),
        phase_sizes,
        color=left_colors,
        edgecolor="black",
        linewidth=0.9,
    )
    ax_left.set_xlabel("No. of colocalizations", fontsize=24, fontweight="bold", labelpad=4)
    ax_left.tick_params(axis="x", labelsize=22)
    for label in ax_left.get_xticklabels():
        label.set_fontweight("bold")
    ax_left.tick_params(axis="y", left=False, labelleft=False)
    # ax_left.set_yticks(np.arange(len(PHASES)))
    # ax_left.set_yticklabels([PHASE_LABELS[p] for p in PHASES], fontsize=10)
    # ax_left.tick_params(axis="y", left=True, labelleft=True)

    for spine in ["top", "right"]:
        ax_left.spines[spine].set_visible(False)

    # ax_top.set_title("Longitudinal overlap of ARG–MGE colocalizations across FMT phases in MDRB", fontsize=18, fontweight="bold", pad=10)

    fig.savefig(out_prefix.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), dpi=dpi, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    # plt.show()

"""

"""

def _extract_initializer_block(cpp_text: str, map_name: str) -> str:
    """
    Pull out the {...} initializer content for a given map variable name.
    Works for your simple initializer-list style.
    """
    # match: <map_name> = { ... };
    m = re.search(rf"\b{re.escape(map_name)}\b\s*=\s*\{{(.*?)\}}\s*;", cpp_text, flags=re.S)
    if not m:
        raise ValueError(f"Could not find initializer for '{map_name}'")
    return m.group(1)

def parse_cpp_int_to_string_map(cpp_text: str, map_name: str) -> Dict[int, str]:
    block = _extract_initializer_block(cpp_text, map_name)
    entries = re.findall(r"\{\s*(\d+)\s*,\s*\"([^\"]+)\"\s*\}", block)
    return {int(k): v for k, v in entries}

def parse_cpp_int_to_bool_map(cpp_text: str, map_name: str) -> Dict[int, bool]:
    block = _extract_initializer_block(cpp_text, map_name)
    entries = re.findall(r"\{\s*(\d+)\s*,\s*(true|false)\s*\}", block, flags=re.I)
    return {int(k): (v.lower() == "true") for k, v in entries}

def print_unique_arg_snp_summary(
    df: pd.DataFrame,
    id_to_group: Dict[int, str],
    id_to_requires_snp_confirmation: Dict[int, bool],
) -> None:
    unique_args = set(df["MEGARes group"].astype(str).str.strip().dropna())
    group_to_snp_status: Dict[str, set] = {}

    for arg_id, arg_group in id_to_group.items():
        if arg_id not in id_to_requires_snp_confirmation:
            continue
        group_to_snp_status.setdefault(arg_group, set()).add(
            id_to_requires_snp_confirmation[arg_id]
        )

    snp_confirmed_args = {
        arg for arg in unique_args
        if arg in group_to_snp_status and False in group_to_snp_status[arg]
    }
    snp_not_confirmed_args = {
        arg for arg in unique_args
        if arg in group_to_snp_status and True in group_to_snp_status[arg]
    }
    unmapped_args = unique_args - set(group_to_snp_status)

    print("\nUnique ARG SNP-confirmation summary")
    print(f"Unique ARG groups in input: {len(unique_args)}")
    print(f"SNP confirmed / usable ARG groups: {len(snp_confirmed_args)}")
    print(f"SNP not confirmed / requires SNP confirmation ARG groups: {len(snp_not_confirmed_args)}")
    print(f"ARG groups missing from SNP map: {len(unmapped_args)}")

    print("\nUnique ARG groups by disease")
    for disease in DISEASES:
        disease_args = set(
            df.loc[
                df["Disease_type"].astype(str).str.strip().str.lower() == disease.lower(),
                "MEGARes group",
            ]
            .astype(str)
            .str.strip()
            .dropna()
        )
        disease_snp_confirmed = disease_args & snp_confirmed_args
        disease_snp_not_confirmed = disease_args & snp_not_confirmed_args
        disease_unmapped = disease_args & unmapped_args
        print(
            f"{disease}: total={len(disease_args)}, "
            f"snp_confirmed={len(disease_snp_confirmed)}, "
            f"snp_not_confirmed={len(disease_snp_not_confirmed)}, "
            f"unmapped={len(disease_unmapped)}"
        )

# ---- USAGE ----
# C++ code (both maps) into a string:
cpp_text = r"""

const std::unordered_map<int, std::string> argIdMap = {
    {1, "A16S"},
    {2, "AAC3"},
    {3, "AAC6-PRIME"},
    {4, "ACN"},
    {5, "ACRA"},
    {6, "ACRB"},
    {7, "ACRD"},
    {8, "ACRR"},
    {9, "ACRS"},
    {10, "AMPCR"},
    {11, "AMPH"},
    {12, "ANT3-DPRIME"},
    {13, "ANT6"},
    {14, "APH3-DPRIME"},
    {15, "APH3-PRIME"},
    {16, "APH6"},
    {17, "ARIR"},
    {18, "ARNT"},
    {19, "ARR"},
    {20, "ARSA"},
    {21, "ARSB"},
    {22, "ARSBM"},
    {23, "ARSCM"},
    {24, "ARSDM"},
    {25, "ARSP"},
    {26, "ARSRM"},
    {27, "ASMA"},
    {28, "ASR"},
    {29, "BACA"},
    {30, "BAES"},
    {31, "BCR"},
    {32, "BHSA"},
    {33, "BLA1"},
    {34, "BLAEC"},
    {35, "BLAZ"},
    {36, "CADX"},
    {37, "CAP16S"},
    {38, "CATA"},
    {39, "CDEA"},
    {40, "CEPA"},
    {41, "CEPAB"},
    {42, "CFX"},
    {43, "CHAA"},
    {44, "CLS"},
    {45, "CMY"},
    {46, "COMR"},
    {47, "COPA"},
    {48, "CORA"},
    {49, "CORB"},
    {50, "CORC"},
    {51, "CPXAR"},
    {52, "CRP"},
    {53, "CTX"},
    {54, "CUEO"},
    {55, "CUER"},
    {56, "CUID"},
    {57, "CUSA"},
    {58, "CUSB"},
    {59, "CUSC"},
    {60, "CUSR"},
    {61, "CUSS"},
    {62, "CUTA"},
    {63, "CUTC"},
    {64, "DFRA"},
    {65, "DFRC"},
    {66, "DFRE"},
    {67, "DFRF"},
    {68, "DHA"},
    {69, "DNAK"},
    {70, "DSBA"},
    {71, "DSBC"},
    {72, "EATAV"},
    {73, "EFMA"},
    {74, "EMEA"},
    {75, "EMMDR"},
    {76, "EMRB"},
    {77, "EMRD"},
    {78, "EMRK"},
    {79, "EMRR"},
    {80, "EMRY"},
    {81, "EPTA"},
    {82, "EPTB"},
    {83, "ERMB"},
    {84, "ERMC"},
    {85, "ERMF"},
    {86, "ERMX"},
    {87, "EVGA"},
    {88, "FABG"},
    {89, "FABI"},
    {90, "FECD"},
    {91, "FETA"},
    {92, "FETB"},
    {93, "FIEF"},
    {94, "FLOR"},
    {95, "FOLP"},
    {96, "FOSA"},
    {97, "FUSA"},
    {98, "GADA"},
    {99, "GADC"},
    {100, "GADX"},
    {101, "GLPF"},
    {102, "GLPT"},
    {103, "GYRA"},
    {104, "GYRB"},
    {105, "GYRBA"},
    {106, "HDEA"},
    {107, "HNS"},
    {108, "ICLR"},
    {109, "ILES"},
    {110, "KDEA"},
    {111, "KDPE"},
    {112, "KEXD"},
    {113, "KMRA"},
    {114, "KPC"},
    {115, "KPN"},
    {116, "KPNE"},
    {117, "KPNF"},
    {118, "KPNO"},
    {119, "LEN"},
    {120, "LMRD"},
    {121, "LNUA"},
    {122, "LNUC"},
    {123, "LNUG"},
    {124, "LPDT"},
    {125, "LPTD"},
    {126, "LSA"},
    {127, "MARA"},
    {128, "MARR"},
    {129, "MCR"},
    {130, "MDE"},
    {131, "MDFA"},
    {132, "MDTA"},
    {133, "MDTB"},
    {134, "MDTJ"},
    {135, "MDTK"},
    {136, "MDTM"},
    {137, "MECA"},
    {138, "MEFA"},
    {139, "MEFE"},
    {140, "MERA"},
    {141, "MERD"},
    {142, "MERR"},
    {143, "MERR1"},
    {144, "MERR2"},
    {145, "MERT"},
    {146, "MGTA"},
    {147, "MLS23S"},
    {148, "MNTH"},
    {149, "MNTP"},
    {150, "MNTR"},
    {151, "MODA"},
    {152, "MODC"},
    {153, "MPHA"},
    {154, "MPRF"},
    {155, "MSBA"},
    {156, "MSRD"},
    {157, "MURA"},
    {158, "MVRC"},
    {159, "NCRA"},
    {160, "NCRAN"},
    {161, "NFSA"},
    {162, "NHAA"},
    {163, "NHAB"},
    {164, "NIKC"},
    {165, "NIRA"},
    {166, "NMPC"},
    {167, "NORA"},
    {168, "O23S"},
    {169, "OMP36"},
    {170, "OMP37"},
    {171, "OMPA"},
    {172, "OMPD"},
    {173, "OMPF"},
    {174, "OMPFB"},
    {175, "OMPK36"},
    {176, "OQXA"},
    {177, "OQXB"},
    {178, "ORN"},
    {179, "OXA"},
    {180, "OXY"},
    {181, "OXYRKP"},
    {182, "PARC"},
    {183, "PARE"},
    {184, "PAREF"},
    {185, "PATB"},
    {186, "PBP4B"},
    {187, "PCOA"},
    {188, "PCOE"},
    {189, "PHOB"},
    {190, "PHOR"},
    {191, "PITA"},
    {192, "PMRF"},
    {193, "PSTC"},
    {194, "PSTS"},
    {195, "PTSL"},
    {196, "QACEDELTA1"},
    {197, "QACG"},
    {198, "QNRB"},
    {199, "QNRD"},
    {200, "QNRS"},
    {201, "RAMR"},
    {202, "RCNR"},
    {203, "RMTC"},
    {204, "ROBA"},
    {205, "RPOB"},
    {206, "RPOS"},
    {207, "RPSA"},
    {208, "RPSL"},
    {209, "RRS"},
    {210, "RRSA"},
    {211, "RRSC"},
    {212, "RRSH"},
    {213, "SAT"},
    {214, "SDEY"},
    {215, "SHV"},
    {216, "SILA"},
    {217, "SILB"},
    {218, "SILC"},
    {219, "SILE"},
    {220, "SILF"},
    {221, "SILP"},
    {222, "SILS"},
    {223, "SITABCD"},
    {224, "SMDA"},
    {225, "SMVA"},
    {226, "SODA"},
    {227, "SODB"},
    {228, "SOXRB"},
    {229, "SUGE"},
    {230, "SULI"},
    {231, "SULII"},
    {232, "SULIII"},
    {233, "TEHA"},
    {234, "TEM"},
    {235, "TERB"},
    {236, "TERD"},
    {237, "TERW"},
    {238, "TERZ"},
    {239, "TET16S"},
    {240, "TET40"},
    {241, "TETA"},
    {242, "TETA46"},
    {243, "TETB"},
    {244, "TETD"},
    {245, "TETM"},
    {246, "TETO"},
    {247, "TETQ"},
    {248, "TETR"},
    {249, "TETW"},
    {250, "TETX"},
    {251, "TOLC"},
    {252, "TUFAB"},
    {253, "UGD"},
    {254, "UHPT"},
    {255, "VANC"},
    {256, "VANHA"},
    {257, "VANRA"},
    {258, "YBTP"},
    {259, "YCHH"},
    {260, "YDEI"},
    {261, "YDEO"},
    {262, "YDEP"},
    {263, "YGJH"},
    {264, "YJAA"},
    {265, "YJCG"},
    {266, "YODB"},
    {267, "YOGI"},
    {268, "ZINT"},
    {269, "ZNTA"},
    {270, "ZNTR"},
    {271, "ZNUA"},
    {272, "ZNUC"},
    {273, "ZRAS"},
    {274, "ZUPT"},
    {275, "ZUR"}
};

const std::unordered_map<int, bool> argIDSNPConfirmation = {
    {1, true},
    {2, false},
    {3, false},
    {4, false},
    {5, false},
    {6, true},
    {7, false},
    {8, true},
    {9, false},
    {10, true},
    {11, false},
    {12, false},
    {13, false},
    {14, false},
    {15, false},
    {16, false},
    {17, false},
    {18, false},
    {19, false},
    {20, false},
    {21, false},
    {22, false},
    {23, false},
    {24, false},
    {25, false},
    {26, false},
    {27, false},
    {28, false},
    {29, false},
    {30, false},
    {31, false},
    {32, false},
    {33, false},
    {34, false},
    {35, false},
    {36, false},
    {37, true},
    {38, false},
    {39, false},
    {40, false},
    {41, false},
    {42, false},
    {43, false},
    {44, true},
    {45, false},
    {46, false},
    {47, false},
    {48, false},
    {49, false},
    {50, false},
    {51, false},
    {52, false},
    {53, false},
    {54, false},
    {55, false},
    {56, false},
    {57, false},
    {58, false},
    {59, false},
    {60, false},
    {61, false},
    {62, false},
    {63, false},
    {64, false},
    {65, true},
    {66, false},
    {67, false},
    {68, false},
    {69, false},
    {70, false},
    {71, false},
    {72, true},
    {73, false},
    {74, false},
    {75, false},
    {76, false},
    {77, false},
    {78, false},
    {79, false},
    {80, false},
    {81, false},
    {82, false},
    {83, false},
    {84, false},
    {85, false},
    {86, false},
    {87, false},
    {88, true},
    {89, false},
    {90, false},
    {91, false},
    {92, false},
    {93, false},
    {94, false},
    {95, true},
    {96, false},
    {97, true},
    {98, false},
    {99, false},
    {100, false},
    {101, false},
    {102, true},
    {103, true},
    {104, true},
    {105, true},
    {106, false},
    {107, false},
    {108, false},
    {109, true},
    {110, false},
    {111, false},
    {112, false},
    {113, false},
    {114, false},
    {115, false},
    {116, false},
    {117, false},
    {118, false},
    {119, false},
    {120, false},
    {121, false},
    {122, false},
    {123, false},
    {124, false},
    {125, false},
    {126, false},
    {127, false},
    {128, false},
    {129, false},
    {130, false},
    {131, false},
    {132, false},
    {133, false},
    {134, false},
    {135, false},
    {136, false},
    {137, false},
    {138, false},
    {139, false},
    {140, false},
    {141, false},
    {142, false},
    {143, false},
    {144, false},
    {145, false},
    {146, false},
    {147, true},
    {148, false},
    {149, false},
    {150, false},
    {151, false},
    {152, false},
    {153, false},
    {154, false},
    {155, false},
    {156, false},
    {157, true},
    {158, false},
    {159, false},
    {160, false},
    {161, false},
    {162, false},
    {163, false},
    {164, false},
    {165, false},
    {166, false},
    {167, false},
    {168, true},
    {169, true},
    {170, false},
    {171, false},
    {172, false},
    {173, true},
    {174, true},
    {175, true},
    {176, false},
    {177, false},
    {178, false},
    {179, false},
    {180, false},
    {181, false},
    {182, true},
    {183, true},
    {184, true},
    {185, false},
    {186, false},
    {187, false},
    {188, false},
    {189, true},
    {190, false},
    {191, false},
    {192, false},
    {193, false},
    {194, false},
    {195, true},
    {196, false},
    {197, false},
    {198, false},
    {199, false},
    {200, false},
    {201, true},
    {202, false},
    {203, false},
    {204, false},
    {205, true},
    {206, false},
    {207, true},
    {208, true},
    {209, true},
    {210, true},
    {211, true},
    {212, true},
    {213, false},
    {214, false},
    {215, false},
    {216, false},
    {217, false},
    {218, false},
    {219, false},
    {220, false},
    {221, false},
    {222, false},
    {223, false},
    {224, false},
    {225, false},
    {226, false},
    {227, false},
    {228, false},
    {229, false},
    {230, false},
    {231, false},
    {232, false},
    {233, false},
    {234, false},
    {235, false},
    {236, false},
    {237, false},
    {238, false},
    {239, true},
    {240, false},
    {241, false},
    {242, false},
    {243, false},
    {244, false},
    {245, false},
    {246, false},
    {247, false},
    {248, true},
    {249, false},
    {250, false},
    {251, false},
    {252, true},
    {253, false},
    {254, true},
    {255, false},
    {256, false},
    {257, false},
    {258, false},
    {259, false},
    {260, false},
    {261, false},
    {262, false},
    {263, false},
    {264, false},
    {265, false},
    {266, false},
    {267, false},
    {268, false},
    {269, false},
    {270, false},
    {271, false},
    {272, false},
    {273, false},
    {274, false},
    {275, false}
};)
"""

id_to_group = parse_cpp_int_to_string_map(cpp_text, "argIdMap")
id_to_requires_snp_confirmation = parse_cpp_int_to_bool_map(cpp_text, "argIDSNPConfirmation")

snp_confirmed_groups = {
    id_to_group[i]
    for i, requires_snp_confirmation in id_to_requires_snp_confirmation.items()
    if not requires_snp_confirmation and i in id_to_group
}

print("SNP-confirmed groups:", len(snp_confirmed_groups))
# print(sorted(list(snp_confirmed_groups))[:], "...")


if not DATA_PATH.exists():
    raise FileNotFoundError(f"Could not find CSV file: {DATA_PATH}")

df = pd.read_csv(DATA_PATH)
df = df.copy()
df["MEGARes group"] = df["MEGARes group"].astype(str).str.strip()
print_unique_arg_snp_summary(df, id_to_group, id_to_requires_snp_confirmation)

df_snp = df[df["MEGARes group"].isin(snp_confirmed_groups)].copy()

print("Rows before:", len(df), "Rows after SNP filter:", len(df_snp))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Generate UpSet plots for all three disease types
for disease in DISEASES:
    counts_coloc, set_sizes_coloc = build_coloc_pattern_counts_and_set_sizes(df_snp, disease=disease)
    if counts_coloc.empty:
        print(f"No data for {disease}. Skipping.")
        continue
    plot_upset_like_clear(
        counts_coloc,
        set_sizes_coloc,
        top_n=31,
        out_prefix=OUTPUT_DIR / f"{disease.lower()}_upset_coloc",
    )
    print(f"Saved UpSet plot for {disease}")
