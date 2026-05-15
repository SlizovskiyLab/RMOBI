import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np

GRAPH_PATH = '../docs/json/graph2.json'
OUTPUT_DIR = Path("images/heatmaps")

# Keep requested x-axis class order fixed across all heatmaps.
MGE_ORDER = ['Plasmids', 'Likely IS/TE', 'ICE', 'Prophages', 'Virus']

MGE_PATTERNS = {
    'Plasmids': ['plasmid'],
    'Likely IS/TE': ['likely is/te', 'is/te', 'transpos', 'insertion sequence'],
    'ICE': ['ice', 'integrative conjugative'],
    'Prophages': ['prophage', 'proph'],
    'Virus': ['virus', 'viral']
}

RESISTANCE_GROUPS = [
    ('betalactams', 'Betalactams'),
    ('aminoglycosides', 'Aminoglycosides'),
    ('mls', 'MLS'),
    ('tetracyclines', 'Tetracyclines'),
    ('trimethoprim', 'Trimethoprim'),
    ('sulfonamides', 'Sulfonamides')
]


def load_graph_nodes(path):
    with open(path, 'r') as f:
        graph_data = json.load(f)
    return graph_data.get('nodes', [])


def normalize_mge(val):
    if pd.isna(val):
        return 'Unclassified'
    val_lower = str(val).lower()
    for group, patterns in MGE_PATTERNS.items():
        if any(pattern in val_lower for pattern in patterns):
            return group
    return 'Unclassified'


def normalize_disease_order(nodes):
    disease_totals = {}
    for node in nodes:
        counts = node.get('diseaseCounts', {})
        for disease, count in counts.items():
            disease_totals[disease] = disease_totals.get(disease, 0) + float(count)

    preferred_order = ['MDRB', 'Melanoma', 'rCDI']
    present_preferred = [d for d in preferred_order if d in disease_totals]
    remaining = sorted([d for d in disease_totals if d not in preferred_order])
    return present_preferred + remaining, disease_totals


def extract_node_records(nodes):
    records = []
    for node in nodes:
        arg_group = str(node.get('argGroup', ''))
        mge_group = normalize_mge(node.get('mgeGroup'))
        disease_counts = node.get('diseaseCounts', {})

        for disease, count in disease_counts.items():
            records.append(
                {
                    'argGroup': arg_group,
                    'mgeGroup_normalized': mge_group,
                    'disease': disease,
                    'count': float(count),
                }
            )

    return pd.DataFrame(records)


def resistance_match(arg_group_value, target_key):
    arg_lower = str(arg_group_value).lower()

    if target_key == 'betalactams':
        return 'betalactam' in arg_lower
    if target_key == 'aminoglycosides':
        return 'aminoglycosid' in arg_lower
    if target_key == 'mls':
        return 'mls' in arg_lower
    if target_key == 'tetracyclines':
        return ('tetracycline' in arg_lower) or ('tettracycline' in arg_lower)
    if target_key == 'trimethoprim':
        return 'trimethoprim' in arg_lower
    if target_key == 'sulfonamides':
        return 'sulfonamide' in arg_lower

    return False


def truncated_cmap(base_cmap: str, min_val: float = 0.30, max_val: float = 0.80):
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"trunc_{base_cmap}",
        plt.get_cmap(base_cmap)(np.linspace(min_val, max_val, 256)),
    )


BLUE_CMAP = truncated_cmap("Blues", 0.0, 0.99)

def plot_heatmap(df_long, disease_order, resistance_key, resistance_label):
    mask = df_long['argGroup'].apply(lambda x: resistance_match(x, resistance_key))
    df_filtered = df_long[mask].copy()

    if len(df_filtered) == 0:
        heatmap_counts = pd.DataFrame(0.0, index=disease_order, columns=MGE_ORDER)
    else:
        heatmap_counts = (
            df_filtered
            .groupby(['disease', 'mgeGroup_normalized'])['count']
            .sum()
            .unstack(fill_value=0)
            .reindex(index=disease_order, columns=MGE_ORDER, fill_value=0)
        )

    # Remove MGE classes that are zero across all diseases for this resistance group.
    nonzero_columns = heatmap_counts.columns[heatmap_counts.sum(axis=0) > 0]
    heatmap_counts = heatmap_counts.loc[:, nonzero_columns]

    if heatmap_counts.shape[1] == 0:
        plt.figure(figsize=(10, 6))
        plt.text(
            0.5,
            0.5,
            f'No non-zero MGE classes for {resistance_label}',
            ha='center',
            va='center',
            fontsize=18,
            fontweight='bold'
        )
        plt.axis('off')

        out_prefix = OUTPUT_DIR / f"heatmap_graph2_disease_mge_{resistance_key}"
        filename = out_prefix.with_suffix(".png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.savefig(out_prefix.with_suffix(".svg"), bbox_inches='tight')
        plt.close()

        print(f"Heatmap saved as '{filename}'")
        print(f"\nNo non-zero MGE classes to plot for {resistance_label}.")
        print("\n" + "="*80 + "\n")
        return

    row_totals = heatmap_counts.sum(axis=1).replace(0, np.nan)
    heatmap_pct = (heatmap_counts.div(row_totals, axis=0) * 100).fillna(0)

    plt.figure(figsize=(14, 8))
    ax = sns.heatmap(
        heatmap_pct,
        annot=True,
        fmt='.1f',
        cmap=BLUE_CMAP,
        annot_kws={'size': 28, 'weight': 'bold'},
        cbar_kws={'label': 'Percentage (%)', 'shrink': 0.8},
        linewidths=0.1,
        linecolor='gray'
    )
    
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=26, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=28, fontweight='bold')
    
    cbar = ax.collections[0].colorbar
    cbar.set_label('Percentage (%)', fontsize=28, fontweight='bold')
    cbar.ax.tick_params(labelsize=28)

    # plt.title(f'{resistance_label}: Disease vs MGE Class', fontsize=18, fontweight='bold')
    plt.tight_layout()

    out_prefix = OUTPUT_DIR / f"heatmap_graph2_disease_mge_{resistance_key}"
    filename = out_prefix.with_suffix(".png")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.savefig(out_prefix.with_suffix(".svg"), bbox_inches='tight')
    plt.close()

    print(f"Heatmap saved as '{filename}'")
    # print(f"\nHeatmap Data (Disease vs MGE for {resistance_label} - % within each disease):")
    # print(heatmap_pct)
    # print(f"\nMatched nodes for {resistance_label}: {len(df_filtered)}")
    # print(f"Total disease-count assignments ({resistance_label}): {df_filtered['count'].sum():.0f}")
    # print("\n" + "="*80 + "\n")

nodes = load_graph_nodes(GRAPH_PATH)
df_long = extract_node_records(nodes)
disease_order, disease_totals = normalize_disease_order(nodes)

print("Diseases inferred from nodes[].diseaseCounts with aggregate counts:")
for disease in disease_order:
    print(f"- {disease}: {int(disease_totals[disease])}")
print("\n" + "="*80 + "\n")

for resistance_key, resistance_label in RESISTANCE_GROUPS:
    plot_heatmap(df_long, disease_order, resistance_key, resistance_label)
