import json
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Load the data
with open('../docs/json/graph2.json', 'r') as f:
    data = json.load(f)

# Define disease types
disease_types = ['MDRB', 'Melanoma', 'rCDI']

cohort_sizes = {
    'MDRB': 14,
    'Melanoma': 12,
    'rCDI': 34,
}

# Map graph2.json argGroup values to display categories
ARGGROUP_MAP = {
    'betalactams': 'Betalactams',
    'Aminoglycosides': 'Aminoglycosides',
    'MLS': 'MLS',
    'Tetracyclines': 'Tetracyclines',
    'Trimethoprim': 'Trimethoprim',
    'Sulfonamides': 'Sulfonamides',
}

display_order = [
    'Betalactams',
    'Aminoglycosides',
    'MLS',
    'Tetracyclines',
    'Trimethoprim',
    'Sulfonamides'
]

earthy_colors = {
    'Betalactams': '#1b5a9d',      # brick red
    'Aminoglycosides': '#3c7eb7',  # burnt sienna
    'MLS': '#92beda',              # light coral
    'Tetracyclines': '#9a2d29',    # brick red
    'Trimethoprim': '#c33e36',     # deep sky blue
    'Sulfonamides': '#e25e49'      # light turquoise
}


def autopct_threshold(threshold=3.0):
    """Show percentage text only for wedges above a threshold."""
    def _autopct(pct):
        return f'{pct:.1f}%' if pct >= threshold else ''
    return _autopct

def plot_pie_chart(disease):
    """
    Create a pie chart showing colocalization distribution across resistance groups.
    
    Parameters:
    disease (str): Disease type (MDRB, Melanoma, rCDI)
    """
    # Get data for this disease from graph2.json nodes
    nodes = data['nodes']
    cohort_size = cohort_sizes[disease]

    total_patients = sum(n['diseaseCounts'].get(disease, 0) for n in nodes)
    total_colocalizations = sum(1 for n in nodes if disease in n['diseaseCounts'])

    if total_colocalizations == 0:
        print(f"No nodes found for {disease}; skipping.")
        return

    # Aggregate patient counts per target resistance category.
    counts = {}
    for n in nodes:
        count = n['diseaseCounts'].get(disease, 0)
        if count == 0:
            continue
        category = ARGGROUP_MAP.get(n['argGroup'])
        if category:
            counts[category] = counts.get(category, 0) + count

    resistance_counts = pd.Series(counts).reindex(display_order).dropna()
    plotted_colocalizations = sum(
        1 for n in nodes
        if disease in n['diseaseCounts'] and ARGGROUP_MAP.get(n['argGroup'])
    )
    plotted_total = int(resistance_counts.sum())
    plotted_percentage = (plotted_total / total_patients) * 100 if total_patients else 0

    if resistance_counts.empty or plotted_total == 0:
        print(f"No selected resistance-group nodes found for {disease}; skipping.")
        return

    colors = [earthy_colors[group] for group in resistance_counts.index]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(13, 9))
    
    # Create donut chart
    wedges, texts, autotexts = ax.pie(
        resistance_counts.values,
        labels=None,
        autopct=autopct_threshold(3.0),
        startangle=90,
        pctdistance=0.78,
        colors=colors,
        textprops={'fontsize': 20, 'weight': 'bold'},
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.3, 'width': 0.42}
    )

    center_text = (
        f"{disease}\n"
        f"{plotted_percentage:.1f}% selected\n"
        f"n={plotted_colocalizations} colocalizations\n"
        f"{cohort_size} patients"
    )
    ax.text(
        0,
        0,
        center_text,
        ha='center',
        va='center',
        fontsize=20,
        fontweight='bold'
    )

    ax.legend(
        wedges,
        resistance_counts.index,
        title='Resistance Group',
        loc='center left',
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
        fontsize=20,
        title_fontsize=22
    )
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(20)
        autotext.set_weight('bold')
    
    # Set title
    # ax.set_title(f'{disease} - Colocalization Distribution by Resistance Group', 
    #              fontsize=28, weight='bold', pad=24)
    
    plt.tight_layout()
    
    Path('images').mkdir(exist_ok=True)
    
    # Save as PNG
    png_path = f'images/piecharts/pie_resistance_groups_{disease.lower()}.png'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {png_path}")
    
    # Save as SVG
    svg_path = f'images/piecharts/pie_resistance_groups_{disease.lower()}.svg'
    plt.savefig(svg_path, format='svg', bbox_inches='tight')
    print(f"Saved: {svg_path}")
    print(
        f"{disease}: plotted {plotted_colocalizations}/{total_colocalizations} colocalizations "
        f"from selected resistance groups; "
        f"patient counts: {plotted_total}/{int(total_patients)} ({plotted_percentage:.1f}%); "
        f"cohort size = {cohort_size} patients"
    )
    
    plt.close()

# Generate pie charts for all diseases
for disease in disease_types:
    plot_pie_chart(disease)

print("\nDone! Generated pie charts for all disease types.")
