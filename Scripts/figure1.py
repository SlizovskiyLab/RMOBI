import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_styled_box(ax, x, y, width, height, text, color='#E3F2FD', shape='rect'):
    center_x = x + width/2
    center_y = y + height/2
    
    if shape == 'rect':
        p = patches.FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.2", 
                                   fc=color, ec='#455A64', linewidth=2)
    elif shape == 'cylinder':
        p = patches.Rectangle((x, y), width, height, fc=color, ec='#455A64', linewidth=2)
        ax.plot([x, x+width], [y+height-0.5, y+height-0.5], color='#455A64', lw=1)
    elif shape == 'circle':
        p = patches.Circle((center_x, center_y), width/2, fc=color, ec='#455A64', linewidth=2)
    
    ax.add_patch(p)
    ax.text(center_x, center_y, text, ha='center', va='center', 
            fontsize=9, fontweight='bold', color='#263238')
    return center_x, center_y

def draw_arrow(ax, x1, y1, x2, y2, style='->'):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), 
                arrowprops=dict(arrowstyle=style, lw=1.5, color='#37474F'))

def create_figure():
    fig, ax = plt.subplots(figsize=(12, 14))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis('off')
    
    ax.add_patch(patches.Rectangle((2, 68), 96, 30, fc='#F5F5F5', ec='#B0BEC5', linestyle='--'))
    ax.text(5, 94, "Panel A: Data Processing & Mobilome Discovery (TELCoMB)", fontsize=12, fontweight='bold')

    draw_styled_box(ax, 5, 75, 12, 10, "263 Fecal\nSamples", color='#FFF3E0', shape='cylinder')
    draw_arrow(ax, 18, 80, 23, 80)

    draw_styled_box(ax, 23, 75, 12, 10, "Host\nRemoval", color='white')

    draw_arrow(ax, 35, 80, 40, 86)
    
    draw_styled_box(ax, 40, 81, 14, 10, "MetaSPAdes\n(Assembly)", color='white')
    
    draw_styled_box(ax, 40, 69, 14, 10, "Databases\n(ARG/MGE)", color='#FFF3E0', shape='cylinder')

    draw_arrow(ax, 54, 86, 60, 86)
    
    draw_arrow(ax, 54, 74, 60, 81) 
    
    draw_arrow(ax, 54, 74, 60, 74) 
    
    draw_styled_box(ax, 60, 81, 14, 10, "Minimap2\n(Alignment)", color='white')
    
    draw_styled_box(ax, 60, 69, 14, 10, "ISFinder\n(Annotation)", color='#FFCCBC')

    draw_arrow(ax, 74, 86, 80, 80)
    draw_arrow(ax, 74, 74, 80, 80)
    
    draw_styled_box(ax, 80, 75, 12, 10, "Colocalization\nList (.csv)", color='#E8F5E9')

    draw_arrow(ax, 86, 75, 86, 56) 

    ax.add_patch(patches.Rectangle((2, 35), 96, 30, fc='#F5F5F5', ec='#B0BEC5', linestyle='--'))
    ax.text(5, 61, "Panel B: RMOBI Graph Construction (Flow: Right -> Left)", fontsize=12, fontweight='bold')

    draw_styled_box(ax, 80, 45, 12, 10, "Input\nCSV", color='#E8F5E9')
    draw_arrow(ax, 79, 50, 72, 50) 

    draw_styled_box(ax, 60, 45, 12, 10, "1. Create\nNodes", color='#FFF9C4', shape='circle')
    draw_arrow(ax, 59, 50, 54, 50) 

    draw_styled_box(ax, 42, 45, 12, 10, "2. Coloc.\nEdges", color='#FFF9C4')
    draw_arrow(ax, 41, 50, 36, 50) 

    draw_styled_box(ax, 24, 45, 12, 10, "3. Temporal\nEdges", color='#FFF9C4')
    draw_arrow(ax, 23, 50, 18, 50) 

    draw_styled_box(ax, 5, 45, 13, 10, "Graph\nJSON", color='#E8F5E9')

    draw_arrow(ax, 11, 44, 11, 23) 

    ax.add_patch(patches.Rectangle((2, 2), 96, 30, fc='#F5F5F5', ec='#B0BEC5', linestyle='--'))
    ax.text(5, 28, "Panel C: Interactive Visualization (D3.js)", fontsize=12, fontweight='bold')

    draw_styled_box(ax, 5, 12, 13, 10, "Load\nJSON", color='#E8F5E9')
    draw_arrow(ax, 19, 17, 26, 17) 

    ax.add_patch(patches.Rectangle((26, 5), 40, 20, fc='white', ec='black', lw=2))
    ax.text(46, 22, "CoNet Dashboard", ha='center', fontweight='bold')
    
    ax.plot([36, 46], [10, 15], 'k-', lw=1)
    ax.plot([46, 56], [15, 12], 'k-', lw=1)
    ax.add_patch(patches.Circle((36, 10), 1.5, fc='#EF5350'))
    ax.add_patch(patches.Circle((46, 15), 1.5, fc='#42A5F5'))
    ax.add_patch(patches.Circle((56, 12), 1.5, fc='#EF5350'))

    draw_arrow(ax, 67, 17, 72, 17) 
    draw_styled_box(ax, 72, 10, 20, 10, "Filters & Export\n(SVG/Stats)", color='#F8BBD0')

    plt.tight_layout()
    plt.savefig('figure1_zigzag_final_v2.png')
    plt.show()

if __name__ == "__main__":
    create_figure()