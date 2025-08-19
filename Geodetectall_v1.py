# %%
### This python code analyses the porosity of any structure in CIF format ###
####       Written By Alejandro DIAZ-MARQUEZ, CNRS Montpellier           ####
###                      Version 1.0 20250801                             ###

from __future__ import division
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import itertools
import math
import mdtraj as md
import sys
np.set_printoptions(threshold=sys.maxsize)
import warnings
warnings.filterwarnings("ignore")
import networkx as nx
from scipy.spatial.distance import pdist
from scipy.stats import kendalltau
import time
from ase.io import read, write
import os
import alphashape
import trimesh
from itertools import combinations
from scipy.spatial import cKDTree
import matplotlib.animation as animation
from PIL import Image
import shutil


# Show local directory
current_directory = os.getcwd()
print (current_directory)


# Define the input and output paths
folder_path = "./CIFs"
all_folder = "all"

# Load radius database
rad = np.loadtxt("radius",dtype='str',delimiter=" ")

# Radius of the particle insertion for detecting void
size = 0.8

# Precision of the 3D grid, distance between the grid points
g = 0.8

# If you want to expand the CIF system to detect bordering pores
expand = True

# If you want to generate a GIF about porosity evolution
video = True 

# Activate graph theory analysis
Graph_theory = True

# Minimun number of dots to consider as a pore, lower level can cause errors in alphashape
community_minimum_size = 50

def process_trajectory(subfolder_path, pdb_file, radius_data, box_size,size,g, community_minimum_size):
    traj1 = md.load(pdb_file)
    traj1.xyz = 10 * traj1.xyz  # Convert to angstrom
    box_size = 10 * traj1.unitcell_lengths[0]

    # Get coordinates and correlate with radius
    coords = []
    for atom in traj1.topology.atoms:
        atom_type = atom.element.symbol
        for j in range(len(radius_data)):
            if atom_type == radius_data[j, 1]:
                coord = traj1.xyz[0][atom.index]  # Get coordinates
                coords.append((coord, float(radius_data[j, 3]), atom_type))

    # Calculate box volume
    angles1 = traj1.unitcell_angles[0]
    alpha_rad, beta_rad, gamma_rad = map(np.radians, angles1)
    volume_box = (
        box_size[0] * box_size[1] * box_size[2]
        * np.sqrt(1 - np.cos(alpha_rad)**2 - np.cos(beta_rad)**2 - np.cos(gamma_rad)**2
                  + 2 * np.cos(alpha_rad) * np.cos(beta_rad) * np.cos(gamma_rad))
    )

    output_text_file = os.path.join(subfolder_path, "box_info.txt")
    with open(output_text_file, "w") as text_file:
        text_file.write(f"Box dimensions: {box_size}\n")
        text_file.write(f"Volume of the unit cell: {volume_box:.2f} Å³\n")

    # Prepare expanded coordinates
    coord_array, radius_array = zip(*[(c[0], c[1] / 100) for c in coords])
    coord_array = np.array(coord_array)
    radius_array = np.array(radius_array)
    print (len(coord_array))

    if expand:
        # Expand for periodic boundary conditions
        expanded_coords, expanded_radii = [], []
        for bx, by, bz in itertools.product(range(-1, 2), repeat=3):
            expanded_xyz = coord_array + np.array([bx * box_size[0], by * box_size[1], bz * box_size[2]]).T
            expanded_coords.append(expanded_xyz)
            expanded_radii.append(radius_array)

        coord_array = np.concatenate(expanded_coords)
        radius_array = np.concatenate(expanded_radii)

    # Define box limits (2 * box size)
    limits = np.array(box_size) + 2*size
    limits0 = np.array([0,0,0]) - 2*size

    # Filter coordinates that exceed limits
    mask = np.all((coord_array >= limits0) & (coord_array <= limits), axis=1)

    # Apply the mask to filter both coord_array and radius_array
    coord_array = coord_array[mask]
    radius_array = radius_array[mask]

    print ("doing points")
    # Generate grid points
    g2 = 2*g
    dots = [
        np.array([x, y, z])
        for x, y, z in itertools.product(np.arange(0, box_size[0], g),
                                         np.arange(0, box_size[1], g),
                                         np.arange(0, box_size[2], g))
        if not np.any(np.linalg.norm(coord_array - np.array([x, y, z]), axis=1) <= radius_array + size)
    ]
    dots = np.unique(dots, axis=0)

    # Save grid and plot
    grid_file = os.path.join(subfolder_path, "grid.xyz")
    np.savetxt(grid_file, dots, fmt="%.2f")

    # Community detection and alpha shape
    G_filtered = nx.Graph()
    for dot in dots:
        G_filtered.add_node(tuple(dot))

    print ("doing edges")
    # Add edges with periodicity
    kd_tree = cKDTree(dots)
    pairs = kd_tree.query_pairs(r=2 * g)
    for p1, p2 in pairs:
        dist = np.linalg.norm(dots[p1] - dots[p2])
        G_filtered.add_edge(tuple(dots[p1]), tuple(dots[p2]), weight=1 / dist)

    communities = nx.community.greedy_modularity_communities(G_filtered, resolution=0.1)
    # Find communities larger than community_minimum_size and store them directly in filtered_communities
    filtered_communities = [community for community in communities if len(community) >= community_minimum_size]
    num_communities = len(filtered_communities)
    print (num_communities)

    cmap = plt.get_cmap('viridis')  # You can choose any other colormap
    # Define colors and iterate through each community
    colors = [cmap(i / num_communities) for i in range(num_communities)]
    total_volume = 0
    for i, community in enumerate(filtered_communities):
        points = np.array(list(community))
        print (len(points))
        alpha_shape = alphashape.alphashape(points, alpha=0.5)
        volume = alpha_shape.volume
        total_volume += alpha_shape.volume
        color = colors[i % len(colors)]
        # Save community shape
        mesh = trimesh.Trimesh(vertices=alpha_shape.vertices, faces=alpha_shape.faces, vertex_colors=color)
        mesh.export(os.path.join(subfolder_path, f"community_{i}.obj"))

        if Graph_theory:
            # Calculate metrics for each community
            subgraph = G_filtered.subgraph(community)
            print ("clustering")
            clustering_coefficient = nx.algorithms.cluster.average_clustering(subgraph)
            print ("assortativity")
            assortativity = nx.degree_assortativity_coefficient(subgraph)

            print ("eccentricity")
            diameter = nx.diameter(subgraph)  # Directly computes the diameter

            print ("betweness")
            betweenness_centrality = nx.betweenness_centrality(subgraph, normalized=True)
            avg_betweenness_centrality = np.mean(list(betweenness_centrality.values()))
        else:
            clustering_coefficient = np.nan
            assortativity = np.nan
            diameter = np.nan
            avg_betweenness_centrality = np.nan         

        with open(output_text_file, "a") as text_file:
            text_file.write(f"Community {i} Metrics:\n")
            text_file.write(f"  - Number dots {len(community)}\n")
            text_file.write(f"  - Volume: {volume:.2f} Å³\n")
            text_file.write(f"  - Clustering Coefficient: {clustering_coefficient:.4f}\n")
            text_file.write(f"  - Assortativity: {assortativity:.4f}\n")
            text_file.write(f"  - Eccentricity: {diameter * g:.2f} Å\n")  # Scaled by 'g'
            text_file.write(f"  - Avg. Betweenness Centrality: {avg_betweenness_centrality:.4f}\n")
            text_file.write("-" * 40 + "\n")  # Separator for readability
    if video:
        # Define the z positions and the graph nodes
        z_positions = np.unique([node[2] for node in G_filtered.nodes])  # Unique z-values
        z_positions.sort()  # Ensure sorted order

        # Create a directory for saving frames
        frame_dir = os.path.join(subfolder_path, "frames")
        os.makedirs(frame_dir, exist_ok=True)

        # Define colormap for communities
        cmap = plt.get_cmap('tab20')

        # Generate frames for each z-position
        for i, z in enumerate(z_positions):
            fig, ax = plt.subplots(figsize=(12, 12))
            ax.set_title(f"z = {z:.2f} Å",fontsize=32)
            
            for idx, community in enumerate(filtered_communities):
                # Extract points in the current community
                points = np.array(list(community))

                # Filter points based on the current z-position
                mask = np.isclose(points[:, 2], z, atol=0.5)  # Adjust tolerance
                filtered_points = points[mask]

                # Plot filtered points if any exist
                if len(filtered_points) > 0:
                    ax.scatter(filtered_points[:, 0], filtered_points[:, 1], 
                            s=10, color=cmap(idx), label=f"Comm {idx + 1}")

            # Set axis limits and aspect ratio
            ax.set_xlim(0, box_size[0])
            ax.set_ylim(0, box_size[1])
            ax.set_aspect('equal')

            ax.set_xlabel('X-axis Label', fontsize=16)
            ax.set_ylabel('Y-axis Label', fontsize=16)

            # Set tick label size
            ax.tick_params(axis='both', which='major', labelsize=14)

            # Save frame as an image
            frame_file = os.path.join(frame_dir, f"frame_{i:03d}.png")
            plt.savefig(frame_file, dpi=300)
            plt.close(fig)  # Close to free memory

        print("Frames generated successfully!")

        # Parameters
        video_file = os.path.join(subfolder_path, "community_evolution_z_axis.gif")
        frame_rate = 2  # FPS
        frame_size = (600, 600)  # Match figure size (6x6 inches at 100 dpi)

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 12))

        # Load frames
        frames = []
        for i in range(len(z_positions)):
            frame_file = os.path.join(frame_dir, f"frame_{i:03d}.png")
            img = Image.open(frame_file).resize(frame_size)  # Resize for consistency
            frames.append(img)

        # Animation function
        def animate(i):
            ax.clear()  # Clear the previous frame
            ax.imshow(frames[i])
            ax.axis('off')

        ani = animation.FuncAnimation(fig, animate, frames=len(frames), interval=1000 / frame_rate)

        # Save as video
        #ani.save(video_file, fps=frame_rate)
        ani.save(video_file, writer='pillow', fps=frame_rate)

        print(f"Video saved as {video_file}")
    
# Process each CIF file
for file_name in sorted(os.listdir(folder_path)):
    if file_name.endswith(".cif"):
        base_name = os.path.splitext(file_name)[0]
        subfolder_path = os.path.join(all_folder, base_name)
        # Remove the base_name folder if it exists
        if os.path.exists(subfolder_path):
            shutil.rmtree(subfolder_path)
        os.makedirs(subfolder_path, exist_ok=True)
        print (file_name)
        cif_file = os.path.join(folder_path, file_name)
        structure = read(cif_file)
        pdb_file = os.path.join(subfolder_path, f"{base_name}.pdb")
        write(pdb_file, structure)

        box_size = structure.get_cell_lengths_and_angles()[:3]
        process_trajectory(subfolder_path, pdb_file, rad, box_size,size,g, community_minimum_size)



