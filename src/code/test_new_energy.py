from TPP.API.energy import Energy
from TPP.API.visualizer import draw_heatmap
import json
import numpy as np
import pandas as pd

def get_protein_pairs_matrix(cliques):
    arr = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]
    ref = {
        "GLY": 0,
        "PRO": 1,
        "ASP": 2,
        "GLU": 3,
        "LYS": 4,
        "ARG": 5,
        "HIS": 6,
        "SER": 7,
        "THR": 8,
        "ASN": 9,
        "GLN": 10,
        "ALA": 11,
        "MET": 12,
        "TYR": 13,
        "TRP": 14,
        "VAL": 15,
        "ILE": 16,
        "LEU": 17,
        "PHE": 18,
        "CYS": 19
    }
    for clique in cliques:
        for i in range(len(clique)):
            for j in range(i + 1, len(clique)):
                if ref[clique[i]] == ref[clique[j]]:
                    arr[ref[clique[i]]][ref[clique[j]]] += 1
                else:
                    arr[ref[clique[i]]][ref[clique[j]]] += 1
                    arr[ref[clique[j]]][ref[clique[i]]] += 1
    return np.array(arr)


def generate_cliques(data):
    cliques = []
    for i in data:
        for j in range(i[1]):
            cliques.append(eval(i[0]))
    return cliques

with open("clique_data.json", "r") as file:
    data = json.load(file)
    total_data = data["total"]
    water_data = data["water"]
    interface_data = data["interface"]
    hydrophobic_data = data["hydrophobic"]

    total_cliques = generate_cliques(total_data)
    water_cliques = generate_cliques(water_data)
    interface_cliques = generate_cliques(interface_data)
    hydrophobic_cliques = generate_cliques(hydrophobic_data)
    total_map = get_protein_pairs_matrix(total_cliques)
    water_map = get_protein_pairs_matrix(water_cliques)
    interface_map = get_protein_pairs_matrix(interface_cliques)
    hydrophobic_map = get_protein_pairs_matrix(hydrophobic_cliques)

    E_total = Energy()
    E_water = Energy()
    E_interface = Energy()
    E_hydrophobic = Energy()

    E_total.update_static_total_pairs_table(total_map)
    E_water.update_static_total_pairs_table(water_map)
    E_interface.update_static_total_pairs_table(interface_map)
    E_hydrophobic.update_static_total_pairs_table(hydrophobic_map)

    E_total.update_epair_values2()
    total_e_pair_table = E_total.get_static_epair_table()

    E_total.update_epair_values2(layer_map=E_water.get_static_total_pairs_table())
    water_e_pair_table = E_total.get_static_epair_table()

    E_total.update_epair_values2(layer_map=E_interface.get_static_total_pairs_table())
    interface_e_pair_table = E_total.get_static_epair_table()

    E_total.update_epair_values2(layer_map=E_hydrophobic.get_static_total_pairs_table())
    hydrophobic_e_pair_table = E_total.get_static_epair_table()

    AAs = ["G", "P", "D", "E", "K", "R", "H", "S", "T", "N", "Q", "A", "M", "Y", "W", "V", "I", "L", "F", "C"]

    draw_heatmap("E_Total", total_e_pair_table, AAs, AAs, "gist_rainbow_r")
    draw_heatmap("E_Water", water_e_pair_table, AAs, AAs, "gist_rainbow_r")
    draw_heatmap("E_Interface", interface_e_pair_table, AAs, AAs, "gist_rainbow_r")
    draw_heatmap("E_Hydrophobic", hydrophobic_e_pair_table, AAs, AAs, "gist_rainbow_r")



# Analyze unusual groupings from pairwise plot start with pos to pos find chain ids, res ids, etc to further analyze later
# write energy function draft
# think about how to visualize higher order potentials/ cliques

        #1 axis pair of cliquyes lowest in sequence for 4 body clique is one idea