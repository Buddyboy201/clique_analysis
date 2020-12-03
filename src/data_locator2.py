import json
import pandas as pd


# goal: create excel file out of residue pairings within cliques data

def get_cliques_with_interaction(aai, aaj, cliques):
    interactions = []
    for clique in cliques:
        if aai in clique[0] and aaj in clique[0]:
            interactions.append(clique)
    return interactions

def get_all_interactions(aai, aaj, all_cliques):
    interactions = {2:[], 3:[], 4:[], 5:[], 6:[]}
    for size in all_cliques:
        data = all_cliques[size]
        interactions[int(size)] = get_cliques_with_interaction(aai, aaj, data)
    return interactions

def get_all_interaction_pairs(all_cliques, ref = ["GLY", "PRO", "ASP", "GLU", "LYS", "ARG", "HIS", "SER", "THR", "ASN", "GLN", "ALA", "MET", "TYR", "TRP", "VAL", "ILE", "LEU", "PHE", "CYS"]):
    interactions = {}
    for i in range(len(ref)):
        for j in range(i, len(ref)):
            aai = ref[i]
            aaj = ref[j]
            arr = sorted([aai, aaj])
            val = "{}|{}".format(arr[0], arr[1])
            interactions[val] = get_all_interactions(aai, aaj, all_cliques)
    return interactions


def get_layer_cliques_for_res_group(interactions):
    total_cliques = {2: [], 3: [], 4: [], 5: [], 6: []}
    for combo in interactions:
        for size in interactions[combo]:
            if len(interactions[combo][size]) > 0:
                for i in interactions[combo][size]:
                    total_cliques[size].append(i)
    return total_cliques

def get_dataframe_for_cliques(total_cliques):
    all_total_cliques = data = {"size_2": [], "freq_2": [], "size_3": [], "freq_3": [], "size_4": [], "freq_4": [], "size_5": [], "freq_5": [], "size_6": [], "freq_6": []}
    longest = 0
    for size in total_cliques:
        total_cliques[size] = sorted(total_cliques[size], reverse=True, key=lambda x: x[1])
        longest = max(longest, len(total_cliques[size]))
    for size in total_cliques:
        total_cliques[size] += [["", 0]] * (longest - len(total_cliques[size]))

    for size in total_cliques:
        if len(total_cliques[size]) > 0:
            for i in total_cliques[size]:
                all_total_cliques["size_{}".format(size)].append(i[0])
                all_total_cliques["freq_{}".format(size)].append(i[1])
    df = pd.DataFrame(all_total_cliques)
    return df


def generate_cliques_data_for_res_group(cliques, ref):
    return get_dataframe_for_cliques(get_layer_cliques_for_res_group(get_all_interaction_pairs(cliques, ref=ref)))

def get_res_group_spreadsheet(global_cliques, ref, file_name):
    all_cliques = generate_cliques_data_for_res_group(global_cliques["total"], ref)
    water_cliques = generate_cliques_data_for_res_group(global_cliques["water"], ref)
    interface_cliques = generate_cliques_data_for_res_group(global_cliques["interface"], ref)
    hydrophobic_cliques = generate_cliques_data_for_res_group(global_cliques["hydrophobic"], ref)
    with pd.ExcelWriter("{}.xlsx".format(file_name)) as writer:
        all_cliques.to_excel(writer, sheet_name="all_layers")
        water_cliques.to_excel(writer, sheet_name="water")
        interface_cliques.to_excel(writer, sheet_name="interface")
        hydrophobic_cliques.to_excel(writer, sheet_name="hydrophobic")

with open("size_sorted_clique_data.json", "r") as data_file:
    size_ref = json.load(data_file)
    get_res_group_spreadsheet(size_ref, ["ARG", "LYS"], "sorted_layerwise_ARG_LYS_2_cliques")
