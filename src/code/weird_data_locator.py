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
        for j in range(i+1, len(ref)):
            aai = ref[i]
            aaj = ref[j]
            arr = sorted([aai, aaj])
            val = "{}|{}".format(arr[0], arr[1])
            interactions[val] = get_all_interactions(aai, aaj, all_cliques)
    return interactions



def get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, layer):
    total_cliques = {2: [], 3: [], 4: [], 5: [], 6: []}
    for combo in polar_nonpolar_interactions:
        for size in polar_nonpolar_interactions[combo][layer]:
            if len(polar_nonpolar_interactions[combo][layer][size]) > 0:
                for i in polar_nonpolar_interactions[combo][layer][size]:
                    total_cliques[size].append(i)
    return total_cliques

def get_dataframe_for_polar_nonpolar(total_cliques):
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

with open("size_sorted_clique_data.json", "r") as data_file:
    size_ref = json.load(data_file)
    all_cliques = get_all_interaction_pairs(size_ref["total"])
    water_cliques = get_all_interaction_pairs(size_ref["water"])
    interface_cliques = get_all_interaction_pairs(size_ref["interface"])
    hydrophobic_cliques = get_all_interaction_pairs(size_ref["hydrophobic"])
    arr = {"total": all_cliques, "water": water_cliques, "interface": interface_cliques, "hydrophobic": hydrophobic_cliques}
    #with open("pairwise_interactions_within_cliques.json", "w") as dump_file:
    #    json.dump(arr, dump_file)


    print()



    polar = ["GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS"]
    non_polar = ["ALA", "ILE", "LEU", "MET", "PHE", "VAL", "PRO", "GLY"]
    combos = []
    for aai in polar:
        for aaj in non_polar:
            arr = sorted([aai, aaj])
            val = "{}|{}".format(arr[0], arr[1])
            combos.append(val)

    polar_nonpolar_interactions = {}
    for i in combos:
        total = all_cliques[i]
        water = water_cliques[i]
        interface  = interface_cliques[i]
        hydrophobic = hydrophobic_cliques[i]
        ref = {"total":total, "water":water, "interface":interface, "hydrophobic": hydrophobic}
        polar_nonpolar_interactions[i] = ref
    #print(polar_nonpolar_interactions["ALA|GLN"]["total"][3])


    #hydrophobic
    #arg_lysine
    #asp_glu
    #get all arginines
    #one or more arginines
    # look at acidic bases against themselves


    val = "ARG|LYS"
    total = all_cliques[val]
    water = water_cliques[val]
    interface = interface_cliques[val]
    hydrophobic = hydrophobic_cliques[val]
    arg_lys_interactions = {}
    ref = {"total": total, "water": water, "interface": interface, "hydrophobic": hydrophobic}
    arg_lys_interactions[val] = ref
    hydrophobic_df = get_dataframe_for_polar_nonpolar(get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, "hydrophobic"))
    with pd.ExcelWriter("arg_lys_hydro.xlsx") as writer:
        hydrophobic_df.to_excel(writer, sheet_name="hydrophobic")





    total_df = get_dataframe_for_polar_nonpolar(get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, "total"))
    water_df = get_dataframe_for_polar_nonpolar(get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, "water"))
    interface_df = get_dataframe_for_polar_nonpolar(get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, "interface"))
    hydrophobic_df = get_dataframe_for_polar_nonpolar(get_layer_cliques_for_polar_nonpolar(polar_nonpolar_interactions, "hydrophobic"))

    #with pd.ExcelWriter("sorted_clique_rankings_by_layer_for_polar_hydrophobic_interactions.xlsx") as writer:
    #    total_df.to_excel(writer, sheet_name="all_layers")
    #    water_df.to_excel(writer, sheet_name="water")
    #    interface_df.to_excel(writer, sheet_name="interface")
    #    hydrophobic_df.to_excel(writer, sheet_name="hydrophobic")



    # want:
    #   -top 50 global non polar polar combinations
    #   -top 50 layer wise non polar polar combinations
    #   -
