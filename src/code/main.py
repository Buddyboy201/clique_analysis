import json
import pandas as pd

def format_results_full_clique_list(layer, ref, file_name):
    ranking = ref[layer]
    top = [[], [], [], [], [], []]

    if ranking.get("2") is None:
        ranking["2"] = []
    for i in ranking["2"]:
        top[0].append(i)

    if ranking.get("3") is None:
        ranking["3"] = []
    for i in ranking["3"]:
        top[1].append(i)

    if ranking.get("4") is None:
        ranking["4"] = []
    for i in ranking["4"]:
        top[2].append(i)

    if ranking.get("5") is None:
        ranking["5"] = []
    for i in ranking["5"]:
        top[3].append(i)

    if ranking.get("6") is None:
        ranking["6"] = []
    for i in ranking["6"]:
        top[4].append(i)

    all_cliques = []
    for i in ranking.values():
        all_cliques += i
    for i in sorted(all_cliques, reverse=True, key=lambda x: x[1]):
        top[5].append(i)


    max_length = max([len(i) for i in top])
    top[0] += [["", 0]] * (max_length - len(top[0]))
    top[1] += [["", 0]] * (max_length - len(top[1]))
    top[2] += [["", 0]] * (max_length - len(top[2]))
    top[3] += [["", 0]] * (max_length - len(top[3]))
    top[4] += [["", 0]] * (max_length - len(top[4]))
    top[5] += [["", 0]] * (max_length - len(top[5]))

    return top

def get_dataframe(top):
    data = {"size_2": [], "freq_2": [], "size_3": [], "freq_3": [], "size_4": [], "freq_4": [], "size_5": [], "freq_5": [], "size_6": [], "freq_6": [], "size_all": [], "freq_all": []}
    for clique, freq in top[0]:
        data["size_2"].append(str(clique))
        data["freq_2"].append(freq)
    for clique, freq in top[1]:
        data["size_3"].append(str(clique))
        data["freq_3"].append(freq)
    for clique, freq in top[2]:
        data["size_4"].append(str(clique))
        data["freq_4"].append(freq)
    for clique, freq in top[3]:
        data["size_5"].append(str(clique))
        data["freq_5"].append(freq)
    for clique, freq in top[4]:
        data["size_6"].append(str(clique))
        data["freq_6"].append(freq)
    for clique, freq in top[5]:
        data["size_all"].append(str(clique))
        data["freq_all"].append(freq)

    df = pd.DataFrame(data)
    return df




with open("size_sorted_clique_data.json", "r") as data_file:
    size_ref = json.load(data_file)
    top_total = format_results_full_clique_list("total", size_ref, "")
    top_water = format_results_full_clique_list("water", size_ref, "")
    top_interface = format_results_full_clique_list("interface", size_ref, "")
    top_hydrophobic = format_results_full_clique_list("hydrophobic", size_ref, "")

    df_total = get_dataframe(top_total)
    df_water = get_dataframe(top_water)
    df_interface = get_dataframe(top_interface)
    df_hydrophobic = get_dataframe(top_hydrophobic)

    with pd.ExcelWriter("sorted_clique_rankings_by_layer_full_ranking.xlsx") as writer:
        df_total.to_excel(writer, sheet_name="all_layers")
        df_water.to_excel(writer, sheet_name="water")
        df_interface.to_excel(writer, sheet_name="interface")
        df_hydrophobic.to_excel(writer, sheet_name="hydrophobic")






# share top 50 by layer and diff sizes and how many unique counts there are (e.g. ["ALA, "LEU, "VAL"] is one)
# reimplement Energy() class with layers included and updated vladimir formula
#









#deep learning for de novo protein design session 4B tuesday aug 4

