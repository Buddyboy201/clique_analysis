import json
import pandas as pd

def format_results(layer, ref, file_name):
    ranking = ref[layer]
    top = [[], [], [], [], [], []]

    if ranking.get("2") is None:
        ranking["2"] = []
    for i in ranking["2"][:50]:
        top[0].append(i)

    top[0] += [["", 0]]*(50-len(top[0]))

    if ranking.get("3") is None:
        ranking["3"] = []
    for i in ranking["3"][:50]:
        top[1].append(i)

    top[1] += [["", 0]] * (50 - len(top[1]))

    if ranking.get("4") is None:
        ranking["4"] = []
    for i in ranking["4"][:50]:
        top[2].append(i)

    top[2] += [["", 0]] * (50 - len(top[2]))

    if ranking.get("5") is None:
        ranking["5"] = []
    for i in ranking["5"][:50]:
        top[3].append(i)

    top[3] += [["", 0]] * (50 - len(top[3]))

    if ranking.get("6") is None:
        ranking["6"] = []
    for i in ranking["6"][:50]:
        top[4].append(i)

    top[4] += [["", 0]] * (50 - len(top[4]))

    all_cliques = []
    for i in ranking.values():
        all_cliques += i
    for i in sorted(all_cliques, reverse=True, key=lambda x: x[1])[:50]:
        top[5].append(i)

    top[5] += [["", 0]] * (50 - len(top[5]))

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
    top_total = format_results("total", size_ref, "")
    top_water = format_results("water", size_ref, "")
    top_interface = format_results("interface", size_ref, "")
    top_hydrophobic = format_results("hydrophobic", size_ref, "")

    df_total = get_dataframe(top_total)
    df_water = get_dataframe(top_water)
    df_interface = get_dataframe(top_interface)
    df_hydrophobic = get_dataframe(top_hydrophobic)

    with pd.ExcelWriter("sorted_clique_rankings_by_layer.xlsx") as writer:
        df_total.to_excel(writer, sheet_name="all_layers")
        df_water.to_excel(writer, sheet_name="water")
        df_interface.to_excel(writer, sheet_name="interface")
        df_hydrophobic.to_excel(writer, sheet_name="hydrophobic")













#deep learning for de novo protein design session 4B tuesday aug 4