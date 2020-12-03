import json

#top 50 all cliques in category
#top 50 size n till n=2 in category

def update_ref(data, ref):
    for clique_str, count in data:
        clique = eval(clique_str)
        if ref.get(len(clique)) is None:
            ref[len(clique)] = []
        ref[len(clique)].append((clique, count))
    return ref

with open("clique_data.json", "r") as data_file:
    data = json.load(data_file)
    size_ref = {"total":{}, "water":{}, "interface":{}, "hydrophobic":{}}

    size_ref["total"] = update_ref(data["total"], {})
    size_ref["water"] = update_ref(data["water"], {})
    size_ref["interface"] = update_ref(data["interface"], {})
    size_ref["hydrophobic"] = update_ref(data["hydrophobic"], {})
    with open("size_sorted_clique_data.json", "w") as dump_file:
        json.dump(size_ref, dump_file)

