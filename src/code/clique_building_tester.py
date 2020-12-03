import sys
import os
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from TPP.API.top_pro_pack import TPP_Engine
from TPP.API.centroid_protein import CentroidProtein

pdb_dir = r"C:\test_proteins\Menv_color\Menv_color"
out_dir = r"C:\test_proteins\Menv_log_out\Menv_log_out"

#for file_name in [f for f in os.listdir(pdb_dir) if ".pdb" in f]:
    #os.rename(pdb_dir + "\\" + file_name, pdb_dir + "\\" + file_name[file_name.find("nor_")+len("nor_"):file_name.find("_r")]+".pdb")

#for file_name in [f for f in os.listdir(out_dir) if ".out" in f and "id" not in f]:
    #os.rename(out_dir + "\\" + file_name, out_dir + "\\" + file_name[:file_name.find("_Menv")]+".out")

pdb_file_names = [f for f in os.listdir(pdb_dir) if ".pdb" in f]
out_file_names = [f for f in os.listdir(out_dir) if ".out" in f and "id" not in f]

pdb_file_names.sort()
out_file_names.sort()
bad_pos = 390
pdb_file_names.pop(bad_pos)

def get_filtered_out_lines(out_file):
    with open(out_file, "rt") as file:
        lines = file.readlines()
        return [[i for i in line.split(" ") if i != ""] for line in lines if line.split(" ")[0].strip(" ") == "2016Menv"]



min_hydrophobic_residues = 14

water_cliques = []
interface_cliques = []
hydrophobic_cliques = []

def get_layer(res, layer_ref):
    return layer_ref[res.resid + 1]

def get_layer_resid(resid, layer_ref):
    return layer_ref[resid + 1]


def filter_water(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [1, 6]]) == 0:
            new_cliques.append(clique)
    return new_cliques


def filter_interface(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [2, 5]]) == 0:
            new_cliques.append(clique)
    return new_cliques


def filter_hydrophobic(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [3, 4]]) == 0:
            new_cliques.append(clique)
    return new_cliques


def sort_cliques_and_convert_to_str(cliques):
    for i in range(len(cliques)):
        # cliques[i].sort(lambda x: x.name)
        cliques[i] = [res.name for res in cliques[i]]
        cliques[i].sort()
        cliques[i] = repr(cliques[i])
    return cliques

def get_clique_with_names_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([i.name for i in clique])

def get_clique_with_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.resid) for i in clique])

def get_clique_with_old_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.old_resid) for i in clique])

def get_clique_layer_info_only(clique, layer_ref):
    clique.sort(key=lambda x: x.name)
    resids = [res.resid for res in clique]
    return ";".join([str(get_layer_resid(resid, layer_ref)) for resid in resids])


def update_clique_counts(cliques, hashmap):
    for clique in cliques:
        # print(clique)
        if hashmap.get(repr(clique)) is None:
            hashmap[clique] = 0
        hashmap[clique] += 1
    return hashmap

engine = create_engine("sqlite:///clique_traceback_v2.db", echo=True)

meta = MetaData()
cliques_table = Table(
    'cliques', meta,
    Column("id", Integer, primary_key=True),
    Column("size", Integer),
    Column("clique", String),
    Column("resid", String),
    Column("oldresid", String),
    Column("layerinfo", String),
    Column("pdbname", String)
)
meta.create_all(engine)

conn = engine.connect()
def insert_clique_into_db(clique, pdb_name, layer_ref, conn, table):
    ins = table.insert().values(size=len(clique), clique=get_clique_with_names_only(clique),
                                resid=get_clique_with_resid_only(clique), oldresid=get_clique_with_old_resid_only(clique), layerinfo=get_clique_layer_info_only(clique, layer_ref), pdbname=pdb_name)
    result = conn.execute(ins)
    return result

good_proteins = []
bad_proteins = []
exceptions = []
Engine = TPP_Engine()
clique_counts_total = {}
clique_counts_water = {}
clique_counts_interface = {}
clique_counts_hydrophobic = {}
log_num = 21
file = open("layer_clique_ranking_log_{}.txt".format(log_num), "w")
for out_file_name, pdb_file_name in list(zip(out_file_names, pdb_file_names)):
    pdb_name = pdb_file_name[:pdb_file_name.find(".pdb")]
    P = Engine.init_protein("", pdb_name, pdb_dir + "\\" + pdb_file_name)
    content = get_filtered_out_lines(out_dir + "\\" + out_file_name)
    hydrophobic_count = 0
    layer_ref = {}
    for line in content:
        res = line[2].strip(" ")
        id = int(line[1].strip(" "))
        layer = int(line[4].strip(" "))
        layer_ref[id] = layer
        if layer == 3 or layer == 4:
            hydrophobic_count += 1
    # print(hydrophobic_count)
    if type(P) is Exception:
        print(str(P))
        exceptions.append(pdb_name + " " + str(P))
    elif hydrophobic_count >= min_hydrophobic_residues:
        print("Triggered")
        good_proteins.append(P)
        cliques = P.centroid_cliques
        for clique in P.centroid_cliques:
            insert_clique_into_db(clique, P.name, layer_ref, conn, cliques_table)
        # print(cliques)
        #cliques_water = filter_water(P.centroid_cliques, layer_ref)
        # print(cliques_water)
        #cliques_interface = filter_interface(P.centroid_cliques, layer_ref)
        #cliques_hydrophobic = filter_hydrophobic(P.centroid_cliques, layer_ref)
        #cliques = sort_cliques_and_convert_to_str(cliques)
        #cliques_water = sort_cliques_and_convert_to_str(cliques_water)
        #cliques_interface = sort_cliques_and_convert_to_str(cliques_interface)
        #cliques_hydrophobic = sort_cliques_and_convert_to_str(cliques_hydrophobic)
        #clique_counts_total = update_clique_counts(cliques, clique_counts_total)
        #clique_counts_water = update_clique_counts(cliques_water, clique_counts_water)
        #clique_counts_interface = update_clique_counts(cliques_interface, clique_counts_interface)
        #clique_counts_hydrophobic = update_clique_counts(cliques_hydrophobic, clique_counts_hydrophobic)

    else:
        print("baddddd")
        bad_proteins.append(P)
file.write("bad_proteins:\n")
for i in bad_proteins:
    file.write("\t" + i.name + "\n")
file.write("exceptions:\n")
for i in exceptions:
    file.write("\t" + i + "\n")
file.close()
#sorted_total_counts = sorted(clique_counts_total.items(), key=lambda i: i[1], reverse=True)
#sorted_water_counts = sorted(clique_counts_water.items(), key=lambda i: i[1], reverse=True)
#sorted_interface_counts = sorted(clique_counts_interface.items(), key=lambda i: i[1], reverse=True)
#sorted_hydrophobic_counts = sorted(clique_counts_hydrophobic.items(), key=lambda i: i[1], reverse=True)
#data = {}
#data["total"] = sorted_total_counts
#data["water"] = sorted_water_counts
#data["interface"] = sorted_interface_counts
#data["hydrophobic"] = sorted_hydrophobic_counts
# with open("clique_data.json", "w") as data_file:
# json.dump(data, data_file)