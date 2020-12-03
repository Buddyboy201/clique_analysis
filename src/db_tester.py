from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String
engine = create_engine('sqlite:///clique_traceback_v2.db', echo = True)
meta = MetaData()
cliques_table = Table(
    'cliques', meta,
    Column("id", Integer, primary_key=True),
    Column("size", Integer),
    Column("clique", String),
    Column("resid", String),
    Column("layerinfo", String),
    Column("pdbname", String)
)
conn = engine.connect()
stmt_duplicate = "SELECT clique FROM cliques WHERE size={} AND instr(clique, '{}');"
stmt_all = "SELECT clique FROM cliques WHERE size={} AND (instr(clique, '{}') AND instr(clique, '{}'));"

#stmt = "SELECT * FROM cliques WHERE size=3 AND (instr(clique, 'ARG') AND instr(clique, 'GLU') AND instr(clique, 'LYS'));"

#for row in conn.execute(stmt):
   # print(row)

all_hash = {"size_2":{}, "size_3":{}, "size_4":{}, "size_5":{}}

def fill_hash(conn, stmt, hash):
    for row in conn.execute(stmt):
        val = row[0]
        if hash.get(val) is None: hash[val] = 0
        hash[val] += 1

def fill_hash_duplicate(conn, stmt, hash, exclude, res):
    for row in conn.execute(stmt):
        val = row[0]
        flag = True
        for i in exclude:
            if i in val: flag = False
        if flag and val.count(res) > 1:
            if hash.get(val) is None: hash[val] = 0
            hash[val] += 1

basic = ["ARG", "LYS", "HIS"]
def get_pairs(residues):
    new_list = []
    for i in range(len(residues)):
        for j in range(i+1, len(residues)):
            new_list.append(list(sorted([residues[i], residues[j]])))
    return new_list

for size in range(2, 6):
    for combo in get_pairs(basic):
        stmt = stmt_all.format(size, combo[0], combo[1])
        if combo[0] != combo[1]:
            #new_basic = basic
            #new_basic.remove(combo[0])
            #fill_hash_duplicate(conn, stmt_duplicate.format(size, combo[0]), all_hash, new_basic, combo[0])
            fill_hash(conn, stmt, all_hash["size_{}".format(size)])

for size in range(2, 6):
    for res in basic:
        stmt = stmt_duplicate.format(size, res)
        new_basic = basic
        new_basic.remove(res)
        fill_hash_duplicate(conn, stmt, all_hash["size_{}".format(size)], new_basic, res)
print(all_hash)
