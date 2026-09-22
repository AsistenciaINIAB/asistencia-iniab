import json, sys
d = json.load(open(sys.argv[1]))
pae = d.get("pae") or d.get("predicted_aligned_error")
n = len(pae)
def block_mean(ra, rb, ca, cb):
    vals = [pae[i][j] for i in range(ra-1,rb) for j in range(ca-1,cb)]
    return round(sum(vals)/len(vals),1)
print("max_pae:", d.get("max_pae"))
print("PAE intra dominio periplasmico (70-362):", block_mean(70,362,70,362))
print("PAE intra nucleo TM (363-542):", block_mean(363,542,363,542))
print("PAE inter-dominio (peri x TM):", block_mean(70,362,363,542))
