import json, sys
f = sys.argv[1]
d = json.load(open(f))
plddt = d["plddt"]
n = len(plddt)
def mean(a,b):
    seg = plddt[a-1:b]
    return round(sum(seg)/len(seg),1)
print("Longitud:", n)
print("pLDDT global:", round(sum(plddt)/n,1))
print("Conector N-terminal 1-69:", mean(1,69))
print("  -- subtramo 30-60 (LCR):", mean(30,60))
print("Dominio periplasmico 70-362:", mean(70,362))
print("Nucleo TM 363-542:", mean(363,542))
for p in [388,358,360,451,427]:
    print(f"pLDDT residuo {p}: {round(plddt[p-1],1)}")
lo50 = [i+1 for i,v in enumerate(plddt) if v<50]
lo70 = sum(1 for v in plddt if v<70)
print("Residuos pLDDT<50:", len(lo50), "->", lo50[:40])
print("Residuos pLDDT<70:", lo70)
