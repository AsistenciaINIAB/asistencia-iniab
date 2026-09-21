import importlib, sys
print("python", sys.version.split()[0])
for m in ["numpy", "matplotlib", "Bio"]:
    try:
        mod = importlib.import_module(m)
        print(m, getattr(mod, "__version__", "?"))
    except Exception as e:
        print(m, "MISSING", e)
