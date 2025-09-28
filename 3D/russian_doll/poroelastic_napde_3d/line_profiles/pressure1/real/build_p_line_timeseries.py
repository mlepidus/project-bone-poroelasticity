# build_p_line_timeseries.py
import glob, csv, re
import numpy as np

files = sorted(glob.glob("p_real_*.csv"))
assert files, "Non trovo file line_profile_*.csv"

S = None
rows = []  # (time, p_at_s array)

def read_csv(fname):
    with open(fname, newline='') as f:
        rdr = csv.DictReader(f)
        arc, pres, times = [], [], []
        for r in rdr:
            arc.append(float(r.get("Arc Length") or r.get("arc_length") or r["Arc Length"]))
            # prova campi Pressure possibili
            p = (r.get("Pressure") or r.get("p") or r.get("pressure"))
            pres.append(float(p))
            t = r.get("Time") or r.get("time")
            if t is not None:
                times.append(float(t))
        return np.array(arc), np.array(pres), (None if not times else times[0])

# leggi tutti i file
t_vals = []
p_mat = []
for f in files:
    s, p, t = read_csv(f)
    if S is None:
        S = s
    else:
        if len(s)!=len(S) or np.max(np.abs(s - S)) > 1e-10:
            raise RuntimeError(f"Arc Length non coerente tra file. Controlla Resolution/linea. {f}")
    if t is None:
        # ricava t dal nome file (numeri come #### o ####.##)
        m = re.findall(r"(\d+\.?\d*)", f)
        if not m: raise RuntimeError(f"Impossibile ricavare il tempo da {f}")
        t = float(m[-1])
    t_vals.append(t)
    p_mat.append(p)

# ordina per tempo
order = np.argsort(np.array(t_vals))
t_sorted = np.array(t_vals)[order]
P = np.vstack(p_mat)[order,:]  # shape (Nt, Ns)

# scrivi pivot CSV: header time,s0,s1,... (s in metri)
out = "p_real_tot.csv"
with open(out, "w", newline='') as f:
    w = csv.writer(f)
    header = ["time"] + [f"{si:.12g}" for si in S]
    w.writerow(header)
    for tk, row in zip(t_sorted, P):
        w.writerow([f"{tk:.12g}"] + [f"{val:.12g}" for val in row])

print(f"OK: scritto {out} con {len(t_sorted)} tempi e {len(S)} campioni s.")
