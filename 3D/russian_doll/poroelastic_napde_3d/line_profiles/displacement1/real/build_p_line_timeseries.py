# build_u_line_timeseries_with_coords.py
import glob, csv, re
import numpy as np

files = sorted(glob.glob("d_real_*.csv"))
assert files, "Non trovo file 'disp_ideal_*.csv' nella cartella corrente."

def read_csv_one_timestep(fname):
    """Legge un CSV (un timestep) e ritorna:
       (u0, u1, u2, time_value, arc_length, X, Y, Z)
       dove X,Y,Z possono essere None se non presenti.
    """
    u0, u1, u2, times = [], [], [], []
    s_list, x_list, y_list, z_list = [], [], [], []
    with open(fname, newline='') as f:
        rdr = csv.DictReader(f)
        # possibili nomi colonne
        # displacement
        k_u0, k_u1, k_u2 = "u:0", "u:1", "u:2"
        # tempo
        k_t_candidates = ["Time", "time"]
        # arc length
        k_s_candidates = ["Arc Length", "arc_length", "ARC_LENGTH"]
        # coordinate
        kx_candidates = ["Points:0", "X", "x"]
        ky_candidates = ["Points:1", "Y", "y"]
        kz_candidates = ["Points:2", "Z", "z"]

        for r in rdr:
            # time (se ripetuto per riga, prendo il primo alla fine)
            t = None
            for kt in k_t_candidates:
                if kt in r and r[kt] not in (None, ""):
                    t = float(r[kt]); break
            if t is not None:
                times.append(t)

            # u components
            u0.append(float(r[k_u0]))
            u1.append(float(r[k_u1]))
            u2.append(float(r[k_u2]))

            # arc length (se presente)
            s_val = None
            for ks in k_s_candidates:
                if ks in r and r[ks] not in (None, ""):
                    s_val = float(r[ks]); break
            s_list.append(s_val)

            # coordinates (se presenti)
            def read_first(keys):
                for k in keys:
                    if k in r and r[k] not in (None, ""):
                        return float(r[k])
                return None
            x_list.append(read_first(kx_candidates))
            y_list.append(read_first(ky_candidates))
            z_list.append(read_first(kz_candidates))

    # tempo del file: se presente prendi il primo; altrimenti dal nome file
    if times:
        tval = times[0]
    else:
        m = re.findall(r"(\d+\.?\d*)", fname)
        if not m:
            raise RuntimeError(f"Impossibile ricavare il tempo da '{fname}' (né da colonna Time né dal nome file).")
        tval = float(m[-1])

    # converti in array
    U0 = np.array(u0, dtype=float)
    U1 = np.array(u1, dtype=float)
    U2 = np.array(u2, dtype=float)
    S  = np.array(s_list, dtype=float) if any(v is not None for v in s_list) else None

    def to_array_or_none(lst):
        return (np.array(lst, dtype=float)
                if any(v is not None for v in lst) else None)
    X = to_array_or_none(x_list)
    Y = to_array_or_none(y_list)
    Z = to_array_or_none(z_list)

    return U0, U1, U2, tval, S, X, Y, Z

# Lettura di tutti i file
t_list = []
U0_list, U1_list, U2_list = [], [], []
S_ref = None
X_ref = Y_ref = Z_ref = None
Ns = None

for f in files:
    u0, u1, u2, t, S, X, Y, Z = read_csv_one_timestep(f)

    if Ns is None:
        Ns = len(u0)
        # salva il primo set di coordinate disponibile
        if X is not None and Y is not None and Z is not None:
            X_ref, Y_ref, Z_ref = X, Y, Z
        elif S is not None:
            S_ref = S
    else:
        if len(u0) != Ns or len(u1) != Ns or len(u2) != Ns:
            raise RuntimeError(f"Numero di campioni non coerente in '{f}'.")

        # verifica coerenza coordinate tra file
        def maxdiff(a, b):
            return np.max(np.abs(a - b)) if (a is not None and b is not None) else 0.0

        if X_ref is not None and X is not None:
            if maxdiff(X_ref, X) > 1e-10 or maxdiff(Y_ref, Y) > 1e-10 or maxdiff(Z_ref, Z) > 1e-10:
                raise RuntimeError(f"Coordinate X,Y,Z non coerenti tra file. Controlla linea/risoluzione. {f}")
        elif S_ref is not None and S is not None:
            if maxdiff(S_ref, S) > 1e-10:
                raise RuntimeError(f"Arc Length non coerente tra file. Controlla linea/risoluzione. {f}")

        # se nel primo file non c'erano XYZ ma qui sì, prendile
        if X_ref is None and X is not None and Y is not None and Z is not None:
            X_ref, Y_ref, Z_ref = X, Y, Z

    t_list.append(t)
    U0_list.append(u0)
    U1_list.append(u1)
    U2_list.append(u2)

# Ordina per tempo
t_arr = np.array(t_list)
order = np.argsort(t_arr)
t_sorted = t_arr[order]
U0 = np.vstack(U0_list)[order, :]  # (Nt, Ns)
U1 = np.vstack(U1_list)[order, :]
U2 = np.vstack(U2_list)[order, :]

# Costruisci le etichette col "meglio" che hai: priorità XYZ, altrimenti arc_length
if X_ref is not None and Y_ref is not None and Z_ref is not None:
    labels = [f"x={X_ref[k]:.12g} y={Y_ref[k]:.12g} z={Z_ref[k]:.12g}" for k in range(Ns)]
elif S_ref is not None:
    labels = [f"s={S_ref[k]:.12g}" for k in range(Ns)]
else:
    # fallback neutro
    labels = [f"sample_{k}" for k in range(Ns)]

def write_component_csv(outname, times, U, labels):
    with open(outname, "w", newline="") as f:
        w = csv.writer(f)
        header = ["time"] + labels
        w.writerow(header)
        for tk, row in zip(times, U):
            w.writerow([f"{tk:.12g}"] + [f"{val:.12g}" for val in row])

# Scrivi i tre file (uno per componente)
write_component_csv("disp_u0_timeseries.csv", t_sorted, U0, labels)
write_component_csv("disp_u1_timeseries.csv", t_sorted, U1, labels)
write_component_csv("disp_u2_timeseries.csv", t_sorted, U2, labels)

print("OK: scritto disp_u0_timeseries.csv, disp_u1_timeseries.csv, disp_u2_timeseries.csv")
print(f"    Nt = {len(t_sorted)} time step, Ns = {U0.shape[1]} campioni per linea")
