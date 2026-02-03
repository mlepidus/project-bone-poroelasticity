import subprocess
import re
import matplotlib.pyplot as plt
import numpy as np
import os

# --- CONFIGURAZIONE ---
DATA_FILE = "data_2d_square_convergence.txt"
EXECUTABLE = "./main_coupled"
# Discretizzazioni da testare (es. 4 cubi per lato, 8, 16...)
MESH_SIZES = [4, 8, 16] 
# Se il calcolo diventa troppo lento, fermati a 16 o 24. 
# 32 in 3D significa 32^3 * dof = enorme per un portatile.

def update_config_file(filename, N):
    """Legge il file .txt e aggiorna la discretizzazione spaziale"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    with open(filename, 'w') as f:
        for line in lines:
            # Sostituisce i parametri di discretizzazione
            if "spatialDiscretizationX" in line:
                f.write(f"        spatialDiscretizationX = {N}\n")
            elif "spatialDiscretizationY" in line:
                f.write(f"        spatialDiscretizationY = {N}\n")
            elif "spatialDiscretizationZ" in line:
                f.write(f"        spatialDiscretizationZ = {N}\n")
            # Assicuriamoci che non usi mesh esterne
            elif "meshExternal" in line:
                f.write("        meshExternal = none\n")
            else:
                f.write(line)

def run_simulation():
    """Lancia il codice C++ e cattura l'output"""
    print(f"Running {EXECUTABLE} with {DATA_FILE}...")
    try:
        # Esegue il comando e cattura lo stdout
        result = subprocess.run([EXECUTABLE, "-f", DATA_FILE], capture_output=True, text=True)
        output = result.stdout
        
        # Cerca la riga dell'errore usando RegEx
        # Cerca qualcosa tipo: Error: [pressure, displacement] = [0.123, 0.456]
        # Adatta questa regex se il tuo output è leggermente diverso
        match = re.search(r"Error:\s*\[.*\]\s*=\s*\[([0-9\.eE\-\+]+),\s*([0-9\.eE\-\+]+)\]", output)
        
        if match:
            err_p = float(match.group(1))
            err_u = float(match.group(2))
            return err_p, err_u
        else:
            print("Errore: Non ho trovato la riga 'Error: ...' nell'output.")
            print("Output parziale:\n", output[-500:]) # Stampa le ultime righe
            return None, None
            
    except Exception as e:
        print(f"Errore durante l'esecuzione: {e}")
        return None, None

# --- MAIN LOOP ---
h_values = []
errors_p = []
errors_u = []

print("=== INIZIO TEST DI CONVERGENZA ===")
# Assicuriamoci che sia compilato
subprocess.run(["make", "coupled"])

for N in MESH_SIZES:
    print(f"\nTesting mesh N = {N} (h = 1/{N})...")
    
    # 1. Aggiorna il file data
    update_config_file(DATA_FILE, N)
    
    # 2. Lancia la simulazione
    e_p, e_u = run_simulation()
    
    if e_p is not None:
        h = 1.0 / N
        h_values.append(h)
        errors_p.append(e_p)
        errors_u.append(e_u)
        print(f"  -> h={h:.4f} | Err P: {e_p:.2e} | Err U: {e_u:.2e}")

# --- CALCOLO ORDINE E GRAFICO ---
if len(h_values) > 1:
    # Calcola pendenze (ordine di convergenza) tra l'ultimo e il penultimo punto
    slope_p = np.log(errors_p[-2]/errors_p[-1]) / np.log(h_values[-2]/h_values[-1])
    slope_u = np.log(errors_u[-2]/errors_u[-1]) / np.log(h_values[-2]/h_values[-1])
    
    print("\n=== RISULTATI FINALI ===")
    print(f"Ordine di convergenza stimato (Pressione): {slope_p:.2f}")
    print(f"Ordine di convergenza stimato (Spostamento): {slope_u:.2f}")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.loglog(h_values, errors_p, '-o', label=f'Pressure Error (Slope ~{slope_p:.1f})')
    plt.loglog(h_values, errors_u, '-s', label=f'Displacement Error (Slope ~{slope_u:.1f})')
    
    # Linee di riferimento (pendenza 1 e 2)
    ref_h = np.array(h_values)
    plt.loglog(ref_h, ref_h * (errors_p[0]/ref_h[0]), '--', color='gray', alpha=0.5, label='Order 1 (Reference)')
    plt.loglog(ref_h, ref_h**2 * (errors_p[0]/ref_h[0]**2), ':', color='gray', alpha=0.5, label='Order 2 (Reference)')
    
    plt.xlabel('Mesh size h')
    plt.ylabel('L2 Error Norm')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.title('Convergence Analysis')
    plt.savefig('convergence_plot.png')
    print("Grafico salvato come 'convergence_plot.png'")
    # plt.show() # Decommenta se hai X11 o display grafico
else:
    print("Non abbastanza dati per il plot.")



# per problema disaccoppiato:
# print("\n=== RISULTATI ===")
    
# # Calcolo Slope Pressione
# if len(h_values) >= 2:
#     slope_p = np.polyfit(np.log(h_values), np.log(errors_p), 1)[0]
#     print(f"Ordine Pressione: {slope_p:.2f}")
    
# # Calcolo Slope Spostamento (Gestione ZeroError)
# if all(e < 1.0e-15 for e in errors_u):
#      print("Ordine Spostamento: PERFETTO (Errore Zero)")
#          # Per il plot, sostituiamo gli zeri con un numero piccolissimo per non rompere il loglog
#      err_u = [1e-16 for _ in errors_u] 
# elif len(h_values) >= 2:
#     slope_u = np.polyfit(np.log(h_values), np.log(errors_u), 1)[0]
#     print(f"Ordine Spostamento: {slope_u:.2f}")

# # Plot
# plt.figure()
# plt.loglog(h_values, errors_p, '-o', label=f'Pressure (Order {slope_p:.2f})')
    
# # Riferimento ordine 2
# ref = [errors_p[0] * (h/h_values[0])**2 for h in h_values]
# plt.loglog(h_values, ref, 'k--', alpha=0.5, label='Order 2 Ref')
    
# plt.xlabel('Mesh Size h')
# plt.ylabel('L2 Error')
# plt.grid(True, which="both", ls="-")
# plt.legend()
# plt.title('Convergence Analysis (Decoupled Steady State)')
# plt.savefig('convergence_plot.png')
# print("Grafico salvato: convergence_plot.png")