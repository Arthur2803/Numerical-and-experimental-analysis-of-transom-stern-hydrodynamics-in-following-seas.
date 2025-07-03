#Denne filen tar inn en tidsserie med vannhøyder på akterspeilet og returnerer en graf som viser basedrag tidsserie

import numpy as np
import capytaine as cpt
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
from apread import APReader
import pandas as pd
from pathlib import Path

_cached_wetted_area = None

widths = np.array([0, 0.04642, 0.0599, 0.0687, 0.07933, 0.08625, 0.09375, 0.09859, 0.10401, 0.09898, 0.093752])  # m
heights = np.array([0, 0.001177, 0.00392, 0.00865, 0.02201, 0.03408, 0.0466, 0.0591, 0.07415, 0.0882, 0.10254]) 

def extract_variable_from_filename(file_name, variable):
    parts = file_name.split('-', 3)
    if len(parts) < 4:
        raise ValueError("Invalid file name format.")
    date = parts[0]
    velocity = parts[1].replace("_", ".")
    period = parts[2].replace("_", ".")
    steepness_raw, _, test_part = parts[3].partition('#')
    test = test_part.replace(".bin", "")
    if '(' in steepness_raw and ')' in steepness_raw:
        i = steepness_raw.find('(')
        j = steepness_raw.find(')', i)
        steepness_val = steepness_raw[:i]
        info_val = steepness_raw[i+1:j]
    elif '-' in steepness_raw:
        steepness_val, info_val = steepness_raw.split('-', 1)
    else:
        steepness_val = steepness_raw
        info_val = None
    var = variable.lower()
    if var == "date":
        return date
    elif var == "velocity":
        try:
            return float(velocity)
        except ValueError:
            return velocity
    elif var == "period":
        try:
            return float(period)
        except ValueError:
            return period
    elif var == "steepness":
        try:
            return float(steepness_val)
        except ValueError:
            return steepness_val
    elif var == "info":
        return info_val
    elif var in ("test", "testnumber"):
        try:
            return int(test)
        except ValueError:
            return test
    else:
        raise ValueError("Invalid variable requested.")

def get_speedaligned_channel(file_list, plot_channel="FX_total", align_channel="Speed", root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):
    # hvis en fil putter i liste.
    if isinstance(file_list, str):
        file_list = [file_list]
    for file_name in file_list:
        matches = list(Path(root_folder).rglob(file_name))
        if not matches:
            print(f"File not found: {file_name}")
            continue
        full_path = matches[0]
        try:
            parts = file_name.split("-")
            velocity_str = parts[1]  # eks "0_6"
            target_speed = -float(velocity_str.replace("_", "."))
        except Exception as e:
            print(f"Error extracting target speed from {file_name}: {e}")
            continue
        df = bin_to_dataframe(str(full_path))
        condition = df[align_channel].abs() >= abs(target_speed)
        if not condition.any():
            print(f"Target speed {target_speed} not reached in {file_name}")
            continue
        start_time = condition.idxmax()

        end_time = condition[::-1].idxmax()
        df_aligned = df.loc[start_time:end_time].copy()
        df_aligned.index -= start_time
        return df_aligned[plot_channel]
    return None

def bin_to_dataframe(bin_file):
    reader = APReader(bin_file)
    data = {str(ch).split('"')[1]: ch.data for ch in reader.Channels}
    time_key = [key for key in data if "Time" in key][0]
    df = pd.DataFrame(data)
    df.set_index(time_key, inplace=True)
    return df


def wetted_area(mesh_name="2112open.stl", rotation_center=None):
    global _cached_wetted_area
    if _cached_wetted_area is not None:
        return _cached_wetted_area
   
    mesh = cpt.load_mesh(mesh_name,file_format="stl")

    rc = rotation_center or mesh.center_of_buoyancy

    body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=rc),
        center_of_mass=rc  # use COB if CG unknown
    )
    # 4) compute hydrostatics and return wet surface area
    hydro = body.compute_hydrostatics(rho=1000, g=9.81)
    print(f"Hydrostatics: {hydro}")
    _cached_wetted_area = float(hydro["wet_surface_area"])
    return _cached_wetted_area


def calc_Friction_coeff_speedaligned(file_name, dt=0.005, cutoff_hz=0.05, root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):
    #Takes the drag force and returns C_F timeseries
    #measured_drag = get_speedaligned_channel(file_name, "FX_total", root_folder=root_folder)
   
    #nyquist = 0.5 / dt
    #norm_cutoff = cutoff_hz / nyquist
    #b, a = butter(N=3, Wn=norm_cutoff, btype="low")
    #measured_drag_filtered = filtfilt(b, a, measured_drag)
    #wetted_area = wetted_area("2112open.stl")
    boat_speed = extract_variable_from_filename(file_name, "velocity")

    
    #CF = measured_drag_filtered / (0.5 * 1000 * (boat_speed**2) *wetted_area )
    CF = 0.075/(np.log10(boat_speed*1*1/1e-6)-2)**2
    print(f"CF: {CF}")
    
    return CF


def calc_Dynamic_Basedrag_coeff(file_name, dt=0.005, cutoff_hz=20, root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):
    water1_unfiltered = get_speedaligned_channel(file_name, "Rel_WP1", root_folder=root_folder)
    water2_unfiltered = get_speedaligned_channel(file_name, "Rel_WP2", root_folder=root_folder)
    #
    nyquist = 0.5 / dt
    norm_cutoff = cutoff_hz / nyquist
    b, a = butter(N=3, Wn=norm_cutoff, btype="low")
    water1 = filtfilt(b, a, water1_unfiltered)
    water2 = filtfilt(b, a, water2_unfiltered)

    
    # Average water height
    water_height = (water1 + water2) / 2 + 0.01

    # Hull shape (static)
    widths = np.array([0, 0.04642, 0.0599, 0.0687, 0.07933, 0.08625, 0.09375, 0.09859, 0.10401, 0.09898, 0.093752])  # m
    heights = np.array([0, 0.001177, 0.00392, 0.00865, 0.02201, 0.03408, 0.0466, 0.0591, 0.07415, 0.0882, 0.10254])  

    # Calculate wetted area of transom time series
    wetted_area_transom = []
    for h in water_height:
        h = np.clip(h,0, heights[-1])  # Clip within hull limits
        wetted_widths = []
        wetted_heights = []
        for i in range(len(heights) - 1):
            if heights[i+1] >= h:
                w_interp = widths[i] + (widths[i+1] - widths[i]) * (h - heights[i]) / (heights[i+1] - heights[i])
                wetted_widths += [widths[i], w_interp]
                wetted_heights += [heights[i], h]
                break
            wetted_widths += [widths[i]]
            wetted_heights += [heights[i]]
        wetted_widths = [0] + wetted_widths
        wetted_heights = [0] + wetted_heights
        area = 2 * np.trapezoid(wetted_widths, wetted_heights)
        wetted_area_transom.append(area)
    wetted_area_transom = np.array(wetted_area_transom)
    
    
    #CAPYTAINE FINNE VÅTT AREAL
    wett_area = wetted_area("2112open.stl")
    wetted_area_minus_transom = wett_area - wetted_area_transom

    #CF = calc_Friction_coeff_speedaligned(file_name, dt=dt, cutoff_hz=cutoff_hz, root_folder=root_folder)
    
    CF = 0.002
    
    C_BD = 0.029*(((wetted_area_transom/wetted_area_minus_transom)**3)/CF)**0.5
    mean_C_BD = np.mean(C_BD)
    print(f"Mean C_BD: {mean_C_BD}")
    
    # Take second half
    mid_idx = len(C_BD) // 2
    second_half = C_BD[mid_idx:]

    # Return mean
    return second_half


def plot_CF_and_CBD(file_name, dt=0.005, cutoff_hz=20, root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):
    #CF = calc_Friction_coeff_speedaligned(file_name, dt=dt, cutoff_hz=cutoff_hz, root_folder=root_folder)
    CF = 0.002
    C_BD = calc_Dynamic_Basedrag_coeff(file_name, dt=dt, cutoff_hz=cutoff_hz, root_folder=root_folder)
    
    time = np.arange(len(C_BD)) * dt

    CF_const = CF[0] if hasattr(CF, "__len__") else CF
    
    plt.figure(figsize=(12, 6))
    plt.plot(time, CF_const, label='C_F', color='blue')
    plt.plot(time, C_BD, label='C_BD', color='orange')
    plt.xlabel('Time [s]')
    plt.ylabel('Coefficient')
    plt.title('C_F and C_BD over Time')
    plt.legend()
    plt.grid()
    plt.show()
    
def plot_totalDrag_and_BaseDrag(file_name, dt=0.05, cutoff_hz=20, root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):
    totalDragUnfiltered = get_speedaligned_channel(file_name, "FX_total", root_folder=root_folder)
    C_BD = calc_Dynamic_Basedrag_coeff(file_name, dt=dt, cutoff_hz=cutoff_hz, root_folder=root_folder)
    
    nyq = 0.5 / dt
    b, a = butter(3, cutoff_hz/nyq, btype="low")
    totalDrag = filtfilt(b, a, totalDragUnfiltered)
    
    mid_idx = len(totalDrag) // 2
    totalDrag_second_half = totalDrag[mid_idx:]
    
    water1_unfiltered = get_speedaligned_channel(file_name, "Rel_WP1", root_folder=root_folder)
    water2_unfiltered = get_speedaligned_channel(file_name, "Rel_WP2", root_folder=root_folder)
    #
    nyquist = 0.5 / dt
    norm_cutoff = cutoff_hz / nyquist
    b, a = butter(N=3, Wn=norm_cutoff, btype="low")
    water1 = filtfilt(b, a, water1_unfiltered)
    water2 = filtfilt(b, a, water2_unfiltered)

    
    # Average water height
    water_height = (water1 + water2) / 2 +0.01

    # Hull shape (static)
    widths = np.array([0, 0.04642, 0.0599, 0.0667, 0.07633, 0.08425, 0.09175, 0.09859, 0.10401, 0.09898, 0.093752])
    heights = np.array([0, 0.001177, 0.00392, 0.00865, 0.02201, 0.03408, 0.0466, 0.0591, 0.07415, 0.0882, 0.10254])   

    # Calculate wetted area of transom time series
    wetted_area_transom = []
    for h in water_height:
        wetted_widths = []
        wetted_heights = []
        for i in range(len(heights) - 1):
            if heights[i+1] >= h:
                w_interp = widths[i] + (widths[i+1] - widths[i]) * (h - heights[i]) / (heights[i+1] - heights[i])
                wetted_widths += [widths[i], w_interp]
                wetted_heights += [heights[i], h]
                break
            wetted_widths += [widths[i]]
            wetted_heights += [heights[i]]
        wetted_widths = [0] + wetted_widths
        wetted_heights = [0] + wetted_heights
        area = 2 * np.trapezoid(wetted_widths, wetted_heights)
        wetted_area_transom.append(area)
    wetted_area_transom = np.array(wetted_area_transom)
    
    mid_idx = len(wetted_area_transom) // 2
    wetted_second_half = wetted_area_transom[mid_idx:]
    
    
    D_B = C_BD * 0.5 * 1000 * extract_variable_from_filename(file_name, "velocity")**2 * wetted_second_half
    
    time = np.arange(len(C_BD)) * dt

    plt.figure(figsize=(12, 6))
    plt.plot(time, totalDrag_second_half, label='Total Drag', color='blue')
    plt.plot(time, D_B, label='D_B', color='orange')
    plt.xlabel('Time [s]')
    plt.ylabel('Coefficient')
    plt.title('TotalDrag and Base drag over Time')
    plt.legend()
    plt.grid()
    plt.show()
    

def plot_Re_transition(
    thresholds=[1e4, 1e5, 1e6], L=1.0, nu=1e-6, speeds=None, num_x=2000):
    if speeds is None:
        speeds = np.linspace(0, 1.5, 1500)
    x = np.linspace(0, L, num_x)

    plt.figure(figsize=(6,4))
    for thr in thresholds:
        x_trans = []
        for V in speeds:
            Re = V * x / nu
            idx = np.argmax(Re > thr)
            x_t = x[idx] if Re[idx] > thr else np.nan
            x_trans.append(x_t)
        plt.plot(speeds, x_trans, label=f'Re > {thr:.0e}')

    plt.xlabel("Boat speed V [m/s]")
    plt.ylabel("Transition location xₜ from bow [m]")
    plt.title("Laminar–Turbulent Transition for Various Re Thresholds")
    plt.legend()
    plt.grid(True)
    plt.show()


#### ????? #####
def prohaska_method(R, V, S, L, rho=1000.0, nu=1e-6, n=5):

    Ct = R / (0.5 * rho * V**2 * S)
    Re = V * L / nu
    Cf = 0.075 / (np.log10(Re) - 2)**2
    
    # compute Froude number
    Fn = V / np.sqrt(9.81 * L)
    
    # build Prohaska variables
    X = Fn**n / Cf
    Y = Ct / Cf
    
    # linear regression Y = a + b X
    a, b = np.polyfit(X, Y, 1)
    k = a - 1.0  # form-factor
    m = b        # wave-resistance pseudocoeff
    
    # Plot Prohaska diagram
    plt.figure(figsize=(6,4))
    plt.plot(X, Y, 'o', label='data')
    plt.plot(X, a + b*X, '-', label=f'fit: intercept={a:.3f}, slope={b:.3f}')
    plt.xlabel(r'$F_n^n / C_F$')
    plt.ylabel(r'$C_T / C_F$')
    plt.title("Prohaska's Method Plot")
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return k, m, Fn, X, Y


def calc_R_hekk(file_name, dt=0.005, cutoff_hz=20, root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater", rho=1000.0, g=9.81):

    w1_raw = get_speedaligned_channel(file_name, "Rel_WP1", root_folder=root_folder)
    w2_raw = get_speedaligned_channel(file_name, "Rel_WP2", root_folder=root_folder)
    nyq = 0.5 / dt
    b, a = butter(3, cutoff_hz/nyq, btype="low")
    w1 = filtfilt(b, a, w1_raw)
    w2 = filtfilt(b, a, w2_raw)
    
    
    T = (w1 + w2) / 2         # transom-dybde over tid

    # --- 2) statisk skrogprofil (hull shape) ---
    widths = np.array([0, 0.04642, 0.0599, 0.0687, 0.07933, 0.08625, 0.09375, 0.09859, 0.10401, 0.09898, 0.093752])  # m
    heights = np.array([0, 0.001177, 0.00392, 0.00865, 0.02201, 0.03408, 0.0466, 0.0591, 0.07415, 0.0882, 0.10254])

    # --- 3) for hver T[i], regn ∫₀ᵀ z·B(z) dz via trapz ---
    R_hekk = np.zeros_like(T)
    for i, h in enumerate(np.clip(T, 0, heights[-1])):
        if h == 0:
            continue
        j = np.searchsorted(heights, h) - 1
        # bygg z- og B-vektor opp til h
        z = np.concatenate(( [0], heights[1:j+1], [h] ))
        B_tail = widths[j] + (widths[j+1]-widths[j]) * (h-heights[j])/(heights[j+1]-heights[j])
        B = np.concatenate(( [0], widths[1:j+1], [B_tail] ))
        R_hekk[i] = rho * g * np.trapezoid(z * B, z)
    return R_hekk



def plot_transom_wetted_area(file_name,
                             dt=0.005,
                             cutoff_hz=20,
                             root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):

    # 1) les vannhøyder
    w1 = get_speedaligned_channel(file_name, "Rel_WP1", root_folder=root_folder)
    w2 = get_speedaligned_channel(file_name, "Rel_WP2", root_folder=root_folder)

    # 2) low-pass filter
    nyq = 0.5 / dt
    b, a = butter(3, cutoff_hz/nyq, btype="low")
    h1 = filtfilt(b, a, w1)
    h2 = filtfilt(b, a, w2)

    # 3) gjennomsnittlig dybde over transom
    water_height = (h1 + h2) / 2 + 0.01

    # 4) statisk profil z vs. B(z)
    widths = np.array([0, 0.04642, 0.0599, 0.0687, 0.07933, 0.08625, 0.09375, 0.09859, 0.10401, 0.09898, 0.093752])  # m
    heights = np.array([0, 0.001177, 0.00392, 0.00865, 0.02201, 0.03408, 0.0466, 0.0591, 0.07415, 0.0882, 0.10254])

    # 5) beregn transom-våtområde over tid
    wetted_transom = []
    for h in np.clip(water_height, 0, heights[-1]):
        j = np.searchsorted(heights, h) - 1
        z = np.concatenate(([0], heights[1:j+1], [h]))
        B_tail = widths[j] + (widths[j+1] - widths[j]) * (h - heights[j]) / (heights[j+1] - heights[j])
        B = np.concatenate(([0], widths[1:j+1], [B_tail]))
        A_t = 2 * np.trapz(B, z)
        wetted_transom.append(A_t)
    wetted_transom = np.array(wetted_transom)

    # 6) plott
    t = np.arange(len(wetted_transom)) * dt
    plt.figure(figsize=(10, 5))
    plt.plot(t, wetted_transom, label='A_transom')
    plt.xlabel('Tid [s]')
    plt.ylabel('Transom våt flate [m²]')
    plt.title('Transom våt flate over tid')
    plt.grid(True)
    plt.legend()
    plt.show()



#filename= "0901-1_1-0_75-40#1.bin"
#filename = "1601-1_1-0_8-40(1kgAkter)#1.bin"
#filename = "1401-0_7-0_65-40#4.bin"
filename = "1001-1_1-0_675-40#1.bin"
#R_ts = calc_R_hekk(filename, dt=0.005, cutoff_hz=20)
#t = np.arange(len(R_ts)) * 0.005
#plt.plot(t, R_ts)
#plt.xlabel("Tid [s]")
#plt.ylabel(r"$R_{\rm hekk}$ [N]")
#plt.show()

def plot_CBD_dynamisk(file_name,
             dt=0.005,
             cutoff_hz=20,
             root_folder=r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater"):

    # beregn C_BD-tidsserie (andre halvdel)
    C_BD = calc_Dynamic_Basedrag_coeff(
        file_name,
        dt=dt,
        cutoff_hz=cutoff_hz,
        root_folder=root_folder
    )

    # tid for hver C_BD-verdi
    t = np.arange(len(C_BD)) * dt

    # plott
    plt.figure(figsize=(10, 5))
    plt.plot(t, C_BD, label='C_BD', linewidth=1)
    plt.xlabel('Tid [s]')
    plt.ylabel('Dynamic basedragkoeff. $C_{BD}$')
    plt.title('Tidsserie av dynamic base-drag-koeffisient')
    plt.grid(True)
    plt.legend()
    plt.show()


SENSOR_OFFSET = 0.01    # [m], sensorkalibrering mot vannlinje
CF_CONST      = 0.002 

def compute_transom_area_at_height(h_raw, offset=SENSOR_OFFSET):
    """
    Beregn transom-våtområde [m²] for rå vannhøyde h_raw [m] fra sensor.
    """
    h = np.clip(h_raw + offset, 0.0, heights[-1])
    wetted_w = [0.0]
    wetted_z = [0.0]
    for i in range(len(heights)-1):
        if heights[i+1] >= h:
            w0, w1 = widths[i], widths[i+1]
            z0, z1 = heights[i], heights[i+1]
            w_interp = w0 + (w1 - w0)*(h - z0)/(z1 - z0)
            wetted_w += [w0, w_interp]
            wetted_z += [z0, h]
            break
        wetted_w.append(widths[i])
        wetted_z.append(heights[i])
    area = 2.0 * np.trapz(wetted_w, wetted_z)
    return max(area, 0.0)

def CBD_at_given_height( h_raw):
    """
    Beregn base-drag koeff C_BD for én rå vannhøyde h_raw.
    """
    #CF = calc_Friction_coeff_speedaligned(file_name)
    CF = 0.002
    A_trans = compute_transom_area_at_height(h_raw)
    A_total = wetted_area("2112open.stl")
    A_hull  = A_total - A_trans
    ratio   = A_trans/A_hull if A_hull>0 else np.nan
    return 0.029 * np.sqrt((ratio**3)/CF)

def plot_CBD_vs_heights(file_name, num_points=100):
    """
    Plott C_BD som funksjon av rå vannhøyder fra 0 til max-sensorverdi.
    """
    h_vals = np.linspace(0.0, heights[-1]-SENSOR_OFFSET, num_points)
    cbd_vals = [CBD_at_given_height(h) for h in h_vals]
    plt.figure(figsize=(8,5))
    plt.plot(h_vals, cbd_vals, '-', label=r'$C_{BD}(h)$')
    plt.xlabel('Rå vannhøyde h [m]')
    plt.ylabel('Base-drag koeff. $C_{BD}$')
    plt.title('C_{BD} som funksjon av vannhøyde')
    plt.grid(True)
    plt.legend()
    plt.show()


# EKSEMPEL PÅ BRUk
if __name__ == "__main__":
    filename = "1001-1_1-0_675-40#1.bin"
    # Tidsserieplott:           # eksisterende tidsserie-plott

    # Statisk CBD vs. høyde:
    #plot_CBD_vs_heights(filename)

    dO=compute_transom_area_at_height(0.02, offset=0.01)
    print(f"Transom area at height 0.02 m: {dO:.4f} m²")
    
# Sensorkalibrerings-offset og friksjonskoeffisient


