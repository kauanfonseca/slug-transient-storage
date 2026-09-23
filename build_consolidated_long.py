import pandas as pd
import numpy as np

events = pd.read_csv("data_derived/events.csv")
btc = pd.read_csv("data_derived/btc_conservative.csv")
mtsm = pd.read_csv("data_derived/master_tsm.csv")

HYD_COLS = [
    "reach_length_m", "reach_mean_width_m", "discharge_Ls", "water_velocity_ms",
    "stream_depth_m", "discharge_source", "discharge_method", "nacl_mass_g",
    "amonium_added_g", "phophate_added_g", "ph", "temp_c",
    "has_logger", "has_nutrients",
    "flag_discharge_invalid", "flag_mixing_incomplete", "flag_no_baseline", "flag_weak_signal",
]
ID_COLS = ["event_id", "stream", "date", "station", "slug_label"]

hyd_lookup = events.set_index("event_id")[HYD_COLS]
id_lookup = events.set_index("event_id")[["stream", "date", "station", "slug_label"]]

DATA_COLS = ["spc_uscm", "spc_corr_uscm", "nacl_mgL", "nacl_mgL_grab",
             "conc_ugL", "background_ugL", "conc_corr_ugL", "conc_corr_coprec_ugL"]

def blank_frame(n):
    return pd.DataFrame({c: [np.nan] * n for c in DATA_COLS})

rows = []

# --- 1. logger rows -----------------------------------------------------
logger_sub = btc[btc.time_since_release_s >= 0].copy()
if len(logger_sub):
    d = blank_frame(len(logger_sub))
    d["spc_uscm"] = logger_sub["spc_uscm"].values
    d["spc_corr_uscm"] = logger_sub["spc_corr_uscm"].values
    d["nacl_mgL"] = logger_sub["nacl_mgL"].values
    d["event_id"] = logger_sub["event_id"].values
    d["time_since_release_s"] = logger_sub["time_since_release_s"].values
    d["source"] = "logger"
    d["solute"] = np.nan
    rows.append(d)

# --- 2. hand-probe NaCl grab rows (from master_tsm, deduplicated) -------
probe_sub = mtsm[mtsm.nacl_mgL_grab.notna() & (mtsm.time_since_release_s >= 0)] \
    .drop_duplicates(subset=["event_id", "time_since_release_s", "nacl_mgL_grab"]).copy()
if len(probe_sub):
    d = blank_frame(len(probe_sub))
    d["nacl_mgL_grab"] = probe_sub["nacl_mgL_grab"].values
    d["event_id"] = probe_sub["event_id"].values
    d["time_since_release_s"] = probe_sub["time_since_release_s"].values
    d["source"] = "nacl_probe_grab"
    d["solute"] = np.nan
    rows.append(d)

# --- 3. nutrient grab rows -----------------------------------------------
nut_sub = mtsm[mtsm.solute.notna()].copy()
if len(nut_sub):
    d = blank_frame(len(nut_sub))
    for c in ["conc_ugL", "background_ugL", "conc_corr_ugL", "conc_corr_coprec_ugL"]:
        d[c] = nut_sub[c].values
    d["event_id"] = nut_sub["event_id"].values
    d["time_since_release_s"] = nut_sub["time_since_release_s"].values
    d["source"] = "nutrient_grab"
    d["solute"] = nut_sub["solute"].values
    rows.append(d)

out = pd.concat(rows, ignore_index=True)

# attach event identity + hydraulics columns
out = out.join(id_lookup, on="event_id")
out = out.join(hyd_lookup, on="event_id")

col_order = ID_COLS + ["time_since_release_s", "source", "solute"] + DATA_COLS + HYD_COLS
out = out[col_order].sort_values(["event_id", "source", "time_since_release_s"]).reset_index(drop=True)

out.to_csv("data_derived/tsm_consolidated_long.csv", index=False)
print("rows:", len(out), " events:", out.event_id.nunique())
print(out["source"].value_counts())
