import pandas as pd
import numpy as np

events = pd.read_csv("data_derived/events.csv")

events["point_id"] = np.where(events["station"].isna(), events["stream"],
                               events["stream"] + "_" + events["station"].fillna(""))

METRIC_COLS = ["discharge_Ls", "water_velocity_ms", "stream_depth_m",
               "reach_length_m", "reach_mean_width_m"]

rows = []
for point_id, grp in events.groupby("point_id"):
    stream = grp["stream"].iloc[0]
    station = grp["station"].iloc[0]
    excluded = grp[grp["flag_discharge_invalid"] == True]
    valid = grp[grp["flag_discharge_invalid"] != True]
    for source in ["logger", "probe", "manual"]:
        sub = valid[valid["discharge_source"] == source]
        if len(sub) == 0:
            continue
        row = {
            "point_id": point_id, "stream": stream, "station": station,
            "discharge_source": source, "n_campaigns": len(sub),
            "event_ids": ";".join(sorted(sub["event_id"])),
        }
        for c in METRIC_COLS:
            vals = sub[c].dropna()
            row[f"mean_{c}"] = vals.mean() if len(vals) else np.nan
            row[f"sd_{c}"] = vals.std(ddof=1) if len(vals) > 1 else np.nan
            row[f"min_{c}"] = vals.min() if len(vals) else np.nan
            row[f"max_{c}"] = vals.max() if len(vals) else np.nan
        row["excluded_events_note"] = (
            f"{len(excluded)} event(s) at this point excluded (flag_discharge_invalid): "
            + ";".join(sorted(excluded["event_id"]))
        ) if len(excluded) else ""
        rows.append(row)
    # a point where EVERY event was excluded still gets one row, all-NA, so it's not silently dropped
    if len(valid) == 0 and len(excluded) > 0:
        row = {"point_id": point_id, "stream": stream, "station": station,
               "discharge_source": "NONE_VALID", "n_campaigns": 0, "event_ids": ""}
        for c in METRIC_COLS:
            row[f"mean_{c}"] = np.nan; row[f"sd_{c}"] = np.nan
            row[f"min_{c}"] = np.nan; row[f"max_{c}"] = np.nan
        row["excluded_events_note"] = (
            f"ALL {len(excluded)} event(s) at this point excluded (flag_discharge_invalid): "
            + ";".join(sorted(excluded["event_id"]))
        )
        rows.append(row)

out = pd.DataFrame(rows)
col_order = (["point_id", "stream", "station", "discharge_source", "n_campaigns", "event_ids"] +
             [f"{stat}_{c}" for c in METRIC_COLS for stat in ["mean", "sd", "min", "max"]] +
             ["excluded_events_note"])
out = out[col_order].sort_values(["point_id", "discharge_source"]).reset_index(drop=True)
out.to_csv("data_derived/hydraulics_by_point.csv", index=False)
print(out.to_string())
