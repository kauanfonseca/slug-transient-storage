import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ARIAL = "Arial"
TITLE_FONT = Font(name=ARIAL, size=13, bold=True)
SECTION_FONT = Font(name=ARIAL, size=11, bold=True, color="FFFFFF")
SECTION_FILL = PatternFill("solid", fgColor="2F5496")
HEADER_FONT = Font(name=ARIAL, size=10, bold=True)
HEADER_FILL = PatternFill("solid", fgColor="D9E1F2")
NOTE_FONT = Font(name=ARIAL, size=10, italic=True, color="7F0000")
WARN_FILL = PatternFill("solid", fgColor="FFC7CE")
WARN_FONT = Font(name=ARIAL, size=11, bold=True, color="9C0006")
BODY_FONT = Font(name=ARIAL, size=10)
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

events = pd.read_csv("data_derived/events.csv")
btc = pd.read_csv("data_derived/btc_conservative.csv")
mtsm = pd.read_csv("data_derived/master_tsm.csv")

TSM_PIPELINE_EVENTS = set(mtsm["event_id"].unique()) | {
    e for e in events["event_id"] if events.loc[events.event_id == e, "has_logger"].iloc[0] or True
}
# An event is "used in the TSM pipeline" if it isn't flagged invalid AND has either a
# logger or probe conservative-tracer series -- matches fit_all_hydraulics()'s own filter.
def used_in_pipeline(row):
    if bool(row.get("flag_discharge_invalid", False)):
        return False
    return bool(row.get("has_logger", False)) or row.get("discharge_source") == "probe"

HYD_FIELDS = [
    ("stream", "Stream code", ""),
    ("date", "Campaign date", ""),
    ("station", "Station", ""),
    ("slug_label", "Slug label", ""),
    ("addition_datetime", "Tracer/nutrient addition datetime", "America/Sao_Paulo"),
    ("reach_length_m", "Reach length", "m"),
    ("reach_mean_width_m", "Reach mean width", "m"),
    ("discharge_Ls", "Discharge (adopted)", "L/s"),
    ("water_velocity_ms", "Water velocity (adopted)", "m/s"),
    ("stream_depth_m", "Stream depth (adopted)", "m"),
    ("discharge_source", "Discharge source used for TSM hydraulics", "logger/probe/manual"),
    ("discharge_method", "Discharge estimation method", ""),
    ("nacl_mass_g", "NaCl tracer mass injected", "g"),
    ("amonium_added_g", "NH4Cl mass added", "g"),
    ("phophate_added_g", "PO4 salt mass added", "g"),
    ("ph", "pH", ""),
    ("temp_c", "Water temperature", "degC"),
    ("ca_mgL_up", "Ca2+ upstream", "mg/L"),
    ("ca_mgL_down", "Ca2+ downstream", "mg/L"),
    ("background_spc_probe", "Background SpC (probe)", "uS/cm"),
    ("nh4_background_ugL", "NH4-N background", "ug/L"),
    ("srp_background_ugL", "SRP background", "ug/L"),
    ("has_logger", "Has continuous conductivity logger?", ""),
    ("has_nutrients", "Has nutrient grabs?", ""),
    ("flag_discharge_invalid", "FLAG: discharge invalid", ""),
    ("flag_mixing_incomplete", "FLAG: mixing incomplete", ""),
    ("flag_no_baseline", "FLAG: no baseline", ""),
    ("flag_weak_signal", "FLAG: weak signal", ""),
    ("flag_sparse_recession", "FLAG: sparse recession", ""),
    ("flag_discharge_alt", "FLAG: alternate discharge used", ""),
]

def style_header_row(ws, row, ncols, start_col=1):
    for c in range(start_col, start_col + ncols):
        cell = ws.cell(row=row, column=c)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.border = BORDER
        cell.alignment = Alignment(vertical="center")

def section_title(ws, row, col, text, span):
    ws.merge_cells(start_row=row, start_column=col, end_row=row, end_column=col + span - 1)
    cell = ws.cell(row=row, column=col, value=text)
    cell.font = SECTION_FONT
    cell.fill = SECTION_FILL
    cell.alignment = Alignment(horizontal="left", vertical="center")

def autosize(ws, col_widths):
    for col, width in col_widths.items():
        ws.column_dimensions[get_column_letter(col)].width = width

wb = Workbook()

# ---------------------------------------------------------------- README ---
readme = wb.active
readme.title = "README"
readme["A1"] = "TSM / nutrient uptake pipeline -- review workbook"
readme["A1"].font = TITLE_FONT
readme["A2"] = ("Companion, human-readable export of the pipeline's derived tables, generated "
                "2026-09-23. One tab per campaign (event_id): hydraulics/campaign metadata, "
                "the conservative-tracer (NaCl) series (logger or hand-probe grabs, whichever "
                "was actually used to calibrate that event's transient-storage hydraulics), and "
                "the nutrient grab series. This workbook is a READ-ONLY convenience view -- the "
                "canonical, version-controlled source data is in data_derived/events.csv, "
                "data_derived/btc_conservative.csv and data_derived/master_tsm.csv (joined on "
                "event_id), documented in README.md and data_derived/dictionary/*.csv. Edit "
                "those, not this file, if a correction is needed.")
readme["A2"].alignment = Alignment(wrap_text=True)
readme.merge_cells("A2:G2")
readme.row_dimensions[2].height = 60
for c in range(1, 8):
    readme.column_dimensions[get_column_letter(c)].width = 16
readme.column_dimensions["A"].width = 24
readme.column_dimensions["B"].width = 24

idx_row = 4
headers = ["event_id", "stream", "date", "used in TSM pipeline?", "conservative-tracer source",
           "has nutrient grabs?", "note"]
for j, h in enumerate(headers, start=1):
    readme.cell(row=idx_row, column=j, value=h)
style_header_row(readme, idx_row, len(headers))

ev_sorted = events.sort_values(["stream", "date"]).reset_index(drop=True)
r = idx_row + 1
for _, row in ev_sorted.iterrows():
    eid = row["event_id"]
    used = used_in_pipeline(row)
    src = row.get("discharge_source")
    src = src if isinstance(src, str) else ("logger" if row.get("has_logger") else "none")
    has_nut = eid in set(mtsm.loc[mtsm["solute"].notna(), "event_id"])
    note = ""
    if bool(row.get("flag_discharge_invalid", False)):
        note = "Excluded: flag_discharge_invalid"
    elif not has_nut:
        note = "No nutrient addition at this station (NaCl/hydraulics only)"
    elif eid == "SR_20231011_single":
        note = "No logger; hydraulics borrowed from SR_20231009_downstream (see handoff.md)"
    readme.cell(row=r, column=1, value=eid)
    readme.cell(row=r, column=2, value=row["stream"])
    readme.cell(row=r, column=3, value=str(row["date"]))
    readme.cell(row=r, column=4, value="yes" if used else "no")
    readme.cell(row=r, column=5, value=src)
    readme.cell(row=r, column=6, value="yes" if has_nut else "no")
    readme.cell(row=r, column=7, value=note)
    for c in range(1, 8):
        readme.cell(row=r, column=c).font = BODY_FONT
        readme.cell(row=r, column=c).border = BORDER
    if not used:
        for c in range(1, 8):
            readme.cell(row=r, column=c).fill = WARN_FILL
    r += 1
readme.freeze_panes = "A5"
readme.column_dimensions["G"].width = 55

# ------------------------------------------------------- per-event sheets ---
def sheet_name_for(eid):
    return eid[:31]

for _, ev in ev_sorted.iterrows():
    eid = ev["event_id"]
    ws = wb.create_sheet(sheet_name_for(eid))
    title = f"{eid}  --  {ev['stream']} campaign on {ev['date']}"
    ws["A1"] = title
    ws["A1"].font = TITLE_FONT
    ws.merge_cells("A1:F1")

    row_cursor = 2
    if bool(ev.get("flag_discharge_invalid", False)):
        ws.cell(row=row_cursor, column=1,
                value="FLAGGED INVALID (flag_discharge_invalid = TRUE) -- not used in the TSM pipeline. Data kept below for reference only.")
        ws.cell(row=row_cursor, column=1).font = WARN_FONT
        ws.cell(row=row_cursor, column=1).fill = WARN_FILL
        ws.merge_cells(start_row=row_cursor, start_column=1, end_row=row_cursor, end_column=6)
        row_cursor += 1
    row_cursor += 1

    # --- Block A: hydraulics / campaign metadata (columns A-B) ---
    hyd_start = row_cursor
    section_title(ws, hyd_start, 1, "CAMPAIGN / HYDRAULICS METADATA", 2)
    r = hyd_start + 1
    ws.cell(row=r, column=1, value="Field"); ws.cell(row=r, column=2, value="Value")
    style_header_row(ws, r, 2)
    r += 1
    field_start = r
    for col, label, unit in HYD_FIELDS:
        val = ev.get(col, "")
        if pd.isna(val):
            val = ""
        label_txt = f"{label} ({unit})" if unit else label
        ws.cell(row=r, column=1, value=label_txt).font = BODY_FONT
        ws.cell(row=r, column=2, value=val).font = BODY_FONT
        ws.cell(row=r, column=1).border = BORDER
        ws.cell(row=r, column=2).border = BORDER
        r += 1
    hyd_end = r - 1

    # --- Block B: conservative-tracer (NaCl) series (columns D onward) ---
    trac_col = 4  # column D
    has_logger = bool(ev.get("has_logger", False))
    is_probe_source = ev.get("discharge_source") == "probe"
    logger_sub = btc[(btc.event_id == eid) & (btc.time_since_release_s >= 0)].sort_values("time_since_release_s")
    probe_sub = mtsm[(mtsm.event_id == eid) & (mtsm.nacl_mgL_grab.notna()) & (mtsm.time_since_release_s >= 0)] \
        .drop_duplicates(subset=["time_since_release_s", "nacl_mgL_grab"]).sort_values("time_since_release_s")

    if len(logger_sub) > 0:
        section_title(ws, hyd_start, trac_col, "NaCl -- CONTINUOUS LOGGER SERIES" +
                      ("" if not is_probe_source else " (recorded, but PROBE grabs below were used for TSM hydraulics)"), 4)
        r2 = hyd_start + 1
        cols = ["time_since_release_s", "spc_uscm", "spc_corr_uscm", "nacl_mgL"]
        col_labels = ["time_since_release_s (s)", "spc_uscm (uS/cm)", "spc_corr_uscm (uS/cm)", "nacl_mgL (mg/L)"]
        for j, lab in enumerate(col_labels):
            ws.cell(row=r2, column=trac_col + j, value=lab)
        style_header_row(ws, r2, len(cols), start_col=trac_col)
        r2 += 1
        for _, lr in logger_sub.iterrows():
            for j, c in enumerate(cols):
                v = lr[c]
                ws.cell(row=r2, column=trac_col + j, value=(None if pd.isna(v) else v)).font = BODY_FONT
            r2 += 1
        logger_end_row = r2 - 1
    elif len(probe_sub) > 0:
        section_title(ws, hyd_start, trac_col, "NaCl -- HAND-PROBE GRAB SERIES (no continuous logger for this event)", 2)
        r2 = hyd_start + 1
        cols = ["time_since_release_s", "nacl_mgL_grab"]
        col_labels = ["time_since_release_s (s)", "nacl_mgL_grab (mg/L)"]
        for j, lab in enumerate(col_labels):
            ws.cell(row=r2, column=trac_col + j, value=lab)
        style_header_row(ws, r2, len(cols), start_col=trac_col)
        r2 += 1
        for _, lr in probe_sub.iterrows():
            for j, c in enumerate(cols):
                v = lr[c]
                ws.cell(row=r2, column=trac_col + j, value=(None if pd.isna(v) else v)).font = BODY_FONT
            r2 += 1
        logger_end_row = r2 - 1
    else:
        ws.cell(row=hyd_start, column=trac_col, value="No NaCl logger or probe series recorded for this event.").font = NOTE_FONT
        logger_end_row = hyd_start

    # If discharge_source == probe but a logger series ALSO exists, still show probe grabs
    # separately (they are what Stage 1 actually used) in a second block beneath the logger one.
    if is_probe_source and len(logger_sub) > 0 and len(probe_sub) > 0:
        pr_start = logger_end_row + 2
        section_title(ws, pr_start, trac_col, "NaCl -- HAND-PROBE GRABS (this is what TSM hydraulics were actually fit to)", 2)
        r2 = pr_start + 1
        cols = ["time_since_release_s", "nacl_mgL_grab"]
        col_labels = ["time_since_release_s (s)", "nacl_mgL_grab (mg/L)"]
        for j, lab in enumerate(col_labels):
            ws.cell(row=r2, column=trac_col + j, value=lab)
        style_header_row(ws, r2, len(cols), start_col=trac_col)
        r2 += 1
        for _, lr in probe_sub.iterrows():
            for j, c in enumerate(cols):
                v = lr[c]
                ws.cell(row=r2, column=trac_col + j, value=(None if pd.isna(v) else v)).font = BODY_FONT
            r2 += 1

    # --- Block C: nutrient grabs (columns I onward) ---
    nut_col = 9  # column I
    nut_sub = mtsm[(mtsm.event_id == eid) & (mtsm.solute.notna())].sort_values(["solute", "time_since_release_s"])
    if len(nut_sub) > 0:
        section_title(ws, hyd_start, nut_col, "NUTRIENT GRABS", 6)
        r3 = hyd_start + 1
        cols = ["solute", "time_since_release_s", "conc_ugL", "background_ugL", "conc_corr_ugL", "conc_corr_coprec_ugL"]
        col_labels = ["solute", "time_since_release_s (s)", "conc_ugL (raw, ug/L)", "background_ugL (ug/L)",
                      "conc_corr_ugL (bg-corrected, ug/L)", "conc_corr_coprec_ugL (SRP co-precip added back, ug/L)"]
        for j, lab in enumerate(col_labels):
            ws.cell(row=r3, column=nut_col + j, value=lab)
        style_header_row(ws, r3, len(cols), start_col=nut_col)
        r3 += 1
        for _, nr in nut_sub.iterrows():
            for j, c in enumerate(cols):
                v = nr[c]
                ws.cell(row=r3, column=nut_col + j, value=(None if pd.isna(v) else v)).font = BODY_FONT
            r3 += 1
    else:
        ws.cell(row=hyd_start, column=nut_col,
                value="No nutrient addition/grabs recorded for this event (NaCl/hydraulics-only campaign).").font = NOTE_FONT

    ws.freeze_panes = "A" + str(field_start)
    widths = {1: 40, 2: 22, 3: 2, 4: 20, 5: 14, 6: 15, 7: 12, 8: 2,
              9: 8, 10: 20, 11: 18, 12: 16, 13: 20, 14: 30}
    autosize(ws, widths)

wb.save("data_derived/tsm_review_workbook.xlsx")
print("saved: data_derived/tsm_review_workbook.xlsx, sheets:", wb.sheetnames)
