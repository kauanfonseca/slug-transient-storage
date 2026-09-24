---
editor_options: 
  markdown: 
    wrap: 72
---

# Slug / Transient-Storage Pipeline — Handoff

Paste as the first message of a new session, with `02_integration.Rmd`
and `03_btc_review.Rmd` attached. Supersedes all earlier handoffs.
Everything is settled unless marked PENDING.

This version's new work is entirely inside Section 9 (co-precipitation)
of `02_integration.Rmd`, mirrored in `R/section09_coprecipitation.Rmd`.
Sections 1-8 and the discharge machinery are unchanged from the last
handoff and are carried forward here without re-verification beyond what
is noted in §4.6/§5.

------------------------------------------------------------------------

## 1. Objective

Kauan Fonseca, doctoral research. NaCl slug additions with paired NH₄-N
and SRP additions in Brazilian karst streams, dry season 2023.

Two deliverables:

1.  **Analysis-ready tables** for inverse transient-storage modelling
    and nutrient uptake. The House (1990) co-precipitation QAQC is now
    closed (§4.6); the next phase is building the storage-zone /
    inverse-TSM uptake calculation itself, in a **new conversation**,
    once Kauan brings the specific method and references he wants to use
    (§6). Do not guess at TASCC/OTIS/inverse-model mechanics before
    that.
2.  **A manuscript** whose novel contribution is combining BTC/TASCC
    with the House (1990) CaCO₃ co-precipitation correction for SRP — no
    precedent in a river breakthrough-curve context. A second,
    independent finding (from a companion project, not this pipeline) is
    that Sw of P correlates with calcium transfer across streams — see
    the circularity caveat in §4.6 and §7 before that finding is written
    up using the *corrected* series.

Each stream is a replicate and none can be discarded, so the approach is
to quantify and flag uncertainty rather than exclude data.

**Pipeline**

| file | role | status |
|----|----|----|
| `01_hobo_qaqc.Rmd` | logger QAQC → `btc_conservative`, `events_hobo` | done |
| `02_integration.Rmd` | nutrients, baseline, pairing, hydraulics, co-precipitation | done, confirm with a clean knit (§5.1) |
| `03_btc_review.Rmd` | hysteresis and visual review | **in progress** — Section 1 built (parity plot + per-event curves with both `conc_corr_ugL` and `conc_corr_coprec_ugL` overlaid); hysteresis loop-geometry analysis not started |
| `R/phreeqc_speciation.R` | phosphate activities, run once | done, fixed this round (§4.6) |
| `R/section09_coprecipitation.Rmd` | standalone Section 9 module, mirrors `02`'s co-precipitation section for the manuscript | kept in sync with `02`'s Section 9 |
| `R/qaqc_discharge_records.R` | campaign discharge records | done |

## 2. Working style — follow these

-   R and tidyverse only. Surgical edits, not rewrites, unless asked.
-   Explicit QAQC flag columns; never silently filter rows.
-   No fabricated citations.
-   Resolve numerical conflicts by physics and by evidence in the raw
    record, never by tuning until a number matches an expectation.
-   Every computed quantity gets a formula and a paragraph explaining
    it, in the style of the discharge sections. Concise, no ornament.
-   Discussion in Portuguese; code, comments and scientific prose in
    English.
-   After editing any chunk: restart R and knit 01 then 02 end to end.
    Partial re-execution silently mixes vintages of objects.
-   **Watch for circularity**, not just unit errors. Any quantity that
    is a deterministic function of another (e.g. `p_coprec_ugL` of
    `delta_ca_mgL`) cannot be used to "confirm" a relationship with that
    other quantity — the relationship is guaranteed by construction, not
    discovered. See §4.6.

## 3. Events

`event_id = stream_YYYYMMDD_{slug_label|station|single}`, 11 in the
pipeline: CB 2023-09-07, CB 2023-10-08, CC 2023-10-06, CD 2023-10-09, RA
2023-09-06 (`_P`, `_N`), RA 2023-10-05 (`_upstream`, `_downstream`), RC
2023-10-07, RS 2023-10-11, SA 2023-10-10 (only non-karst), SR 2023-10-11
(no logger).

Seven received nutrients and carry a calcium pair; those are the
phosphorus comparison set. Nine event_ids actually carry an SRP series
with a defined `p_coprec_ugL` in the current pipeline (§4.6's
mass-balance table lists all nine). The campaign gauged discharge on
other dates too — not in the pipeline, but needed to normalise the
calcium rate.

------------------------------------------------------------------------

## 4. RESOLVED

### 4.1 Bugs in 01 and 02

**`drifting_background()` misaligned end anchor.** In the truncated
branch `y1` was the minimum of the whole window while `t1` was the
window end. CB 2023-10-08 had its minimum at t = +340 s on the rising
limb, so the background fell across the whole tail. Fixed: minimum
searched only after the peak, anchored at its own time, capped so it
never sits above the start anchor. CB 2023-10-08 203 → **317** L/s; CB
2023-09-07 281 → 383.

**Script 01 YAML was `output:a`** — could not knit, so
`btc_conservative.rds` was a stale 9-event vintage.

**`hydraulics` dropped SR** (`left_join` → `full_join`).

**`geometry` had no SR** — built from `events_hobo`, which comes from
the trimmed logger CSVs. Fixed with `GEOMETRY_NO_LOGGER`: 638.9 g, 72 m,
1.5 m.

**`build-events` column collision** produced `.x`/`.y`. Geometry now
enters once via `hydraulics`. Guard:
`names(events) %>% str_subset("\\.(x|y)$")` must be empty.

**`NACL_SLOPE` used before defined.** Own `constants` chunk,
**0.461159**.

**Export wrote before the final build.** `btc_nutrients` is written only
in Section 9.

**`as.Date()` on an integer date.** `as.Date(20231008)` does not error —
it assumes origin 1970 and returns year 57345, so joins returned NA
silently. Use `ymd(as.character(date))`.

### 4.2 Physical findings

**CD 2023-10-09 instrument step.** SpC holds \~432 for 160 s, drops 12
µS/cm in one 10 s interval at constant temperature, holds \~417 until
the tracer arrives at \~800 s. `PRE_STABLE` takes the start anchor from
`c(300, 780)`. 77 → \~46 L/s against 38 from the probe.

**RA 2023-10-05 upstream — incomplete mixing.** Peak ratio 3.61 against
1.03–1.29 elsewhere; the logger sampled the plume core. No valid
discharge.

**RA 2023-09-06 P and N — weak signal.** RA_P peak 3 % above baseline;
noise ±3.5 µS/cm against 24 µS/cm amplitude; 43 % of points zeroed. RA_N
starts at t = +49 s with residual P-slug salt. Both use manual values.

**SR 2023-10-11 discharge = 4.43 L/s — SETTLED, do not reopen.** No
logger; 27 grabs; recession unsampled 1350 s from 1230 s after the peak.
A left-rectangle rule over that gap gives 2.94; the trapezoid gives
4.43. A flat plateau then an abrupt collapse is not how a tracer cloud
passes, so 2.94 is abandoned. Corroborated by the same reach two days
earlier (\~4.3 L/s, different mass, full baseline return). Implied depth
\~33 cm matches field-observed step-pools.

**Nutrient grabs appearing to lead the tracer — investigation PAUSED,
not resolved.** A same-timestamp interpolation check found 11 of 13
event-solute pairs where the nutrient signal rises 1.5–5 minutes before
the conservative tracer. The likely explanation is a field-clock vs.
logger-clock offset — a blind spot the pipeline cannot detect from the
spreadsheet alone, since both series are stamped by the same field
notebook. **Kauan is checking the physical field notebooks
independently; do not touch time-anchoring code on this until he reports
back.**

**SR `addition_time` is not wrong.** t≈0 nearly coinciding with tracer
arrival looked like a possible mis-recorded release time, but is
explained by a sampling-design gap (no nutrient grabs were taken before
first tracer arrival that day), corroborated by two independent field
records. Not the same issue as the paragraph above; no action needed.

**REFUTED — do not re-run.** Sparse-grab interpolation inflating area
(−0.6 %); constant background inflating grabs (made CD worse); clock
offset (translation cannot change ∫C dt); amplitude compression needing
gain (anti-correlates); removing `pmax` (RA_N goes negative);
extrapolating a drift slope from the pre-release (CD predicts 1235 µS/cm
against 400); sensor out of water at the RA tail (0.25 °C, not 25).

**RETRACTED.** Salt-slug-induced CaCO₃ precipitation. The baseline
decline is ordinary diel photosynthesis.

**REJECTED BY KAUAN — do not reintroduce.** An SNR criterion or
`flag_low_snr`. Any composite `autotrophic_signal` from DO: DO confounds
metabolism with reaeration, and in step-pool reaches reaeration alone
can drive it either way.

### 4.3 Discharge source system

`DISCHARGE_SOURCE_OVERRIDE` declares per event which estimate is adopted
and why. Default logger; options `logger`, `probe`, `manual`, `none`.

| event | source | reason |
|----|----|----|
| RA_20231005_upstream | none | incomplete mixing |
| RA_20231005_downstream | probe | start anchor contaminated (logger); flagged, not fallen back to manual |
| RA_20230906_N | probe | no pre-release baseline; flagged (`flag_no_baseline`), not fallen back to manual |
| RA_20230906_P | probe | peak 3 % above baseline; flagged (`flag_weak_signal`), not fallen back to manual |
| SR_20231011_single | probe | no logger; recession unsampled |

Manual values: RA_P 251, RA_N 205, RA_20231005_downstream 201,
RA_20231005_upstream 290 L/s — kept in the data as
`discharge_Ls_manual`, a reference column only; as of Kauan's correction
(§4.8) none of them is the adopted source for any event.

**The probe is not ground truth** — same fixed-background bias, and RA's
between-station spread is 44 %. **The old mass-balance check is
circular**: Q ≡ M/∫C dt, so recomputing recovered mass returns M by
construction. (Contrast with §4.6's new mass-balance check, which is
deliberately built to avoid this: `pct_of_injectate_coprec` uses the
addition record, not the discharge-derived series, as the independent
leg.)

### 4.4 Section 9 — CaCO₃ co-precipitation model (superseded by §4.6)

**Model.** House (1986, 1990),
$n_{P,\mathrm{cop}} = \Delta n_{\mathrm{Ca}}(\sigma N_A d)h$, $h$ from
the two-component Langmuir model on PO₄³⁻ and HPO₄²⁻ activities,
$K_1$/$K_2$ Arrhenius in temperature. d = 2.01e-19 m²/molecule, A₁ =
0.6915, E₁ = 18.2 kJ/mol, A₂ = 4.361e-9, E₂ = 42.6 kJ/mol — all verified
against House (1990) and still correct.

**σ = 5.5e-7 mol/m² is WRONG — superseded, see §4.6.** No literature
support at that order of magnitude; corrected to 5.5e-8 mol/m² (River
Frome, House et al. 1986a; the same value and reasoning Cohen et al.
2013 adopted for the Ichetucknee River).

**PHREEQC speciates the DOWNSTREAM solution only.** Runs once via
`R/phreeqc_speciation.R`. **Two more unit bugs found and fixed this
round — see §4.6**: major ions were passed to PHREEQC 1000× too
concentrated, and PO4/HPO4 activities were used 1e6× too small relative
to House's own µmol/dm³ convention.

**No calcite saturation index.** No carbonate species in the input and
only a campaign-mean alkalinity, so a computed SI would rest on a value
not measured that day. Proposed and removed; do not reintroduce.

**Calcium source is `calcium_reach_change.csv`** — two samples per
event, top and bottom of the reach, both on arrival before release, so
the difference is spatial. Single source for calcium, pH and
conductance. Convention **downstream − upstream, so negative means
calcium lost**; House wants loss positive, so `delta_ca_mgL` is the
negation. Both conventions coexist deliberately, documented at each use.

**Naming.** `conc_corr_coprec_ugL = conc_corr_ugL + p_coprec_ugL`, NA
for NH4-N. **`conc_corr_ugL` is never overwritten** — under inverse
modelling the corrected series *is* the observation vector, so both must
be runnable and the difference attributable.

**Field-sheet error corrected.** `ph` and `temp_c` carried the
*upstream* reading for several events; SR downstream is 7.80 / 26.0, not
7.82 / 27.6. A guard in `R/phreeqc_speciation.R` re-checks every event
against `synoptic_upstream.csv` and warns if any still matches the
upstream value exactly — RS 2023-10-11 currently trips this warning
(ph/temp_c within 0.02/0.2 of the upstream value); **not yet checked
against the raw field sheet, may be a genuine small reach change rather
than a transcription error** (§5.3).

**Log response ratios.**
$\mathrm{LRR}_X = \ln(X_\mathrm{down}/X_\mathrm{up})$ for calcium and
conductance. Log form makes losses and gains symmetric and puts
different scales on one footing. pH is already logarithmic, so the ratio
is on H⁺ activity: $-\ln(10)(\mathrm{pH_{down}} - \mathrm{pH_{up}})$;
negative means alkalinisation. `lrr_ca_excess = lrr_ca − lrr_spc`
isolates what dilution cannot explain. Per-stream aggregation averages
the log ratios (mean of logs = log of geometric mean) and reports the
spread; `same_direction` flags streams whose dates disagree in sign.

**OLD KEY FINDING, needs reconciling — see §5.5.** Across the *full*
14-event calcium survey, most events showed a small calcium *gain*, not
a loss (only SR 2023-10-11 and CD 2023-09-11 unambiguous loss). The
*narrower*, 9-event SRP mass-balance table in §4.6 shows the opposite
pattern — 8 of 9 SRP events show calcium *loss* (`delta_ca_mgL` \> 0 in
House's convention). These are different event sets (14 calcium-survey
dates vs. 9 SRP-paired BTC events) and are not necessarily
contradictory, but nobody has checked whether they actually agree on the
events they share. Do this before writing the manuscript's
co-precipitation framing.

**Analytical precision.** EDTA titration, commercial lab, reported to
0.1 mg/L; several events differ by exactly +0.4. Standard Methods
3500-Ca gives interlaboratory RSD 9.2 % and relative error 1.9 % —
different things; 1.9 % is bias, not a standard deviation. **DECISION:
use 9.2 %**, the only documented value and the conservative one. Now
implemented as Section 9.6 (`ca_delta_sd`, `ca_delta_z`,
`flag_ca_resolved`) — confirmed present in the current file.

### 4.5 Campaign discharge records — `R/qaqc_discharge_records.R`

`discharge_spc_records.csv`, 6590 rows, 20 event-dates.

**READ WITHOUT A LOCALE.** These files use a point decimal. Passing
`locale(decimal_mark = ",")` makes readr treat the point as a thousands
separator and concatenate digits — 1.3833 → 13833, 423.55 → 42355. This
looked like corrupt data for several rounds and was purely a reading
artefact.

Corrections, each declared in a tribble with its reason:

-   **Dates converted automatically when the file was built**, shifting
    months from 08 to 09: CC 02/02 → 02/08, CB 30/09 → **31/08**, RA
    30/09 → **30/08**.
-   **RS 2023-08-29 mass** 1258 → **10291 g**. Gave 153 L/s against 1145
    for the same reach six weeks later from a near-identical curve;
    equal integrals meant the error was entirely in the mass. Now 1253
    L/s, ratio 1.09.
-   **CC 2023-08-02**: two slugs the same day, the first anomalous,
    reach lengthened 126 → 286 m and repeated. Reach length separates
    them in the key; the first is `flag_ancora`, excluded but kept for
    comparison. They give 303 and 264 L/s — agreement that is itself a
    check on the method. The calcium pair corresponds to the 286 m
    reach.
-   **RC 2023-08-10**: logged by both a HOBO and a Hanna. The HOBO
    series is raw conductivity with no temperature compensation, and
    mixing the two gave `retorno = −1.63`. Only the Hanna series is
    used.
-   **`flag_truncado` uses `abs(retorno)`** — a curve ending below
    background is as wrong as one ending above.

**RA 2023-08-04 — adopted 240 L/s.** Lives in its own file,
`discharge_RA_20230804.csv`. First attempt at this site, on a 206 m
reach whose downstream station later became the upstream station of the
sampling reach, containing a pool over 2 m deep. Two problems: the
recession drops straight to zero instead of decaying asymptotically, so
the record was stopped rather than the tracer having finished passing;
and the notebooks disagree on the mass, 6 kg in one and 12 kg in the
other. The integral with 6 kg gives 241 L/s, agreeing with the mean of
the four other RA gaugings (240) — two independent routes to the same
number. The 12 kg reading gives 481, supported by nothing else. **240
adopted, 481 carried in the sensitivity table.** In that file the column
named `Backg_SpC_uScm` is in fact the background-CORRECTED conductance,
not the background.

**Discharge assigned to each calcium sampling** is driven by the calcium
table, same day where it exists, one mapped exception (RS calcium 28/08
uses the 29/08 gauging). Where neither applies the value is NA and
`flag_sem_vazao` is TRUE — never a silent substitution. Output:
`data_derived/discharge_for_calcium.csv`. Consumed by Section 9.7's
areal calcium rate (route 2) — confirmed present and wired up correctly
in the current file.

**Note for the submitted manuscript.** One paper already submitted
reports 481 L/s for RA 2023-08-04 as a site-characterisation value, and
it feeds the calcium precipitation rate used to describe whether the
process is active and its magnitude. It should be corrected to 240 in
revision.

### 4.6 Section 9 continued — three unit bugs, mass-balance verification, and the circularity caveat (this round)

Triggered by Kauan noticing the co-precipitation correction had *no*
visible effect on the SRP parity plot (03 §1.1) and asking whether
PHREEQC's units were right.

**Bug 1 — σ (`SIGMA_P`) was an order of magnitude too high.**
`SIGMA_P <- 5.5e-7` mol/m² had no support in the House/Hartley
literature (published range 5.5e-8 to 2.24e-7 mol/m², clustering into a
field class 0.055–0.095 µmol/m² and a synthetic-solution class 0.13–0.22
µmol/m²). Fixed to **5.5e-8** (River Frome; same value and justification
Cohen et al. 2013 used for the Ichetucknee River, a directly analogous
spring-fed karst system). A new **§9.8 sensitivity table**
(`data_derived/coprec_sigma_sensitivity.csv`) carries the rest of the
field cluster (9.5e-8, Waterston chalk aquifer) and the lab ceiling
(1.5e-7, synthetic Ca(HCO₃)₂ — the value whose \~47× inflation of the
correction and implausible C:P ratios is exactly why Cohen et al. also
rejected it) forward as an explicit envelope.

**Bug 2 — PHREEQC major-ion inputs were 1000× too concentrated.**
`make_solution()` in `R/phreeqc_speciation.R` converted
Ca/Mg/Na/K/SO4/Cl from mg/L to mol/L by dividing only by the molar mass,
skipping the mg→g step (`/1000`). This inflated ionic strength and
ion-pairing enormously, one of two causes of the implausibly small
a_PO4/a_HPO4 PHREEQC was returning. Fixed by adding the missing `/1000`
to each `add()` call. P itself was unaffected (`srp_molL` already did
both steps).

**Bug 3 — a_PO4/a_HPO4 were used 1e6× too small.** House's K1/K2 (Eq.
7–8) and the Langmuir isotherm (Eq. 9) are defined for [PO4³⁻]/[HPO4²⁻]
in **µmol/dm³** — House states this explicitly under Eq. 8. PHREEQC's
`-activities` output is on a mol/L numeric scale, six orders of
magnitude smaller; passed unconverted into `h_langmuir2()`, `h_adsorbed`
collapsed toward \~0 regardless of actual saturation. This was the
larger of the two PHREEQC-side bugs. Fixed:
`h_adsorbed = h_langmuir2(a_PO4 * 1e6, a_HPO4 * 1e6, k1_langmuir, k2_langmuir)`,
in both `02_integration.Rmd`'s `coprec-apply` chunk and its mirror in
`R/section09_coprecipitation.Rmd`.

All three fixes are confirmed present on Kauan's actual working files
(not just a disconnected copy — see the process note below) and
confirmed to knit: Kauan re-ran `R/phreeqc_speciation.R` and reknit `02`
after all three fixes, and the resulting `a_PO4`/`a_HPO4` activities and
`h_adsorbed`/`frac_of_srp` numbers moved by many orders of magnitude
from the pre-fix run, as expected.

**New §9.9 — event-level mass balance.** The per-grab median
`frac_of_srp` (§9.5) is a concentration-domain ratio, not a mass
balance: it weighs every grab equally regardless of discharge or how
much of the curve it represents, and `p_coprec_ugL` tracks SRP
concentration across the curve via `h_adsorbed`. §9.9 instead integrates
by discharge and time, the same construction §6 uses for the NaCl-based
discharge estimate:

$$ m_{P,\mathrm{cop}} = Q \int_0^{T} C_{P,\mathrm{cop}}(t)\,dt
   \qquad
   m_{\mathrm{SRP}} = Q \int_0^{T} C_{\mathrm{corr}}(t)\,dt $$

giving `frac_mass_coprec = m_coprec/m_SRP` (mass-weighted analogue of
`frac_of_srp`) and, independently,
`pct_of_injectate_coprec = 100 * m_coprec / m_P,added`, built only from
the addition record (`added_mass_PO4_g`, `molar_mass_P_nutrient`) and
the coprecipitation model — it never touches `conc_corr_ugL`, so it
cannot inherit an error from the measured SRP series.
`flag_coprec_exceeds_injectate` fires if it ever exceeds 100%
(physically impossible → remaining bug, not a real effect).
`p_coprec_ugL` is integrated signed, not clamped at zero (preserves the
`flag_delta_ca_neg` direction); `conc_corr_ugL` is clamped at zero in
its integral, matching §6's convention for `nacl_mgL_grab`. Output:
`data_derived/coprec_mass_balance.csv`.

Kauan's first run, post-fix:

| event_id | n | discharge_Ls | frac_mass_coprec | pct_of_injectate_coprec |
|----|----|----|----|----|
| SR_20231011_single | 27 | 4.43 | 0.222 | 1.512 |
| CD_20231009_single | 29 | 49.6 | 0.0173 | 0.282 |
| RA_20230906_P | 21 | 251 | 0.00804 | 0.0855 |
| RC_20231007_single | 23 | 501 | 0.00594 | 0.2186 |
| CC_20231006_single | 19 | 194 | 0.00371 | 0.0678 |
| CB_20231008_single | 29 | 317 | 0.00310 | 0.0825 |
| RS_20231011_single | 21 | 1099 | 0.00268 | 0.1154 |
| CB_20230907_single | 13 | 383 | −0.00173 | −0.540 |
| RA_20231005_downstream | 28 | 201 | NA | 0.1409 |

`CB_20230907_single`'s negative values track its negative `delta_ca_mgL`
(reach gained calcium; `flag_delta_ca_neg` — model not parameterised for
that direction, unchanged open question, §7). `RA_20231005_downstream`'s
`NA` `frac_mass_coprec` (but defined `pct_of_injectate_coprec`) is
diagnosed, not a new bug: a `NA` in `conc_corr_ugL` (independent of
`p_coprec_ugL`, which does not depend on measured SRP at all) propagates
through the `mass_srp_exported_ug` integral but not through
`mass_p_coprec_ug`. A new flag, `flag_conc_corr_na_in_window`, names
this so it does not have to be reverse-engineered from a bare `NA`
again.

**Literature cross-check (Cohen et al. 2013; Corman 2024 review, newly
supplied by Kauan).** Cohen et al. (2013), using the *same* σ, found
co-precipitation responsible for \~30% of total P removal in the
Ichetucknee River (assimilation 70%±9%; areal rates 2.7–8.0 mg P m⁻² d⁻¹
co-precipitation vs. 6.4–20.9 assimilation). Our `frac_mass_coprec` is
far below that for 7 of 9 events but lands in the same order of
magnitude for SR (22.2%) — and SR is also the event with by far the
largest `delta_ca_mgL` (10.4 mg/L vs. \<1 mg/L for every other event),
so the model is scaling sensibly with its own driving variable rather
than being uniformly negligible or uniformly dominant. Corman (2024),
citing Bedore et al. (2008), notes that co-precipitation's impact on a
stream's P budget becomes comparatively small under high-P conditions —
directly applicable here, since every one of these events is a P
*injection*, i.e., deliberately elevated far above ambient. Treat this
as a plausibility bound, not a validation of the exact numbers.

**A note on circularity — carries forward to any future work using these
series.** `p_coprec_ugL` (and therefore `conc_corr_coprec_ugL`) is a
deterministic function of `delta_ca_mgL`. Any uptake metric (Sw, Vf,
areal uptake, or otherwise) computed from the *corrected* series
therefore carries a component that is mechanically linear in ΔCa by
construction — independent of any real process. **Comparing a
corrected-series uptake metric back against ΔCa or calcium transfer, and
reporting agreement as confirmation of the co-precipitation mechanism,
is circular.** The evidence that is not circular is (a) the relationship
measured on the *uncorrected* series, and (b) the independent
mass-balance check above. This came up because a companion project (not
this pipeline) found Sw of P correlates with calcium transfer using
*uncorrected* data (real, independent signal — as expected under a
coprecipitation-driven P-limitation hypothesis), but the correction as
originally (buggy) implemented showed no effect, prompting the unit-bug
hunt above; and separately, Kauan reports that in that companion project
Sw's correlation with Ca persists after correction, while Uamb and Vf —
which had *no* relationship uncorrected — gain one only after
correction. That specific pattern (no relationship uncorrected, a
relationship appearing only after a ΔCa-linear correction is added) is
at least as consistent with the correction's own linear structure as
with newly revealed biology, and should be checked against that
possibility (e.g., a synthetic-ΔCa null) before being reported as a
finding. This is a companion-project concern, not something to fix in
this pipeline, but it must travel with the corrected series wherever it
goes next.

**Process note.** Fixes were first made only in a disconnected sandbox
copy of the repo and did not reach Kauan's actual working files; this
was caught when he reran `phreeqc_speciation.R` and got byte-identical
output. All three fixes (plus §9.8, §9.9, and the circularity note) were
then re-applied directly to his real files via targeted find-replace
patches, each verified present via `grep` on the actual device file and
parse-checked as valid R (`parse()` on every extracted chunk, no errors)
— but **never executed**, since PHREEQC and the full data are not
available outside Kauan's machine. The numeric results in the table
above are from Kauan's own console/knit output, not from anything run by
Claude. Keep this distinction in mind: file edits in this pipeline are
sandbox-parse-verified only until Kauan's next real knit confirms them
end to end.

### 4.7 Section 10 rewritten — `master_tsm.csv` replaces `RA_reach_tsm.csv`, no event excluded

Kauan corrected a standing pipeline decision: Section 10 previously
excluded RA_2023-10-05 (both stations) from its one consolidated export,
on the (documented, not silent) grounds that the SRP series that day is
erratic and its NaCl pairing is broken. Kauan's instruction: **keep
every event, exclude none** — consistent with this project's own
standing rule (§2) that was never supposed to allow a silent exclusion
in the first place.

Section 10 is now titled "Consolidated master dataset for inverse TSM"
and builds one table, `master_tsm` / `data_derived/master_tsm.csv` /
`.rds`, across all 11 events — a left join of the full `btc_nutrients`
(every grab, every event, both solutes) against each event's own
`discharge_Ls`/`reach_length_m`/`reach_mean_width_m` **from `events`,
with no override for any event.** An earlier draft of this section
overrode RA_20230906_N/\_P's `discharge_Ls` to a shared reach-level mean
(228 L/s); Kauan corrected that too: **each slug keeps its own
individually-measured discharge** (N = 205 L/s, P = 251 L/s), because
each slug's own velocity/depth (and whatever an inverse model later
derives from its own BTC) were computed from that value in Section 6,
and overwriting it in the master table would silently break that
internal consistency. What the two slugs' reach genuinely shares that
day — one physical river, one survey geometry — is now a separate,
explicitly-labelled reference pair, `discharge_Ls_day_mean` /
`discharge_Ls_day_mean_sd` (228 ± 33 L/s, mean and sd of the two
independent slug-derived estimates), populated only for those two rows
and never substituted into `discharge_Ls` anywhere. Reach geometry (135
m / 9.46 m) needed no override in the first place — `events` already
assigns the same values to both RA dates, confirmed by inspection before
removing the redundant override code.

Every known data-quality problem is now a flag column instead of a
missing row: `flag_discharge_invalid` (RA_20231005_upstream only, no
usable discharge), `flag_mixing_incomplete` (same event, plume-core
sampling), `flag_no_baseline` (RA_20230906_N, RA_20231005_downstream),
`flag_sparse_recession`, and a new `flag_broken_pairing` (both
RA_20231005 stations — the erratic-SRP/broken-NaCl-pairing problem that
used to justify the exclusion). `RA_20231005_upstream`'s `discharge_Ls`
is still genuinely `NA` — that has not changed, and is not fabricated —
but the row itself is present, flagged, for a modelling step to decide
what to do with, rather than never seeing it. Per Kauan,
`RA_20231005_upstream` is expected to sit unused as a reference row for
now — no further action pending on it.

Verified: parse-checked as valid R (161 expressions across 38 chunks, no
errors) and the join logic was independently tested against a small
synthetic stand-in for `events`/`btc_nutrients` in a separate R session,
confirming RA_20230906_N keeps `discharge_Ls == 205`, RA_20230906_P
keeps `discharge_Ls == 251`, both pick up `discharge_Ls_day_mean == 228`
as a reference value only, RA_20231005_upstream keeps `NA` discharge
with its flag and `NA` `discharge_Ls_day_mean`, RA_20231005_downstream
keeps 201 L/s and gains `flag_broken_pairing`, an unrelated event passes
through completely unchanged, and no rows are dropped
(`nrow(master_tsm) == nrow(btc_nutrients)`). **Not yet run against the
real data** — confirm on Kauan's next knit, same caveat as §4.6.

`RA_reach_tsm.csv`/`dict_RA_reach_tsm.csv` no longer exist as pipeline
outputs; superseded by `master_tsm.csv`/`.rds` and
`dict_master_tsm.csv`.

------------------------------------------------------------------------

### 4.8 Discharge source corrected: RA events now adopt probe (not manual); logger/probe both carried in the final dataset; 1:1 comparison plot added

Kauan caught that RA_20230906_N, RA_20230906_P, and
RA_20231005_downstream were adopting the **manual** Excel discharge
instead of the pipeline's own **probe** (handheld Hanna/YSI) trapezoidal
estimate, and was explicit: *"para o RA tem que usar os valores de
discharge devem ser da probe!!! no df final tem que ter os valores do
logger e probe ou uma modo de discriminar qual a fonte do discharge. tb
é bom comparar as estimativas com um plot 1:1."* Asked and confirmed:
switch all three events' adopted source to probe for all three (not a
subset), and keep each event's original physical concern (no baseline /
weak signal / contaminated logger anchor) as an explicit flag rather
than as grounds to fall back to manual.

Changes, all in Section 6 (`hydraulics-combine`) and Section 10
(`master-tsm-export`) of `vignettes/02_integration.Rmd`:

-   `DISCHARGE_SOURCE_OVERRIDE`: the three events are now
    `use_source = "probe"` (was `"manual"`). `discharge_Ls` for each now
    comes from `discharge_Ls_probe`, not `discharge_Ls_manual`.
    `discharge_Ls_manual` stays in the table as a reference-only column
    — see §4.3's updated table.
-   New flag `flag_weak_signal` (RA_20230906_P only): the trapezoidal
    peak is only \~3% above background — previously the stated reason to
    fall back to manual, now an explicit caveat attached to the adopted
    probe value instead of a reason to change source. `flag_no_baseline`
    already existed and already covered RA_20230906_N and
    RA_20231005_downstream's concern — no new column needed for those
    two.
-   New chunk `hydraulics-logger-vs-probe-plot` (Section 6, right after
    the discharge table): a log-log 1:1 scatter of `discharge_Ls_logger`
    vs. `discharge_Ls_probe` for every event with both, coloured by
    adopted `discharge_source`, hover-labelled by `event_id`
    (`ggplotly`). Shows at a glance how far apart the two independent
    trapezoidal estimates sit for every event, not only the three just
    switched.
-   Section 10's `master_hydraulics` now also selects
    `discharge_method`, `reason_discharge`, `discharge_Ls_logger`,
    `discharge_Ls_probe`, `discharge_Ls_manual`, and `flag_weak_signal`
    from `events` — so `master_tsm.csv` alone supports a
    logger-vs-probe-vs-manual comparison, no need to go back to
    `events.csv`.
-   The RA 2023-09-06 day-level reference mean
    (`discharge_Ls_day_mean`/`_sd`, §4.7) is now computed
    **dynamically** from `events$discharge_Ls` for the two slugs
    (`ra_q_lookup` / `RA_Q_VALUES`), instead of the hardcoded
    `c(N = 205, P = 251)` — so it automatically tracks whichever source
    is currently adopted for each slug rather than silently going stale
    the way the hardcoded pair just did. It is very likely no longer
    literally 228 ± 33 L/s now that N/P are probe-based; the actual
    number depends on the real probe trapezoidal integrals, unknown here
    since R/PHREEQC execution never happens in this sandbox — confirm on
    Kauan's knit.

Verified: parse-checked as valid R (39 chunks, 166 expressions, no
errors) and the full source-selection + flag + day-mean logic was
independently tested against a synthetic stand-in for
`hydro_logger`/`hydro_probe`/`events`/`btc_nutrients`, deliberately
using probe values that differ from the old manual numbers (e.g. N:
probe 198 vs. manual 205) specifically to catch any code that might
still be silently keying off the old manual value — confirming
`discharge_source`/`discharge_Ls` correctly switch to probe for the
three events, `flag_weak_signal` is `TRUE` only for RA_20230906_P, the
day-mean reference recomputes from the new probe values (not the old
manual ones), `discharge_Ls_logger`/`discharge_Ls_probe` both survive
into the final joined table even where one is `NA`, and no rows are
dropped. **Not yet run against the real data** — confirm on Kauan's next
knit, same caveat as §4.6/§4.7.

------------------------------------------------------------------------

## 5. PENDING — in order

**5.1 Confirm the clean knit still holds.** The last handoff's 5.1
(`events` was stale; CB 2023-10-08 should read 317 L/s, CD \~46) and
5.3–5.6 (reach width, §9.6 calcium-uncertainty propagation, §9.7 areal
calcium rate, RA_20231005_downstream's stable window) all appear
resolved — the corresponding code exists in the current file and Kauan's
own PHREEQC/knit run this round produced sane discharge values (SR =
4.43 L/s, matching the settled value). Nobody has re-verified the full
checklist end to end since those fixes, though. On the next clean `01` →
`02` knit, re-check: 11 events in `btc_conservative`;
`min(nacl_mgL) == 0`; no `.x`/`.y` in `events`; SR ≈ 4.43 with
`discharge_source == "probe"`; the PHREEQC join guard prints 7 (or the
current expected count), not 0; RC 2023-10-07 (399 vs. 494 L/s, 19%), SA
2023-10-10 (202 vs. 316, 36%), SR 2023-10-10 (5.81, truncated — correct
value \~4.3 from the other sensor), RA 2023-10-05 reach length (137 m
campaign file vs. 135 m in `events`) — these four were flagged
discrepancies in the last handoff and were never confirmed re-checked.

**5.2 RESOLVED — RA 2023-09-06 sample 25 (SRP): the "PHREEQC convergence
failure" was a missing lab value, not a numerical failure.**
`measured_SRP_ug_L` for this one grab (`sample_btc == 25`, t=6780s since
addition, falling limb of RA_20230906_P) is blank in
`nutrient_addition_phosphate_data.csv` -- the SRP lab measurement is
simply missing for this sample (the neighbouring grabs, samples 22-24,
all have valid SRP: 11.86, 11.75, 6.78 ug/L). Because `srp_molL` is `NA`
for this row, `make_solution()`'s `add()` helper
(`if (!is.na(value)) lines <<- c(...)`) never writes a `P` line into
this sample's PHREEQC `SOLUTION` block -- PHREEQC is simply never given
any phosphorus to speciate for this one solution, so the resulting
`a_PO4`/`a_HPO4`/`a_H2PO4` `NA`s are the expected consequence of a
missing input, not an equilibrium/convergence failure. No code change
needed: the existing `p_coprec_ugL`/`conc_corr_coprec_ugL` handling
already produces `NA` correctly for just this one row and excludes it
via the standard `!is.na()` filters -- confirmed directly against
`master_tsm.csv`: all 21 other SRP grabs for RA_20230906_P (samples
3-24, t=720-6300s) have complete `conc_corr_ugL` and
`conc_corr_coprec_ugL` values. RA_20230906_P has a full, usable raw- and
coprec-corrected SRP series; only this one falling-limb point is
missing, same as any other single missing grab.

**5.3 RS 2023-10-11 PHREEQC guard warning.** `ph`/`temp_c` for the
downstream sample are within the guard's 0.02/0.2 thresholds of the
upstream synoptic value (7.57 vs. 7.56; 27 vs. 27 °C). Possibly a
genuine small reach change, possibly a transcription echo of the
upstream value like the ones already caught and fixed for other events
(§4.4). Not checked against the raw field sheet.

**5.4 RESOLVED — Section 10 was rewritten, not just re-worded.** The
stale "negligible (\~1e-8 to 1e-10)" prose no longer exists — the whole
section was replaced (§4.7), not patched. The old worry (RA 2023-10-05's
exclusion resting on the pre-fix bug numbers) is moot: 2023-10-05 is no
longer excluded at all.

**5.5 Reconcile the two calcium-direction findings.** See §4.4's "OLD
KEY FINDING" note: the full 14-event calcium survey found mostly calcium
*gain*; the 9-event SRP mass-balance table (§4.6) shows mostly calcium
*loss*. Check whether the events they share actually agree before
writing this into the manuscript.

**5.6 `flag_conc_corr_na_in_window` fix — pushed, not yet re-run.**
Added to both `02_integration.Rmd` and `R/section09_coprecipitation.Rmd`
this round (§4.6); needs a fresh knit to confirm it actually flags
`RA_20231005_downstream` as expected.

**5.7 Discharge source fix (§4.8) — pushed, not yet re-run.** On the
next knit, check the Section 6 printed table and the new
`hydraulics-logger-vs-probe-plot`: `discharge_source` should read
`"probe"` for RA_20230906_N, RA_20230906_P, and RA_20231005_downstream;
`flag_weak_signal` should be `TRUE` only for RA_20230906_P;
`discharge_Ls_day_mean` in `master_tsm.csv` (Section 10) should no
longer be exactly 228 (it now derives from the two slugs' probe values,
not the old manual 205/251) — note whatever the new number actually is,
since §4.8's text above only states the mechanism, not the resulting
value.

**5.8 RESOLVED (2026-09-23) — the TSM/uptake batch pipeline was not
reproducible across machines/R installations; fixed via per-job
seeding.** Comparing `data_derived/tsm_uptake_summary.csv` before and
after the `ca2d3a2` "handoff updates" commit (its parent `d286ddb`
changed nothing in `R/` or `master_tsm.csv`) showed that re-running the
identical `run_tsm_uptake_all(seed = 1)` on a different machine changed
`uptake_significant`/`lambda_at_bound` for several event x solute rows,
including `RS_20231011_single` NH4-N going from `pct_total_uptake`
0.01% to 26.57%. Root cause: `set.seed(seed)` was called once, before
looping over all events (Stage 1) or all event x solute x correction
jobs (Stage 2), so every job's Latin Hypercube draw depended on
however many random draws every job before it had consumed -- and
`sample()` (used inside `lhs_sample()`, `tsm_calibrate.R`, to break
correlation across LHS parameters) is only guaranteed to reproduce
identically for a given seed on the same R version/platform (R's
default `sample.kind` changed in R 3.6.0). Fixed by giving every job
its own seed, derived deterministically from a stable string key
(`job_seed(seed, key)`, new in `tsm_calibrate.R`) instead of from call
order -- `fit_all_hydraulics()` now seeds each event off its
`event_id`, `run_one_uptake()`/`run_tsm_uptake_all()` off
`"event_id|solute|conc_col"`, and `app_tsm_review.R`'s Refit buttons do
the same (base seed 1), so repeated Refit clicks with identical inputs
now return bit-identical fits instead of silently re-sampling. This
does not fix non-identifiability itself -- it was, unsurprisingly,
exactly the events already flagged `lambda_at_bound`/`lambda_s_at_bound`
that moved the most, since a flat/near-boundary objective function is
what makes an optimiser's outcome sensitive to tiny numerical
differences in the first place. **Any batch number for a
borderline/non-identifiable event from before this fix (i.e. anything
in the current `tsm_uptake_summary.csv` as of 2026-09-23 or earlier)
should be treated as unconfirmed until re-run with the current code**;
see the equivalent note added to vignette 04's "Notes and known
limitations". No `renv.lock`/`sessionInfo()` capture exists yet to
pin down the exact environment difference between the two runs (still
recommended as a follow-up, not yet done).

**5.9 Rafa's exact-solver Stage-1 rewrite (`ded1543`, `0ac9cfb`) --
spot-verified, one open reproducibility gap.** Rafa replaced Stage-1's
finite-volume/ODE solver with an exact analytical one (`tsm_model.R`:
Laplace-domain convolution via `besselI`, adaptive node spacing sized to
both the ADE pulse and the storage kernel -- the earlier version
under-resolved fast-exchange storage zones, e.g.
`RA_20231005_downstream`'s As/A = 0.0017), added pre-arrival baseline
zeroing (`zero_pre_arrival()`, a stopgap for logger drift like
`CB_20230907_single`'s 0 -> 0.7 mg/L creep before arrival), and made
velocity a second Stage-1 candidate: every event is now fit twice
(`v_fixed` = `events$water_velocity_ms`, and `v_fitted` = v as a 4th
free parameter, Q always fixed from dilution gauging), with the model
chosen automatically by AICc (`aicc_min_gain = 2`,
`fit_all_hydraulics()`/`TSM_METHOD` in `run_tsm_uptake.R`). Rationale:
peak-time velocity systematically underestimates true channel velocity
once storage delays the peak. Also fixed a real bug (bound-clamped
internally but the unclamped optimiser vector was what got reported)
and added a hard-coded `logger_rescaled` treatment for `RA_20230906_N`/
`_P` specifically (logger SHAPE + probe-consistent mass, since the
logger recovered only 35%/61% of the injected salt --
`diag_ra_logger.R`).

Verified 2026-09-24 (fresh machine, R 4.3.3, same input CSVs):
`tsm_selftest()` passes cleanly, including the new fast-exchange
convergence check Rafa added (As/A = 0.0017/0.02/0.5, all < 5e-4
relative). Re-running `fit_all_hydraulics()` at production settings
(seed 1, n_lhs 250) for `RC_20231007_single` reproduced the committed
`v_ms`/`D_m2s`/`alpha_1s`/`As_over_A`/`hydraulic_rmse` to ~1e-5 relative
difference -- the `job_seed()` reproducibility fix (§5.8) holds up
under the new solver for a well-identified event.

**Not holding up: `RA_20230906_N`/`_P`.** The same rerun gave
`v`/`D`/`alpha`/`As_over_A` 20-70% off from the committed
`tsm_hydraulics_summary.csv`, and for `RA_20230906_P` the AUTOMATIC
MODEL SELECTION ITSELF FLIPPED -- committed run: `v_fixed` (dAICc =
-1.3, barely below the keep-v_fixed threshold); this rerun: `v_fitted`,
decisively (dAICc = +43.2). These are exactly the two events Rafa's own
`diag_ra_logger.R` already flagged as poorly constrained (probe-grab
fits disagree with each other: v 0.051 vs 0.038, As/A 0.67 vs 0.37).
Read together: for a borderline-identifiable event, which machine/
R-version runs the batch can change which hydraulic MODEL gets used for
the final result, not just its parameter values -- a step beyond the
§5.8 caveat (which was about parameter drift within the same model).
Not re-run for the other 11 events (a time-boxed spot check, not a full
batch rerun). Recommend before treating `RA_20230906_N`/`_P`'s current
numbers as final: either pin their hydraulic model/parameters manually
(the Shiny-review flow already has one open case needing exactly this --
`CB_20230907_single`'s accepted hydraulics used `conservative_source =
logger`, but Kauan had described accepting `probe_grab` mid-session;
not yet confirmed which was intended) or re-run with a wider/more
robust search (more LHS draws or starts) specifically for these two and
confirm the AICc gap is no longer near the threshold.

**Also stale:** the committed `data_derived/tsm_hydraulics_summary.csv`
is missing two columns (`logger_recovery`, `nacl_mass_fit_g`) that the
current `tsm_hydraulics_table()` (as committed in the same `ded1543`)
does produce -- the file was not regenerated after the code's last edit
before that commit. Not touched here (see the reproducibility gap
above -- regenerating it now would just add a third machine's numbers
to the mix); worth a clean re-run once the RA question above is
settled.

------------------------------------------------------------------------


## 6. NEXT — storage-zone / inverse transient-storage modelling

**Start a new conversation for this.** Bring `02_integration.Rmd`,
`03_btc_review.Rmd`, and this handoff, plus **the specific
uptake-calculation method and references Kauan wants to use** — do not
let a new session guess at TASCC/OTIS/inverse-model mechanics before
that; his TASCC/Sw work for the companion project already exists
elsewhere and should not be rebuilt here from assumptions.

**What is ready:**

-   `data_derived/btc_nutrients.csv`/`.rds` — one row per grab, both
    solutes, all events, with `conc_corr_ugL` (background-corrected)
    and, for SRP, `conc_corr_coprec_ugL` (+ co-precipitation correction)
    side by side, paired `nacl_mgL` (conservative tracer),
    `time_since_release_s`, `curve_limb`.
-   `data_derived/btc_conservative.csv`/`.rds` — full-resolution logger
    series per event.
-   `data_derived/events.csv` — event-level discharge (`discharge_Ls`),
    reach geometry, and QAQC flags.
-   `data_derived/coprec_mass_balance.csv` — event-level
    co-precipitation mass balance (§4.6), for deciding how much the
    correction should matter to any given event's uptake estimate.
-   `data_derived/master_tsm.csv`/`.rds` — the consolidated master
    export (§4.7): every grab, every event (all 11, none excluded), both
    solutes, with discharge/geometry and every QAQC flag joined in. This
    is very likely the single file an inverse-TSM tool should read from;
    build a narrower per-event view from it if the method needs one,
    rather than re-deriving discharge/geometry again.

**Carry forward into that conversation:** the circularity caveat
(§4.6/§7) — whatever uptake metric gets computed from
`conc_corr_coprec_ugL`, do not treat its relationship with ΔCa or
calcium transfer as new confirming evidence of the co-precipitation
mechanism.

## 7. Open scientific questions

**ΔCa \< 0 treatment.** Where the reach gains calcium the model returns
a negative $n_{P,\mathrm{cop}}$ and the correction subtracts phosphorus,
outside the direction House parameterised. No precedent in the BTC
literature. Agreed approach: sensitivity analysis across treatments
(exclusion, no correction, zeroing via `pmax`) rather than asserting one
as correct. Still open; §5.5 may change how central this is once the
calcium-direction reconciliation is done.

**RA 2023-09-06 instrument disagreement.** At background before salt
arrival the probe reads 466 while the logger reads 502 (P) and 483 (N).
Every refuted mechanism was tested against this. Unresolved.

**Transient-storage parameters for SR.** Kauan's preference: use the
2023-10-09 high-frequency curve directly for transient storage and
reserve the 2023-10-11 grabs for uptake. Either way **the discharge used
for uptake must be the 2023-10-11 value**, since that is the day of the
nutrient slug.

**Sw_P vs. calcium transfer, and the metric-dependent pattern after
correction.** Companion-project finding (§4.6): Sw's correlation with Ca
persists after the House correction; Uamb and Vf, which show no
relationship uncorrected, gain one only after correction. Flagged as
plausibly explained by the correction's own ΔCa-linear structure rather
than newly revealed biology (§4.6's circularity note) — worth a
synthetic-ΔCa null check before it goes into any write-up. Not something
to resolve in this pipeline; noted here so it is not lost between
conversations.
