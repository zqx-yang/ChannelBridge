# Verification against iMARS

Wang J, Ouyang X, Meng S, ..., Shi T, Zhao Y-L, Ni J (2025). **Rational
multienzyme architecture design with iMARS.** *Cell* 188:1349–1362.
Audit performed 2026-10-23 against the paper and its Tables S1–S7.

ChannelBridge uses the iMARS DO score. This document records how much of that
has been checked against the paper's own published data, what has not, and
what it would take.

## 1. Are we doing the same thing?

Yes, to a closer degree than expected: iMARS uses the same two tools.

| | iMARS | ChannelBridge |
|---|---|---|
| structure | AlphaFold2 via ParaFold, top of 3 models, pLDDT ~90 | AlphaFold2/3, Boltz (prep stage) |
| active site | structure alignment to a homologous complex, else AutoDock Vina | explicit `site` selection in the pair file |
| channels | **CAVER 3.0**, run from PyMOL on a local server | **CAVER 3.0.3**, run from a config |
| channel chosen | "the top-ranked channel in each enzyme" | pair file `caver.cluster`, chosen on convergence + cofactor proximity + burial |
| D_IC | distance between the outermost sphere centres of the two tunnels | same |
| O_IC | angle at the upstream outermost sphere, arms to the upstream second-to-last sphere and to the downstream outermost sphere | same, but the arm is taken over a fixed arc length (see §5) |
| DO | `-0.0079*D_IC + 0.0012*O_IC + 0.7865`, confined to [0,1] | same, coefficient set recorded per run |
| what is varied | the linker of a fusion protein | the pose of a binder-held pair |

The geometry definitions are quoted verbatim from the paper's Methods
("Molecular docking and channel calculation") and Figure 1D legend:

> Segment X indicates the direction of intermediate molecule movement at the
> end of tunnel A. Segment Y indicates the connection between exits of tunnel
> A and tunnel B. D_IC is the length of Segment Y. O_IC is the angle between
> line segments X and Y.

> The top-ranked channel in each enzyme is selected for D_IC and O_IC
> calculations. [...] The line connecting the centers of the outermost sphere
> of upstream and downstream enzymes is used to calculate the D_IC value. For
> the O_IC value, starting from the second-to-last sphere of the upstream
> enzyme, the centers of the outermost spheres of upstream and downstream
> enzymes are connected sequentially.

So our construction is theirs. The one substantive difference is §5.

## 2. What is verified

The paper publishes D_IC, O_IC, DO and the measured product yield for 130
fusion enzymes, with no coordinates. That is enough to verify the whole
**scoring layer** end to end, and `tests/test_imars_verification.py` does so.

**Their yield regression reproduces exactly.** Refitting their Table S2A
(100 4CL-STS fusions) gives

```
Pr = -0.7101*D_IC + 0.1070*O_IC + 55.105     R^2 = 0.7350
paper:  -0.71          0.107        55.10    R^2 = 0.73
```

**The DO derivation reproduces.** The paper says DO is that regression
"normalized to confine the range between 0 and 1", with 0 at D_IC = 100 Å,
O_IC = 0°. Those corners give Pr = −15.90 and +74.36, a range of 90.26, so

```
a = -0.71/90.26  = -0.0078662       paper states -0.0079
b =  0.107/90.26 = +0.0011855       paper states +0.0012
c =  71.0/90.26  = +0.7866164       paper states +0.7865
```

The slopes round to the stated values exactly. The intercept rounds to 0.7866,
not 0.7865 — a 1 × 10⁻⁴ discrepancy in the paper, which shifts every DO
equally and changes no ranking.

**Our implementation matches their numbers.** Scoring their Table S2F (the
independent 30-linker validation set) with the unrounded constants reproduces
their DO column to **5 × 10⁻⁵**.

**DO does predict yield, including at our distances.** On their fit set
r² = 0.735; on the independent validation set r² = 0.738, and that set reaches
down to D_IC 17.4 Å — the fusion with the highest DO (0.827) is also the best
measured (61.3 mg/L).

## 3. Three inconsistent DO columns in the supplement

Worth knowing before trusting any single iMARS table.

| table | what its DO column was computed with | agreement |
|---|---|---|
| **S2F** (30-linker validation) | the unrounded derivation above | 5 × 10⁻⁵ |
| **S2A** (100-fusion fit set) | `-0.00801133 / +0.00142190 / +0.77969380` | 4 × 10⁻⁷ |
| **S2G** (32 predicted pairs) | a mixture, plus a clamp at 0 | no single set fits |

The S2A set does **not** derive from any regression on its own table. It
implies 5.63° of O_IC per Å of D_IC, where the table's own yield regression
gives 6.63; no subset, leave-one-out or log transform reproduces 5.63. It also
predicts the measured yields very slightly worse (r² 0.7323 vs 0.7350). We
treat it as an error in that table, available as `imars-st2a` for reproducing
it but not recommended for scoring.

S2G reuses S2A's values for the 4CL-STS pair and clamps negatives: its row at
D_IC 116.5 Å / O_IC 68.5° reads exactly `0.000000` where the unclamped value
is −0.0488. Hence `do_clamped()`.

**Conclusion: the constants ChannelBridge defaults to are the right ones.**
They are what the paper states in three places (Figure 1H, Methods,
Quantification) and what its own regression derives. Secondary sources
summarising the paper quote the same three constants and the same Pearson
r = 0.86 against yield, which we reproduce as 0.857.

## 4. What is **not** verified

That our **geometry extraction** reproduces theirs. Everything above tests
`DO(D_IC, O_IC)`. Nothing tests `(D_IC, O_IC) ← structure`, because the paper
publishes no coordinates — the Methods note there are "no available
co-crystal structures in the PDB database" for these fusions, and the
predicted models are not deposited.

This is the gap that matters, because it is where the CAVER-build sensitivity
documented in [caver-builds.md](caver-builds.md) lives.

### Route A — ask their server: **ruled out**

iMARS is nominally live at `http://www.imars-lumybio.com/` (HTTP only; HTTPS
refuses the connection), behind a sign-up, advertising "DO score evaluation".

**It does not return results.** A GOx–HRP prediction task was submitted on
2025-08-20. The acknowledgement said:

> Due to current server overload, the estimated queue time is over one month.
> The prediction results will be returned to you via email.

Nothing arrived in the following fourteen months. The contact address has also
changed from the paper's `imars_service@lumybio.com` to a lab Gmail address on
the current site. The web front end is up; the prediction pipeline behind it
is not delivering, so this route cannot be relied on for validation.

### Route B — rebuild their structures (the only viable route)

Everything needed is published:

* **sequences** — Table S2G gives accession numbers for 32 pairs; the
  4CL–STS system used throughout is `AY376729` (At4CL) and `ACC76753`
  (AhSTS).
* **linkers** — Table S3 gives the 50-, 200- and 1,025-linker libraries, and
  Tables S2A/S2F name the linker for every measured fusion.
* **construction** — Table S7B is their script: `enzyme1 + linker + enzyme2`,
  plain concatenation.
* **ground truth** — Tables S2A and S2F give D_IC and O_IC for 130 fusions.

So: build the 100 fusions, predict each with AlphaFold2, run CAVER, extract
D_IC/O_IC with ChannelBridge, and regress our values against their 130. That
is a real test with 130 points, on the same enzyme pair, against numbers
produced by the same two tools.

What it would cost and where it could go wrong:

* 100 fusions of roughly 800 residues: a day or two of GPU time with ColabFold or AlphaFold3.
* **Their CAVER settings are not published.** Probe radius, shell radius,
  shell depth and clustering threshold are all absent. A mismatch there
  changes the sphere chain, which is exactly the sensitivity we measured.
  Mitigation: sweep the settings and report which reproduce their D_IC — that
  turns the unknown into a measurement, and would incidentally recover the
  settings iMARS used.
* **Their starting point** comes from active-site localisation by homology
  alignment or Vina docking. We would use the known 4CL and STS catalytic
  residues instead, which should agree but is a second degree of freedom.
* AlphaFold version drift: they used AlphaFold2 via ParaFold in ~2024. Our
  models will differ somewhat, which puts a floor on the achievable
  agreement. Their own reported pLDDT (~90) gives a sanity target.

A sensible reduced version: do the **30** fusions of Table S2F rather than
all 100. They span D_IC 17–90 Å, they are the set with independent measured
yields, and 30 structures is a few hours rather than days.

#### Start with n = 1

Before committing GPU time, do the single best-characterised fusion:

```
STS-(GGGGN)4-4CL      linker GGGGNGGGGNGGGGNGGGGN   (Table S2F row 25)
iMARS reports:        D_IC 17.386 A   O_IC 148.008 deg   DO 0.825
measured:             61.3 mg/L resveratrol, their best fusion
```

Build `AhSTS + GGGGNGGGGNGGGGNGGGGN + At4CL` from `ACC76753` and `AY376729`,
predict it, run CAVER, extract D_IC/O_IC. One structure, about an hour.

The outcome is informative either way:

* **within a few Å / tens of degrees** — our extraction agrees with theirs on
  a real case, and scaling to 30 points is worth it.
* **far off** — compare before blaming the extraction: check the fold is right
  (their pLDDT was ~90), then sweep the CAVER settings, since those are
  unpublished and are the known sensitivity. A settings sweep that lands on
  their number also recovers the settings iMARS used, which is worth having
  on its own.

Note this fusion sits at D_IC 17.4 Å, which is both their validated floor and
almost exactly where this project's best pose sits (17.0 Å) — so it is the single most
relevant point in their entire dataset for our purposes.

### Route C — what cannot be done

Verifying against the individual enzymes' crystal structures (4CL and STS both
have PDB entries) does not help: the paper reports no tunnel geometry for the
isolated enzymes, only for fusions. There is no published per-enzyme number to
check against.

## 5. Where we knowingly differ

Three deliberate deviations, each recorded so a reviewer can see them.

**The O_IC arm.** iMARS takes P1 as literally "the second-to-last sphere".
We take it a fixed 1.0 Å of arc length inward. On a CAVER build with ~1 Å
sphere spacing these are the same thing, and the reference build had 0.98 Å
spacing —
so this reproduces their intent. On the 0.17 Å-spacing build available here,
taking the literal second-to-last sphere would measure the angle over a
0.26 Å baseline. See [caver-builds.md](caver-builds.md).

**Channel selection.** iMARS always takes CAVER's top-ranked cluster. The
reference run of the worked example chose the FNR channel empirically (cluster 2, burial 18, over cluster 1,
burial 37) because it scored better in the pose search. That is a departure
from the published method and should be stated in any methods text. The pair
file records the choice explicitly.

**Extrapolation range.** The DO regression was fitted on fusion proteins whose
channel exits were **45–91 Å** apart and validated against yields down to
**17.4 Å**. A binder-held pair can be much closer: our pose grid runs 6–36 Å,
and the worked example's best pose sits at 17.0 Å — just below the validated floor. DO is a
linear fit with no support below that; it keeps rising as D_IC falls, which is
an assumption rather than a measurement. `extrapolation_note()` flags it, and
a methods section should quote the fitted range alongside any DO.

## 6. Reproducing this audit

```bash
export CHANNELBRIDGE_IMARS_ST2=/path/to/1-s2.0-S0092867424014715-mmc3.xlsx
pip install openpyxl
pytest tests/test_imars_verification.py -v
```

The workbook is copyrighted and is not redistributed with this repository;
the tests skip without it.
