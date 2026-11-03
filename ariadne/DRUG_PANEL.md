# Expanding the drug panel: what the model needs, drug by drug

Status 31 Oct 2026. Modelled: **8 drugs**. Target asked for: 20.

This file is the shopping list. It exists because "20 drugs" is a UI number, not a science number —
the honest ceiling with today's model is about 12–13, and the rest need either a new model term or an
assumed value. Both options are written out below so the choice is explicit.

---

## 1. What ARIADNE needs from a drug

The model has exactly **two** input channels. Everything a drug does must land in one of them:

| Channel | What feeds it | Formula |
|---|---|---|
| `x0_effect` (excitability) | sodium-channel **inactivated-state block fraction** `f_Na`, **GABA-A potentiation fraction** `g_GABA`, plus any residual `other_x0` | `(f_Na + ρ·g_GABA) / NAV_REF + other_x0`, ρ = g_I/g_E = 1.0 (derived) |
| `coupling_effect` (spread) | **transmitter-release reduction**, either measured directly or derived from Ca²⁺ current block via `1 − (1 − block)^2.6` | used as a multiplier on connectome weights |

So a complete drug entry is **five numbers**:

| # | Quantity | Usual source |
|---|---|---|
| A | molecular weight (g/mol) | any label |
| B | therapeutic plasma reference range (mg/L) | Patsalos et al. ILAE TDM guidelines; FDA label |
| C | plasma protein binding % → free fraction `fu` | same |
| D | **at least one mechanism number**: an IC50 in µM from a named assay, or a % effect at a stated concentration | primary patch-clamp literature |
| E | therapeutic ceiling (top of range, or max licensed dose) | FDA label |

A, B, C and E are easy and come in bulk. **D is the whole problem.**

---

## 2. Modelled now (8)

| Drug | Class | Mechanism source | Status |
|---|---|---|---|
| carbamazepine | Na⁺ blocker | Kahlig 2022 Table 1, VDB IC50 44.37 µM | complete |
| lamotrigine | Na⁺ blocker | Kahlig 2022, 39.09 µM | release reduction 0.3 still **assumed** |
| phenytoin | Na⁺ blocker | Kahlig 2022, 47.78 µM | complete; non-linear PK caveat recorded |
| lacosamide | Na⁺ blocker (slow inactivation) | Kahlig 2022, 269.3 µM | complete; slow-inactivation caveat recorded |
| cenobamate | Na⁺ + GABA-A | Kahlig 2022, 66.71 µM + GABA potentiation 45% | complete |
| valproate | broad | Kahlig 2022 (18% inhibition @ 1 mM) | `other_x0` 0.4 and release 0.2 still **assumed** |
| pregabalin | α2δ / release | — | release reduction 0.6 **assumed**, no measurement |
| levetiracetam | SV2A + N-type Ca²⁺ | Lukyanetz 2002, 18% / 37% Ca²⁺ current block | **lower bound**: SV2A has no published effect size |

---

## 3. Reachable with ordinary PK data only (D already in hand)

| Drug | Have | Need | Blocker |
|---|---|---|---|
| **oxcarbazepine** | VDB IC50 **42.0 µM** (Kahlig 2022) | A, B, C, E | **Real problem.** Kahlig measured the *parent*. In patients the active circulating species is licarbazepine (MHD, the monohydroxy derivative), which is what the reference range refers to. Pairing the parent's IC50 with MHD's concentration would be wrong. Either find an **MHD inactivated-state IC50**, or model parent-only and publish it as a lower bound. Not added yet, deliberately. |

---

## 4. Reachable with one mechanism number plus PK

| Drug | Class | Mechanism number needed | Notes |
|---|---|---|---|
| **phenobarbital** | GABA-A | GABA-A current **potentiation fraction** at a stated concentration | NG217 5.2.6 (focal third-line add-on) and 5.1.3 (GTC). The cleanest GABA-only drug to add — it would test the derived ρ independently of cenobamate. |
| **eslicarbazepine acetate** | Na⁺ | inactivated-state IC50 (µM) for **eslicarbazepine** on a comparable preparation | NG217 5.2.5, same tier as cenobamate — a direct head-to-head. |
| **zonisamide** | Na⁺ + T-type Ca²⁺ | Na⁺ inactivated-state IC50 **and** T-type Ca²⁺ block % at a stated concentration | NG217 5.2.2/5.2.4. The T-type part needs a mapping decision (see §6). |
| **topiramate** | Na⁺ + GABA-A + AMPA | Na⁺ IC50 **and** GABA-A potentiation fraction | AMPA antagonism has no model term; would be an explicit omission, like levetiracetam's SV2A. |
| **clobazam** | GABA-A | potentiation fraction for clobazam **and** for N-desmethylclobazam | NDMC circulates at roughly 8× the parent, so both are needed or the drug is mis-modelled. |
| **clonazepam** | GABA-A | potentiation fraction | Simplest benzodiazepine case. |
| **ethosuximide** | T-type Ca²⁺ | T-type current block % at therapeutic free concentration | Mapping decision needed (§6). Absence-only clinically, so it never competes in the focal ranking. |
| **gabapentin** | α2δ / release | a **measured** Ca²⁺ current block or release reduction | Adding it on pregabalin's assumed 0.6 would add a drug and no evidence. Worth doing only with a real number — which would also replace pregabalin's assumption. |
| **brivaracetam** | SV2A (~20× levetiracetam) + Na⁺ | Na⁺ IC50; SV2A effect size does not exist | **Highest scientific value of anything on this list.** Brivaracetam vs levetiracetam is a direct test of whether the SV2A gap explains ARIADNE ranking levetiracetam last. |

---

## 5. Not modellable without a new model term — list them, do not fake them

| Drug | Mechanism | Term the model would have to grow |
|---|---|---|
| perampanel | AMPA antagonist | AMPA-mediated excitatory coupling (arguably the most defensible addition: it is a coupling reduction) |
| vigabatrin | irreversible GABA-transaminase inhibition | ambient/tonic GABA elevation; acts over days, so no instantaneous concentration→effect map |
| tiagabine | GAT-1 reuptake inhibition | IPSC decay prolongation → effective inhibitory charge |
| felbamate | NMDA antagonist + GABA-A | NMDA term |
| stiripentol | GABA-A + CYP inhibition | its clinical role is a PK interaction with clobazam; outside a monotherapy model |
| rufinamide | prolongs Na⁺ inactivated state | would fit the Na⁺ channel if an IC50 existed; reportedly weak |
| cannabidiol, fenfluramine, piracetam | GPR55/TRPV, 5-HT/σ1, unclear | no quantified mechanism to map |

---

## 6. Two mapping decisions that must be made before adding drugs

1. **T-type calcium (zonisamide, ethosuximide).** The release-reduction path assumes *presynaptic* Ca²⁺ controlling
   transmitter release, with cooperativity 2.6. T-type channels are largely **somatodendritic** and control burst
   firing, not release. Routing T-type block through the coupling term would be mechanistically wrong. It belongs in
   `other_x0` (excitability), and the conversion from "% T-type current blocked" to "x0 shift" needs its own
   justification — otherwise it is just another assumed number wearing a citation.

2. **Prodrugs and active metabolites** (oxcarbazepine→MHD, primidone→phenobarbital, clobazam→NDMC,
   eslicarbazepine acetate→eslicarbazepine). The model takes one free concentration per drug. Each of these needs
   either the metabolite's own IC50 at the metabolite's own concentration, or an explicit statement that the parent
   is being modelled alone.

---

## 7. The efficient ask: four papers, not twenty

Nearly all of columns A, B, C and E for every drug above sit in **one** source, and column D sits in about three:

1. **Patsalos PN et al., "Antiseizure medications: best practice guidelines for therapeutic drug monitoring",
   Epilepsia / Ther Drug Monit (2018 update).** One table: reference range and protein binding for essentially
   every drug listed here. This single paper unlocks B, C and E across the board.
2. **A comparative sodium-channel patch-clamp study covering eslicarbazepine, zonisamide, rufinamide and
   topiramate** on an inactivated-state protocol — ideally the same hNaV1.6 / PatchXpress setup as Kahlig 2022,
   so the numbers are commensurable with the four already in use.
3. **A comparative GABA-A potentiation study** reporting % current potentiation at stated concentrations for
   phenobarbital, clobazam, clonazepam and topiramate. One study covering several of them is worth far more than
   four separate ones, because the assay conditions then match.
4. **A T-type Ca²⁺ IC50 paper** for ethosuximide and zonisamide — only needed if decision 6.1 is resolved.

Plus, if available: **brivaracetam sodium-channel IC50**, and an **MHD / licarbazepine inactivated-state IC50**.

---

## 8. Recommendation

Do not pad to 20. A panel of 12–13 drugs where **every number has a citation** is a stronger exhibit than 20 where
eight are guessed, and the guessed ones would be the first thing a pharmacology judge asks about.

The strongest single addition is **brivaracetam**, because it turns the project's biggest embarrassment
(levetiracetam ranked last) into a designed experiment. The second strongest is **phenobarbital**, because it tests
the derived GABA weight ρ on a drug that is nothing but GABA-A.
