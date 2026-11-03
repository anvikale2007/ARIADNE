# ARIADNE - abstract (reframed, 31 Oct 2026)

The old abstract claimed ARIADNE "ranks antiseizure drugs for a patient". That overclaims. This version
claims exactly what the evidence supports, which is both more defensible and more interesting.

---

## Title

**ARIADNE: translating single-channel drug potencies into whole-brain seizure spread, with the
uncertainty made explicit**

## Background

About a third of people with epilepsy continue to have seizures despite treatment, and antiseizure
medication is still chosen largely by sequential trial and error. The evidence that guides that choice
sits at two disconnected scales. Drug potencies are measured on isolated ion channels and receptors.
Focal seizures, however, are a network phenomenon: they begin in one region and propagate along
white-matter pathways. No routine method connects a drug's measured molecular action to the spread of a
seizure through an individual brain.

## Aim

To build and critically evaluate a framework that translates published drug potencies, at real
unbound plasma concentrations, into their predicted effect on seizure propagation through a patient's
own structural connectome - and to quantify how much of the resulting answer rests on measurement
rather than on assumption.

## Methods

Diffusion MRI was converted to a structural connectome by deterministic tractography (DIPY). A reduced
Epileptor neural-mass model, verified against The Virtual Brain implementation (agreement to ~1e-19),
was run in every region, coupled through the connectome, and calibrated to each brain's own
seizure-spread threshold. Six antiseizure drugs (cenobamate, lamotrigine, valproate, carbamazepine,
pregabalin, levetiracetam) were represented through two mechanisms: reduced excitability, sized from inactivated-state
sodium-channel block and GABA-A potentiation by Hill occupancy at published free concentrations, and
reduced coupling, representing decreased transmitter release. A seizure-type gate implementing the ILAE
2025 classification and NICE NG217 restricts the analysis to focal seizures and flags drugs that may
worsen other seizure types; pharmacogenomic rules exclude unsafe drugs.

Because every drug enters the model through only two quantities, the boundary between propagation and
containment was mapped once over that two-dimensional plane (1701 simulations), making each drug's
threshold a geometric lookup. This enabled Monte Carlo uncertainty propagation over 2000 parameter
draws in two tiers: sampling only ranges with published sources, then additionally sampling parameters
with no source. Analyses were run on two healthy-volunteer scans and two epilepsy patients (IDEAS II,
OpenNeuro ds007401; Penn 3T-7T). A virtual-surgery test against real surgical outcomes was
pre-registered before being run.

## Results

A drug acting purely by reducing transmitter release would need a 175% reduction to match lamotrigine and 281% to match cenobamate - both impossible. Levetiracetam nonetheless works clinically, which identifies the coupling formulation as the model's clearest point of failure.

In the headline brain (Stanford HARDI, 70 regions), containment thresholds were cenobamate 0.48, lamotrigine 0.76, valproate 1.50, carbamazepine 1.59 and pregabalin 2.23 multiples of standard exposure. Thresholds read from the response surface agreed with independent full simulations to within 0.07, the bisection resolution.

Containment thresholds were expressed against each drug's therapeutic ceiling rather than as abstract
dose multiples. Because receptor occupancy saturates, the maximum licensed dose delivers substantially
less than a proportional increase in effect: cenobamate and lamotrigine reach containment within their therapeutic ceiling, while valproate, carbamazepine and pregabalin do not, because a doubling of dose yields only a 1.19-1.43-fold increase in effect.

Uncertainty propagation gave cenobamate first in 98% of draws under measurement uncertainty alone and
87% once unsourced parameters were included. A single unsourced constant - the
weight given to GABA-A potentiation relative to sodium-channel block - determines the leading drug, and
the ranking inverts at a GABA weight of 0.05, far below the value used (0.5), so the ordering is not finely balanced.

A structural limitation was identified and verified: for drugs acting only on excitability, the
predicted ranking is fixed by their pharmacology and is identical in every brain (observed threshold
ratio 3.36-3.40 across four connectomes against 3.35 predicted). Only the absolute dose required is
patient-specific (cenobamate 0.48 to 1.28 across four connectomes). Patient-specific ordering can arise only from drugs combining
mechanisms.

A pre-registered virtual-surgery test against real surgical outcomes did not support the hypothesis (n = 6; AUC 0.44, one-sided p = 0.65), and is reported as a negative result; per-patient calibration did remove the ceiling effect that invalidated an earlier pilot. The available sample is under-powered: with three
patients per group the smallest attainable one-sided p-value is 0.05, so the earlier pilot had zero
statistical power by construction.

## Conclusion

ARIADNE connects molecular pharmacology to network-level seizure dynamics on real brain wiring, and
makes every assumption behind its output measurable rather than implicit. It is not validated for
clinical drug selection, and cannot currently be: no public dataset links a brain scan to which drug a
patient responded to. The predictions that would test it, and the data required, are pre-registered.
The framework is presented as hypothesis-generating, models efficacy only, and does not account for
tolerability, which drives a substantial share of real prescribing decisions.

**Keywords:** epilepsy, connectome, antiseizure medication, computational pharmacology, uncertainty
quantification, pre-registration
