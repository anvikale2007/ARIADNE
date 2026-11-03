# ARIADNE at Khoj, 3 November 2026 - what to say and how to defend it

All numbers below are from the completed analysis of 31 October 2026.

---

## The claim, reframed

Do not say "ARIADNE ranks antiseizure drugs." It overclaims, and a sharp judge will take it apart.
Say this instead:

> **ARIADNE is a mechanism-to-network translation framework.** Given published drug potencies and a
> patient's own brain wiring, it works out how much a drug's measured molecular effect must change
> seizure spread across the whole network - and it states exactly how much of that answer is
> measurement and how much is assumption.

Every criticism then becomes a named parameter you have already quantified, instead of a hole someone
else finds.

---

## The 60-second pitch

> About a third of people with epilepsy keep having seizures despite medication, and choosing the right
> drug is still largely trial and error, one drug at a time over months.
>
> Drug potencies are measured on single ion channels. Seizures are a network problem: they start in one
> place and spread along white-matter tracts. Nobody routinely connects those two scales.
>
> ARIADNE does. It takes a diffusion MRI scan, builds the brain's wiring map, runs a seizure model in
> every region, and asks for each drug: how much of it does it take to stop the seizure spreading in
> *this* brain?
>
> It runs on real epilepsy patients' connectomes, every drug number is traced to a published source,
> and - this is the part I care about - I have measured how much the answer depends on the numbers I
> could not source. That is on the poster, not hidden in an appendix.
>
> It is not a clinical tool and it is not validated for choosing a drug, because no public dataset links
> a brain scan to which drug worked. What I can test, I pre-registered and tested.

---

## Booth order (what to click, in order)

1. **Seizure type** - classify first, because some drugs make certain seizure types worse. Switch to
   *Myoclonic* and show carbamazepine and pregabalin turning red. This shows clinical literacy fast.
2. **Brain wiring → Seizure** - press play, let them watch it spread. This is the hook.
3. **Add a drug** - cenobamate contains it, carbamazepine does not. Then switch the genetic test to
   HLA-B\*15:02 and watch carbamazepine get excluded.
4. **Ranking** - then immediately go to the uncertainty figures. Do not let the bar chart be the last
   word; the probability figure is the grown-up version of it.
5. **Real patients / Evidence** - only if they are still engaged.

Keep step 2-3 under 90 seconds. Most judges leave within three minutes.

---

## The six questions, with answers

Say these out loud until they are automatic. Short answers. Do not ramble.

**1. "What is your GABA weight based on?"**
> Nothing - it is unsourced, and it decides the winner. So I measured how much that matters rather
> than hiding it. The ranking only flips below a GABA weight of 0.05, and we use 0.5 Below that value lamotrigine ranks first; above it, cenobamate.
> The experiment that would settle it is a patch-clamp study measuring GABA-A potentiation against
> sodium-channel block in the same preparation.

**2. "Your virtual surgery AUC was 0.33. That is worse than chance."**
> It was, and the pilot could not have succeeded. With three patients per group the smallest
> achievable one-sided p-value is exactly 0.05, so its statistical power was zero - it was
> uninformative by construction. I then found two defects using pre-surgery diagnostics only: a
> ceiling effect, where every network was already saturated, and partial resections being
> under-counted. I fixed both, wrote the analysis plan down before re-running it, and the result is
> AUC 0.44 - still not supported, and I report it as a negative result. The data needed to do it properly is 40 files from a public dataset.

**3. "Why is levetiracetam missing? It is the most prescribed drug."**  *(It is no longer missing — use this:)*
> It's in, as of today, and it breaks my model. SV2A, its main mechanism, has no published effect size, so I
> modelled the one target that does: it blocks N-type calcium channels by about 18%, and transmitter release
> scales with calcium to the power 2.6. That gives a derived 40% reduction in release - no invented numbers.
> Modelled that way it comes **last of six**, and loses to lamotrigine in 95% of parameter draws. But NICE lists
> lamotrigine and levetiracetam side by side as the two first choices for focal epilepsy. So my model is wrong
> about a first-line drug, and I can say exactly how: either the coupling term is far too weak, or SV2A does
> nearly all the work. Both are testable, and that is the next experiment.

**3b. (old wording, if asked why you did not just invent a number)**
> Because no published number converts SV2A binding into a size of effect, and I will not invent one.
> What I can do is invert the question: a drug acting purely by reducing transmitter release would
> need a 175% reduction in transmitter release to match lamotrigine here - which is impossible. So
> either my coupling term is too weak to carry a real drug effect, or levetiracetam works through
> something I'm not modelling. That's a concrete failure of my model, and it's the next thing I fix.

**4. "You get the same drug order in every brain. So what is patient-specific?"**
> I found that too, and it is a real limitation. For drugs that act only on excitability, the ranking
> is fixed by their pharmacology and comes out the same in every brain - I verified the ratio is 3.35
> across four brains. Only the absolute dose is patient-specific: cenobamate needs 0.47x in one brain
> and 1.28x in another. Order only becomes patient-specific for drugs mixing both mechanisms. That
> tells me exactly what the next version has to model.

**5. "Can a doctor use this?"**
> No, and I would not want one to. It is not validated for drug choice, and it models efficacy only -
> it ignores side effects, which drive a large share of real prescribing decisions. It is a
> hypothesis generator, and the hypothesis is pre-registered.

**5b. "Hasn't this been done? Newcastle have been predicting surgical outcomes from structural networks for years."**
> Yes — Sinha's 2021 thesis under Peter Taylor, in the same patient cohort whose data I'm using. I'm not
> claiming the virtual-surgery test as novel; I borrowed it as a way to check whether my network model carries
> real information, and it currently fails at n=6. What none of that work does is model drugs. That's the layer
> I built: published drug potencies translated into predicted effects on seizure spread. Their method also
> tells me what to fix — they standardise each patient against healthy controls, and I don't, which is probably
> why my measure tracks how much tissue was removed rather than which tissue.

**6. "How is this different from EpiPick?"**
> EpiPick is validated in 425 patients and picks drugs from the clinical profile - seizure type, age,
> comorbidities. That is a different layer. ARIADNE asks, among the drugs that are already safe and
> appropriate, which one best fits this patient's seizure network. I built the clinical filter in
> front of the simulation, so they compose rather than compete.

### If you are asked something you cannot answer
Say: *"I don't know - that's not something I've tested."* Then say what you would do to find out.
That answer scores better than a guess, every time.

---

## The limitations slide (put this on the poster; do not bury it)

- Not validated for choosing a drug. No public dataset links a brain scan to drug response.
- Efficacy only. Tolerability is not modelled, and it drives much real prescribing.
- Some effect sizes are assumed, not measured: the GABA weight, the transmitter-release reductions,
  valproate's lumped GABA/T-type term. Each is quantified in the uncertainty analysis.
- Ranking among pure-excitability drugs is fixed by pharmacology, not by the patient's brain.
- Virtual surgery is a pilot and is currently under-powered.
- Guideline rules are UK (NICE); Indian practice and drug availability differ. Cenobamate had only a
  CDSCO expert-committee recommendation as of June 2026.

---

## Two things not to do

1. **Do not call it diagnostic, clinical, or validated.** One such word and a clinician judge stops
   believing the rest.
2. **Do not build on Monday.** Monday is rehearsal, printing and sleep. A project presented well beats
   a better project presented by someone who was awake until 4am.
