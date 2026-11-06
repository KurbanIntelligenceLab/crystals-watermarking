# Fixed-budget detection after deformation

BasinMark compares 32 selected structures from 256 candidates in each trial. A frozen audit prioritizes fingerprints likely to survive downstream processing. On two independent CrystalFormer groups it detected 385/500 and 490/500 marked collections after bounded strain and re-relaxation. It retains every positive decision of the previous audit on the same deposited collection. The variation between these groups motivated a separate, optional producer stress preview. [Audit design](protocol.md) · [Producer comparison design](producer-protocol.md)

## Producer stress preview

The producer deforms a copy of each candidate, relaxes the copy and checks its fingerprints. It can use spare energy allowance to prefer candidates carrying the mark across the preview. It emits the original geometry. The comparison enforces the same per-rank and mean energy caps, candidate count, and number of green producer records. Individual emitted structures can differ. Scores on the eventual deposits are withheld from both selectors. The unpreviewed budget selector remains the default.

Each row counts detections in 500 trials using the same keys and candidate draws in both arms. Intervals vary keys and draws within the indicated fixed pool. E and F are the two prospective cohorts for this producer change. C and D were held out for the preceding auditor change and are development comparisons for the producer. The quality columns give the change in mean producer ranking energy in meV/atom and the fraction independently labeled stable for the old and preview producers. CrystalFormer uses source phonon labels. MatterGen counts DFT-confirmed stable emissions among all emissions. Worker-seconds per emitted structure include the preview relaxation and fingerprint measurement on all eight candidates on average, including failed attempts. Hardware, concurrency and timing affect these values. Preview and target use separate copies of the coordinates but closely related axial strains and the same model family.

Both prospective groups cleared the prespecified 450/500 target after strain. The budget producer already detected 492/500 on E and 493/500 on F. The preview recovered 4 and 4 additional collections, bringing the totals to 496/500 and 497/500. Ordinary detection was 498/500 in both arms on both groups. Earlier comparisons included small losses with the preview in ordinary relaxation and shear. These measurements support offering a preview for axial strain when its extra cost is acceptable.

| Cohort | Treatment | Budget | Preview | 95% interval, preview | Gained / lost | Preview wrong keys | Δ energy | Stable, budget / preview | Extra worker-s/output |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| E | strain rerelaxed | 492 | 496 | 98.0–99.8% | 4 / 0 | 16/10,000 | +2.9 | 22.7% / 22.7% | 17.1 |
| E | rerelaxed | 498 | 498 | 98.6–100.0% | 0 / 0 | 13/10,000 | +2.9 | 22.7% / 22.7% | 17.1 |
| F | strain rerelaxed | 493 | 497 | 98.3–99.9% | 4 / 0 | 16/10,000 | +4.3 | 25.4% / 25.2% | 15.0 |
| F | rerelaxed | 498 | 498 | 98.6–100.0% | 0 / 0 | 14/10,000 | +4.3 | 25.4% / 25.2% | 15.0 |
| C | strain rerelaxed | 385 | 457 | 88.6–93.7% | 72 / 0 | — | +5.5 | 31.6% / 30.4% | 30.9 |
| C | rerelaxed | 487 | 494 | 97.4–99.6% | 7 / 0 | — | +5.5 | 31.6% / 30.4% | 30.9 |
| D | strain rerelaxed | 490 | 492 | 96.9–99.3% | 3 / 1 | — | +2.8 | 30.1% / 29.6% | 20.6 |
| D | rerelaxed | 495 | 496 | 98.0–99.8% | 1 / 0 | — | +2.8 | 30.1% / 29.6% | 20.6 |
| A | strain rerelaxed | 454 | 472 | 92.0–96.2% | 18 / 0 | — | +4.4 | 36.6% / 35.5% | 25.3 |
| A | rerelaxed | 494 | 492 | 96.9–99.3% | 0 / 2 | — | +4.4 | 36.6% / 35.5% | 25.3 |
| B | strain rerelaxed | 463 | 487 | 95.6–98.6% | 24 / 0 | — | +3.2 | 27.7% / 27.3% | 38.2 |
| B | rerelaxed | 489 | 489 | 96.1–98.9% | 1 / 1 | — | +3.2 | 27.7% / 27.3% | 38.2 |
| dev | strain rerelaxed | 488 | 493 | 97.1–99.4% | 5 / 0 | — | +5.6 | 31.9% / 31.4% | 14.4 |
| dev | rerelaxed | 498 | 497 | 98.3–99.9% | 0 / 1 | — | +5.6 | 31.9% / 31.4% | 14.4 |
| first-confirm | strain rerelaxed | 451 | 466 | 90.6–95.2% | 19 / 4 | — | +3.5 | 34.7% / 34.7% | 15.2 |
| first-confirm | rerelaxed | 489 | 493 | 97.1–99.4% | 4 / 0 | — | +3.5 | 34.7% / 34.7% | 15.2 |
| shear | shear rerelaxed | 470 | 469 | 91.3–95.7% | 6 / 7 | — | +3.5 | 34.7% / 34.7% | 15.2 |
| mattergen | rerelaxed | 500 | 500 | 99.3–100.0% | 0 / 0 | — | +0.9 | 82.2% / 82.3% | 4.0 |
| mattergen | paired dft | 495 | 497 | 98.3–99.9% | 2 / 0 | — | +0.9 | 82.2% / 82.3% | 4.0 |

The prospective strain tests each gave 16 positive audits in the original 10,000 wrong-key trials (0.16%). We checked this with 200 additional independently drawn keys on each of the same 500 preview releases, keeping those first results in the table. The follow-up gave 85/100,000 on E and 93/100,000 on F at the nominal 0.1% threshold. The conditional audit sizes predict 99.9 and 99.9 positive decisions. This follow-up measures variation over fresh keys on these fixed releases. [Per-release control records](preview-controls.json.gz)

## Bounded removal with audit feedback

A third party tries up to eight ordinary or edited collections, with access to the binary audit result. The collection membership stays fixed. The ordinary audit detected 489/500. Applying all energy-eligible shears detected 483/500. Choosing between the proposals gave 474/500, a loss of 15 previously detected collections. A separate force calculation compared both endpoints of each of the 85 distinct emitted structures eligible for an edit. 11 passed the phonon screen at both endpoints, 1 lost a pass, 3 gained one and 70 failed both. Restricting edits to pairs passing at both endpoints gave 487/500. Allowing an edit whenever it preserved an existing pass gave 474/500. These tests assess one shear, one cohort and a finite model-based phonon screen. The reference and deformed geometries both come from the initially released coordinates.

## Cross-catalogue transfer at lower exposure

With 256 genuine outputs available to the observer, cross-catalogue reuse triggered 20/200 full-fingerprint audits and 105/200 no-cell audits. At the previous exposure of 1,130 genuine outputs, both representations gave 200/200. The lower-exposure comparison uses the same catalogues and paired seeds and changes the observation budget. [Earlier catalogue comparison](../v0.3/README.md)

## Running the comparison

Use the replay scripts to regenerate every statistical record from archived measurements. The producer API accepts an optional `preview` policy with one primary/secondary preview pair per candidate. A JSON input to `basinmark select-collection --policy preview` carries a `preview` pair or explicit null on each candidate. Each scored deposited audit uses the fixed survival table from four historical cohorts. [Reproduction instructions and source identities](reproduce.md)
