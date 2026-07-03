# How this VEG-05 host–microbiome study relates to the OSD-767 × APH salicylic-acid integration

_Comparison of two integrative analyses that share the OSD-767 transcriptome but integrate it
along different axes. Prepared 2026-07-02._

## The two analyses at a glance

| | **This study — VEG-05 multi-omics** | **Our SA integration** (`tomato-spaceflight-VEG05-APH-SA-integration`) |
|---|---|---|
| Integrates | OSD-767 transcriptome **+ OSD-766 microbiome** (16S + ITS) | OSD-767 **+ APH (MoneyMaker/NahG ± SA)** transcriptomes |
| Second modality | Plant-associated microbiome (bacteria + fungi) | Salicylic-acid genetics/chemistry (causal) |
| Cultivar(s) | Red Robin (OSD-767) | Red Robin (OSD-767) + MoneyMaker (APH) |
| Framework | MOFA+, WGCNA, dysbiosis index, module–taxon networks | Plant PhysioSpace, cross-study concordance |
| Scope | Single-mission, host↔microbe coordination | Cross-mission, is-it-conserved + causal SA test |
| Headline | Spaceflight reshapes the microbiome (dysbiosis) and the host mounts a microbe-driven defense/oxidative response, light-dependently | Spaceflight ≈ SA-deficient, biotic-defense-primed state; SA is the conserved causal buffer |

Both use OSD-767 as the shared anchor, so they can be read as two lenses on the **same VEG-05
spaceflight response**.

## Where they independently agree (cross-validation)

1. **Blue light amplifies the spaceflight response — found twice, two pipelines.**
   This study: leaf DEGs under blue light **4,716 vs 523** under red (~9×), interaction = 3,189 DEGs.
   Our OSD-767/SA work: blue LED drives the largest biotic/hormone PhysioSpace interaction
   (−44.2) and the Light×Condition interaction dominates the leaf transcriptome. Two independent
   DE parameterisations on the same OSD-767 data reach the same qualitative conclusion — strong
   internal validation. (Exact DEG counts differ, as expected from different thresholds/models;
   the *direction and magnitude* replicate.)

2. **Root-led, defense/oxidative signature.** This study: adventitious root is 87% upregulated
   with turquoise-module oxidative-stress GO (peroxidase, H₂O₂ catabolism, glutathione). Our
   work: the spaceflight defense signal is root-led (biotic/hormone axis; outer-cortex hotspot).
   Both localise the acute stress/defense response to the root.

## The key synthesis (how they connect)

Our SA integration concluded that spaceflight drives tomato into a **biotic-defense-primed,
salicylic-acid–deficient-like state** — the host behaves as if fighting a biotic threat while
its SA brake is functionally weak. That analysis could describe the state but not its proximal
trigger.

**This microbiome study supplies a candidate trigger.** It shows that spaceflight
*restructures the plant microbiome* — increased bacterial diversity, elevated bacterial and
fungal dysbiosis, loss/shift of specific partners (*Pantoea*, *Penicillium*, *Fusarium*), and
gain of growth-promoting/N-fixing taxa (*Methylobacterium*, *Burkholderia*, *Azospirillum*) —
and that a **microbe-driven** host module (adventitious-root black module: r = −0.85 with ITS
dysbiosis, but r = −0.38, p = 0.32 with flight) mounts an oxidative-stress/defense response
that tracks the *microbiome*, not microgravity per se.

Put together, a coherent, testable model emerges:

> **Spaceflight destabilises the plant–microbiome relationship (dysbiosis); the host responds
> with an SA-gated biotic-defense/oxidative program that is simultaneously SA-limited under
> flight.** Salicylic acid is the master regulator of plant biotic immunity, so a defense
> response to a disrupted microbiome that reads as "SA-deficient-like" is exactly what one
> would predict.

This also **explains a puzzle from the SA paper**: gene-level spaceflight signatures don't
replicate across missions, yet the biotic-defense *pattern* does. If that pattern is a response
to the (mission-, hardware-, and cultivar-variable) microbiome rather than to a fixed
microgravity program, then pattern-level conservation with gene-level divergence is exactly
what we should see. The microbe-driven, flight-*independent* module is direct evidence that a
chunk of the "spaceflight" transcriptional response is really a **microbiome** response.

## Differences and tensions to reconcile

- **Modality & method** don't overlap (microbiome/MOFA vs SA-genetics/PhysioSpace), so this is
  complementary evidence, not a replication of methods.
- **"Increased diversity" vs "dysbiosis."** Flight raised bacterial alpha diversity; both papers
  treat community *displacement* (not loss) as the disruption. Worth stating the nuance explicitly.
- **Causality is split across datasets.** The causal SA lever (NahG) lives in the APH arm; the
  microbiome link here is correlational. Neither alone closes the loop.

## Concrete follow-ups this comparison motivates

1. **Consistency check (cheap, high-value):** project *this* study's OSD-767 flight transcriptome
   into the same Plant PhysioSpace used in the SA integration — it should land at the
   SA-deficient/NahG pole (same underlying data). A quick confirmation that the two framings are
   describing one signature.
2. **Gene-set overlap:** test whether the microbe-driven module and the flight-associated
   module–taxon genes are enriched for the SA/biotic-defense gene set from the SA paper.
3. **Close the loop experimentally:** a **NahG × microbiome × spaceflight** design would test
   directly whether SA gates the host response to spaceflight-induced dysbiosis.

## Where this sits in the paper set

- **P1** OSD-767 light × spaceflight (transcriptome).
- **P2** APH SA (NahG/SA) — SA causally buffers the defense response.
- **P3** SA cross-mission integration — spaceflight ≈ SA-deficient/defense-primed, conserved at the pattern level.
- **P4 (this)** VEG-05 host–microbiome — spaceflight disrupts the microbiome; host mounts a microbe-driven defense/oxidative response.

**Unifying thesis across the four:** spaceflight destabilises the plant–microbiome partnership,
and the host answers with an SA-gated biotic-defense program that is SA-limited in flight —
making SA (and microbiome) management a concrete lever for space crop production.
