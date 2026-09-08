# Activity grouping and dominant-system plotting: next integration choices

8 September 2026

This note reviews Philipp's [PR #2](https://github.com/Tijmenvermeij/SSLIP/pull/2), at `222a9de`, against the fitted results on the SSLIP review branch. Ordinary activity plotting has now been separated using his `plotSSLIP_SlipActivities` function/interface. Grouping and dominant-system ranking are not enabled by this checkpoint.

## Grouping needs a defined physical quantity

Let `q_j` be the four-component projected tensor actually used in fitting: the physical tensor, or its normalized version when `normalizeInplane = 1`. The net contribution of a group is:

`H_group = sum(q_j * gamma_j)`

If its tensors are exactly proportional, `q_j = c_j * q_reference`, an equivalent signed coefficient exists:

`gamma_group = sum(c_j * gamma_j)`

The PR accepts approximately equal or opposite projected tensors but then adds their coefficients unchanged ([comparison](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_SlipActivities.m#L133), [sum](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_SlipActivities.m#L165)). For opposite tensors and coefficients `+0.2, -0.2`, the raw sum is zero while the equivalent contribution is `0.4` times the first tensor.

Conversely, equal tensors with those opposite coefficients have zero net deformation but a sum of absolute coefficients of `0.4`. These quantities answer different questions. Absolute coefficient activity also depends on basis normalization and how a degenerate fit distributes activity.

**Recommendation:** begin with exactly indistinguishable projected tensors, accounting for signs and proportionality. Defer approximate grouping: near-equal tensors generally have no exact equivalent coefficient. If such groups are shown later, plot their net gradient contribution explicitly rather than labeling it as the amplitude of a representative physical system.

The PR's crystallographic gate also calls `angle` on Miller objects without disabling symmetry. This does not guarantee parallelism of the particular oriented `b,n` pair; the projected-tensor check is still needed. Rotation is already stored separately in this branch, so no rotation pseudo-systems should enter activity groups.

## Ranking should retain narrow bands

The PR ranks finite absolute coefficients using the [85th percentile](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L94), including valid zero pixels, and applies a threshold of `0.01`. A field with amplitude `0.2` on 10% of pixels has P85 zero. A diffuse field of `0.015` ranks above it, even though the narrow field has greater mean contribution (`0.02` versus `0.015` for unit-length projected tensors).

For an overview of the fitted deformation, a useful alternative is:

`score_j = mean(norm(q_j) * abs(gamma_j))`

This combines strength and occupied area and retains valid zeros. For a fixed decomposition, rescaling each basis tensor and inversely rescaling its coefficient leaves the score unchanged. Rerunning SSLIP with a different basis normalization can still change the decomposition because it changes the minimization objective. Invalid/failed fits are excluded. The score ranks contributions within the fitted decomposition; it does not establish that individual slip systems are uniquely identified.

The following calculations use the saved MTEX 7 Ni/HCP fits and their actual basis tensors. They do not rerun or alter identification:

| Example | Rotation | Top 3 by P85 | Top 3 by mean projected contribution | Systems with P85 ≥ 0.01 |
| --- | --- | --- | --- | --- |
| Ni | off | 1, 8, 11 | 8, 1, 11 | 2/12 |
| Virtual HCP | off | 3, 11, 4 | 3, 7, 19 | 0/24 |
| Ni | on | 1, 8, 11 | 8, 1, 11 | 2/12 |
| Virtual HCP | on | 3, 7, 8 | 3, 7, 19 | 0/24 |

The virtual HCP input was generated with systems `3, 7, 19`. All 24 fitted systems fall below the PR's default P85 threshold, both with and without rotation. The mean-contribution ranking places those three generated systems first in both cases. This supports using mean projected contribution for the initial dominant-system overview, while leaving the full activity maps available.

A percentile over active pixels or a peak value would instead emphasize local band strength and would need an explicit definition of an active pixel. The PR's current activity threshold should not be transferred unchanged to a different metric.

## Selected systems and omitted activity

The dominant plot uses coefficient-channel indices to index the full slip-system list ([ranking](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L95), [labels and traces](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L326)). The adaptation must map each activity row through `opt.NoSs`; otherwise subsets and reordered lists are mislabeled.

Its [leftover map](https://github.com/PhilKro/SSLIP/blob/222a9de2448eb7e152a422d1b3bf0d335fd2bb06/src/plotting/plotSSLIP_DominantActivities.m#L213) adds signed coefficients from potentially unrelated tensors. Replace that with an explicitly named quantity:

- Net omitted gradient magnitude: `norm(sum(q_j * gamma_j))`, including physical cancellation.
- Sum of omitted contribution magnitudes: `sum(norm(q_j) * abs(gamma_j))`, measuring fitted activity without cancellation.

These should not share an unqualified raw-slip-amplitude colorbar with the top-system maps.

**Suggested next implementation:** adapt Philipp's dominant plot with correct original labels, rank by mean projected contribution, and use a consistently labeled contribution magnitude for the omitted systems. Keep automatic grouping disabled until its exact-equivalence behavior is separately checked. Preserve his attribution and use the existing export functions rather than the PR's missing personal export helper.
