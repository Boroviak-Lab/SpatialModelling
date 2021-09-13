# Spatial modelling

Demo scripts for investigating marmoset spatial gene expression. Files with names demo_CSX.m are demo scripts for loading in and visualising date for Carnegie Stage X. We have Four such stages:

CS5, CS6, CS7 and earlyCS6

Within these files there are demos for loading in the spatial data, using GP models for smoothing expression patterns. Expression can then be used to:

1) 3D reconstruction of expression patterns from projected samples in 3D
2) 2D reconstruction of tissues at chosen plane bisections
3) 1D reconstructions along a given trajectory (for extracting AP gradients or for generating virtual 2D gastruloids)

Other demo code is availble:

3DPlots_withAP_EmDisc.m and 3DPlots_withAP_VE.m are used to project expression patterns on to 3D embryos with a focus on the EmDisc and VE respectively. SectionPlots.m generates side by side visualisation of sections with expression projections.

In vitro model mapping:

MapPrimedNiaveCorrelations1.m and MapPrimedNiaveCorrelations2.m together project primed and niave PSCs on to the marmoset embryo (CS5 and CS6).

MapAmnioidCorrelations.m maps PSCs and amnion-like cells from human amnioids (10X) data on to the marmsoset embryo.

