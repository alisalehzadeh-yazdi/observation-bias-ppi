# Observation Bias in Evolving Protein Interaction Networks

Code accompanying:

> Salehzadeh-Yazdi A, Tuncbag N, Gursoy A, Keskin O, Hütt M-T. *A
> quantitative analysis of observation biases in protein interaction
> networks.*

This repository quantifies whether successive releases of public
protein-protein interaction (PPI) databases (STRING, BioGRID, HIPPIE,
IntAct) expand by uniformly discovering unknown biology, or by
preferentially re-visiting already well-connected proteins. It does this by
comparing each empirically observed database update against edge-addition
null/alternative models, using a Jaccard-overlap-based, cross-database
normalized bias score.