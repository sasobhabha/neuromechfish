Anatomical meshes for NeuroMechFish
===================================

`head.stl` and `fin.stl` are derived from the simZFish larval-zebrafish
model (Xiangxiao Liu, Luca Zunino, et al., "Artificial embodied circuits
uncover neural architectures of vertebrate visuomotor behaviors",
Science Robotics 10, eadv4408 (2025); EPFL Biorobotics Laboratory).

Source (Apache License 2.0):
    https://ponyo.epfl.ch/proj/zebrafish/simzfish/
    files:  Mesh/Zebrafish_head.stl  ->  head.stl
            Mesh/Tail fin.stl        ->  fin.stl

Import/transform script:  scripts/import_simzfish_meshes.py
Vendored originals in:    _thirdparty/simzfish/

Copyright (c) 2022-2025 the simZFish authors / EPFL BioRob.
Licensed under the Apache License, Version 2.0 (see
_thirdparty/simzfish/LICENSE). Meshes are re-oriented to the
NeuroMechFish frame (cyclic axis permutation) and rescaled to the
model's 4.15 mm larva; no shape data was otherwise altered.
