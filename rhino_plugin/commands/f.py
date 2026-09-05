import os, compas_singular
from compas_singular.algorithms import SkeletonDecomposition
print(os.path.dirname(compas_singular.__file__))
print(hasattr(SkeletonDecomposition, 'real_neighbors'))

result = surface.discrete_mapping(D, crv_guids=crv_guids, pt_guids=pt_guids)
outer, inners, features, points = result
for f in features:
    print('feature ends:', f[0], f[-1])
print('outer sample near them:', min(outer, key=lambda p: (p[0]-f[0][0])**2 + (p[1]-f[0][1])**2))