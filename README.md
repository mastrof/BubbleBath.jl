# BubbleBath.jl

[![Build Status](https://github.com/mastrof/BubbleBath.jl/workflows/CI/badge.svg)](https://github.com/mastrof/BubbleBath.jl/actions)
[![codecov](https://codecov.io/gh/mastrof/BubbleBath.jl/branch/main/graphs/badge.svg)](https://codecov.io/gh/mastrof/BubbleBath.jl)

Generate loose packings of spheres in orthorhombic domains, in 2 and 3 dimensions.

## Features
* Fill a domain with spheres from a given distribution of radii to reach a target
    packing fraction, or from already-sampled radii.
* Control minimum allowed distance between spheres.
* Decide whether spheres can cross through domain boundaries or not.