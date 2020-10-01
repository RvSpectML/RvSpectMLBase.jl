```@meta
CurrentModule = RvSpectMLBase
```

# RvSpectMLBase

## Getting Started

- [Install Julia 1.5](https://julialang.org/downloads/).  On Penn State's ICS-ACI, it is avaliable at  `/gpfs/group/ebf11/default/julia/bin/julia`.
- Install the [RvSpectMLBase package](https://github.com/RvSpectML/RvSpectMLBase.jl) and it's dependencies.  From julia
```julia
import Pkg
Pkg.add("https://github.com/RvSpectML/RvSpectMLBase.jl")
Pkg.instantiate()
```
- Run the tests
```
> include("test/runtests.jl")
```


## Related Documentation
- [RvSpectMLBase](https://rvspectml.github.io/RvSpectMLBase.jl/stable/)
- [EchelleInstruments](https://rvspectml.github.io/EchelleInstruments.jl/stable/)
- [EchelleCCFs](https://rvspectml.github.io/EchelleCCFs.jl/stable)
- [Scalpels](https://rvspectml.github.io/Scalpels.jl/stable/)


