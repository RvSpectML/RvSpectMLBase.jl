```@meta
CurrentModule = RvSpectMLBase
```
# RvSpectML Internals

As a heads up, these functions and types are more likely to change going forward than functions and types that are exported.  

```@contents
Pages = ["internals.md"]
Depth = 3
```
## Functions

### General purpose
```@autodocs
Modules = [RvSpectMLBase ]
Public = false
Order = [ :function ]
```


### Instrument specific
```@autodocs
Modules = [RvSpectMLBase.TheoreticalInstrument  ]
Public = false
Order = [ :function]
```

### Other
```@autodocs
Modules = [RvSpectMLBase.Pipeline  ]
Public = false
Order = [:function]
```

## Types

### General purpose
```@autodocs
Modules = [RvSpectMLBase ]
Public = false
Order = [:type ]
```

### Instrument specific
```@autodocs
Modules = [ RvSpectMLBase.TheoreticalInstrument  ]
Public = false
Order = [:type ]
```
### Other
```@autodocs
Modules = [RvSpectMLBase.Pipeline  ]
Public = false
Order = [:type]
```
