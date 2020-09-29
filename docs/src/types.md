```@meta
CurrentModule = RvSpectMLBase
```
# Types Exported by RvSpectMLBase

```@contents
Pages = ["types.md"]
Depth = 3
```
## Abstract Types
```@autodocs
Modules = [ RvSpectMLBase ]
Private = false
Order = [:type]
Filter = t -> isabstracttype(t)
```

## General purpose
```@autodocs
Modules = [ RvSpectMLBase ]
Private = false
Order = [:type]
Filter = t -> !isabstracttype(t)
```

## Instrument specific
```@autodocs
Modules = [ RvSpectMLBase.TheoreticalInstrument  ]
Private = false
Order = [:type]
```

## Other
```@autodocs
Modules = [ RvSpectMLBase.Pipeline  ]
Private = false
Order = [:type]
```
