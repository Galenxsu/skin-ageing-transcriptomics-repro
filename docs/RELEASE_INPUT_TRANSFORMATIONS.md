# Public-release input transformations

The Phase 171A input directory used by the author-controlled isolated reproduction contained 55 audit and planning files. The public repository includes only the six machine-readable files required by the Phase 171B scientific pipeline.

Five files were copied byte-for-byte. In `pathway_signature_membership_preview.csv`, the `source` field was changed from an author-workstation absolute path to the corresponding GMT filename (`h.all.v2025.1.Hs.symbols.gmt`, `c2.cp.v2025.1.Hs.symbols.gmt`, or `c5.go.bp.v2025.1.Hs.symbols.gmt`). No query membership, identifier, direction, weight, order or other scientific field was changed.

`data/frozen_phase171a/FREEZE_MANIFEST.csv` is therefore a **public-release manifest**, not a copy of the original 55-file internal freeze manifest. Its hashes refer to the released files after path-only sanitisation.

Two compact Phase 170B result tables are included solely for the prespecified historical comparison performed during final export. They are not inputs to scoring, null generation or classification.
