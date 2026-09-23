from common171b import *
base=R/'90_pre_lincs_compound_expansion'
paths=[base/'2026-08-14_stage01_systematic_expansion'/name for name in ['07_expanded_CHEM_AB_unique_compounds.csv','04_structure_standardization.csv','02_raw_compound_evidence.csv']]+[base/'2026-08-14_stage02_LINCS_coverage_audit/03_compound_to_LINCS_entity_mapping.csv',H/'RESULTS/PRIMARY_36_TEST_RESULTS.csv',H/'RESULTS/CHEMICAL_IDENTITY_AND_ALIAS_AUDIT.csv']
rows=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),registered_utc=now(),content_use='DEFERRED_UNTIL_BLINDED_STATISTICS_FROZEN; header inspection only before that gate') for p in paths]
csvwrite(O/'POST_STATISTICAL_ANNOTATION_SOURCE_REGISTRY.csv',rows)
for p in paths:
    with p.open(encoding='utf-8-sig',newline='') as f:header=next(csv.reader(f))
    print(p.name,header)
