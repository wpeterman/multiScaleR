# Compact summary consumer

analyze.py reads only the four declared -v2 result files and frozen study
identities. It performs no fitting, sampling, raw-chain reconstruction or pooling.
Its final writer requires the complete planned panel and unchanged identities;
missing outcomes remain pending. It never retries a producer or changes a gate.

Each run must retain all original raw, derived, fixed/own-threshold functional,
source and producer gates. Both independent runs must match the declared case,
input hash, settings and distinct seeds. Both complete unique sets of 228 fixed
events must match exactly. All probabilities must agree within the predeclared
three combined MCSEs. Retain every disagreement, unavailable constant indicator,
failed precision decision, missing record and original pre-sampling failure.

Compact R summaries use twelve decimal places. If serialization uncertainty
straddles the existing three-combined-MCSE boundary, report that event unavailable
for the agreement decision, rather than round it into a pass. This handles
finite stored precision; it does not substitute a new statistical cutoff. The
old Laplace comparisons retain their separate historical-rounding qualification.

test-analysis.py exercises the independent consumer with tiny deterministic
records, including each rejected gate, mismatched targets/seeds/event sets,
constant/nonfinite indicators, disagreement and serialization-boundary cases.
These are software tests, not additional sampler fits or scientific datasets.
