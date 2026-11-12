# Public Database Sequence Search & Hit Interpretation

## 1. Recorded Search Parameters
- **Query Sequence**: `query.fasta` (300 bp nucleotide sequence)
- **Search Tool**: NCBI `blastn` (Nucleotide-Nucleotide BLAST)
- **Target Database**: NCBI Nucleotide collection (`nt`)
- **E-value Cutoff**: `1e-5`
- **Scoring Matrix**: Match/Mismatch `2, -3`; Gap Existence `5`, Extension `2`

## 2. Interpretation of Hits
1. **Gene Identity**: The query shows a **100% identity and 100% query coverage** alignment (`E-value = 2e-154`, `Bit Score = 555`) to *Homo sapiens* hemoglobin subunit beta (`HBB`, accession `NM_000518.5`), confirming the query originates from the human beta-globin coding region (exons 1–2).
2. **Statistical Confidence**: An E-value of `2e-154` indicates virtually zero probability that this alignment occurred by random chance.
3. **Evolutionary Conservation**: Hits 3–5 represent orthologous `HBB` transcripts in closely related primates (*Pan troglodytes* at 99.67%, *Gorilla gorilla* at 99.33%, and *Macaca mulatta* at 97.67%), demonstrating strong purifying selection on the oxygen-binding globin domain across primates.

## 3. Repository Files
- `blast_search.py` — Script that queries the NCBI BLAST API, records parameters, and generates the annotation table.
- `query.fasta` — Input query sequence.
- `blast_output.txt` — Raw alignment summary and recorded search parameters.
- `annotation_table.md` — Structured annotation table of the top 5 hits.
