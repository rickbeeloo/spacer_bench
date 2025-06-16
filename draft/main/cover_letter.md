Dear Editor,

Our manuscript titled "Tool Choice drastically Impacts CRISPR Spacer-Protospacer Detection" presents a comprehensive evaluation of sequence search tools commonly used in CRISPR-based host-virus interaction studies, which has important implications for the field of viral ecology.

Our study addresses a critical gap in the field: while CRISPR spacer-protospacer matching is widely used to infer host-virus relationships, the effects of sequence alignment or search tools remain unknown. Here, these tools are systematically evaluated, using both synthetic and real datasets (IMG/VR4). We demonstrate that different tools and parameters exhibit varying abilities to detect multiple matches and handle sequence variations, potentially affecting downstream biological interpretations.

Key findings of our study include:
1. No single tool achieves perfect recall within the range of valid matches (i.e. spacer-protospacer pairs).
2. Tools show marked differences in their ability to handle multiple occurrences of the same spacer.
3. Bowtie1 demonstrates superior performance for alignments with up to 3 mismatches.
4. The common and popular approach, using BLASTn-short, may miss significant numbers of valid matches.

Our work provides practical guidelines for researchers in the field and highlights important considerations for the interpretation of spacer-protospacer matches. This is particularly relevant given the rapid growth of viral and CRISPR spacer databases, where the choice of search strategy can significantly impact our understanding of host-virus relationships.  

The manuscript includes detailed methods, comprehensive benchmarking data, and discussion of biological implications. All code and data are publicly available through our git repository and Zenodo (https://zenodo.org/doi/10.5281/zenodo.15171878).

We believe this work will be of broad interest to researchers in bioinformatics, specifically those working on comparative sequence analysis, and developing tools for microbial ecology, virology, and host-virus interactions. The findings have immediate practical applications for improving host prediction pipelines.

Thank you for considering our manuscript.
Uri Neri, on behalf of the authors.
DOE Joint Genome Institute
Berkeley, CA, USA
uneri@lbl.gov 