# SVChordinator
merge, annotate, and visualize structural variants (as "chords" in a circos plot).

Currently supports the following callers:

- `SvABA`
- `manta`
- `GRIPSS`
- `nanomonsv`
- `SAVANA`
- `Severus`
- `cuteSV`
- `Sniffles2`

input tsv file with paths to vcf must follow format of the `config/example_samples.tsv` file,
in which callers are labeled with the names of the tools
