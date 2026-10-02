# Combined Pangenome Pipeline (RIBAP / Panaroo / PPanGGOLiN)

A single wrapper script, `run_pangenome.sh`, that takes a folder of genome
FASTAs and routes them through **RIBAP**, **Panaroo**, **PPanGGOLiN**, or all
three, via one `--tool` flag.

## 1. Install dependencies

```bash
# Create one environment with everything needed
conda create -n pangenome -c bioconda -c conda-forge \
    prokka nextflow panaroo ppanggolin -y
conda activate pangenome
```

RIBAP itself doesn't need a separate install — `nextflow run hoelzer-lab/ribap`
pulls the pipeline from GitHub automatically the first time you run it.

## 2. Prepare your input

Put one genome assembly per FASTA file (extension `.fasta`, `.fa`, or `.fna`)
into a single directory:

```
genomes/
├── strainA.fasta
├── strainB.fasta
└── strainC.fasta
```

## 3. Pick a tool and run

```bash
chmod +x run_pangenome.sh

# Just RIBAP
./run_pangenome.sh -i genomes/ -o results/ -t ribap

# Just Panaroo
./run_pangenome.sh -i genomes/ -o results/ -t panaroo

# Just PPanGGOLiN
./run_pangenome.sh -i genomes/ -o results/ -t ppanggolin

# Run all three and compare
./run_pangenome.sh -i genomes/ -o results/ -t all -c 16 -s 0.8
```

## 4. Which tool should I pick?

| Your situation | Recommended tool |
|---|---|
| Comparing strains within **one species** | `panaroo` (best noise-cleaning) or `ppanggolin` (fastest, scales well, gives persistent/shell/cloud partitions) |
| Comparing genomes across **multiple species / a genus** | `ribap` — in the original benchmark it recovered ~60% of annotated genes as core at genus level vs. 3–30% for the others at default settings |
| Not sure yet / want to sanity-check | `all` — runs everything once, so you can compare core genome sizes directly |

If you run `panaroo` or `ppanggolin` on cross-species data and see the core
genome collapse toward zero as you add more diverse genomes, that's the
signature RIBAP's ILP refinement step was built to fix — switch to `ribap`
for that dataset rather than continuing to loosen the similarity threshold.

## 5. Output locations

```
results/
├── prokka/<sample>/<sample>.gff      # only created for panaroo/ppanggolin/all
├── ribap/...                          # RIBAP's own Nextflow output structure
├── panaroo/
│   ├── gene_presence_absence.csv
│   └── summary_statistics.txt
├── ppanggolin/run/
│   └── genomes_statistics.tsv
└── logs/                              # per-tool run logs
```

## 6. Tuning notes

- `-s/--similarity` controls the sequence identity threshold for Panaroo and
  PPanGGOLiN. Defaults (~95%) suit single-species comparisons; drop to
  0.7–0.8 for more divergent, cross-species inputs — though per the RIBAP
  benchmark, loosening this only partially closes the gap versus RIBAP's ILP
  approach at genus level.
- `-p/--profile` switches RIBAP between `local` execution and a `slurm` HPC
  profile if you have one configured — useful since RIBAP's ILP step is the
  most compute-intensive of the three.
- Check `nextflow run hoelzer-lab/ribap --help` for the exact input-flag name
  in your installed RIBAP version (it has changed across releases — samplesheet
  vs. glob pattern) and adjust the `run_ribap` function in the script if needed.

## Reference

Lamkiewicz K, Barf LM, Sachse K, Hölzer M. *RIBAP: a comprehensive bacterial
core genome annotation pipeline for pangenome calculation beyond the species
level.* Genome Biology 25, 170 (2024).
https://github.com/hoelzer-lab/ribap
