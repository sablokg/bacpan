#!/usr/bin/env bash
#
# run_pangenome.sh — Combined bacterial pangenome pipeline
# Wraps Prokka annotation, then routes to RIBAP, Panaroo, and/or PPanGGOLiN.
#
# Gaurav Sablok gsablok@proton.me
# ---------------------------------------------------------------------------
# WHY THIS SCRIPT EXISTS
# ---------------------------------------------------------------------------
# RIBAP, Panaroo, and PPanGGOLiN all consume the same starting point (a folder
# of assembled genome FASTAs) but expect different downstream inputs:
#   - RIBAP:      raw genome FASTAs -> it runs Prokka + Roary + ILP internally
#                 via its own Nextflow pipeline.
#   - Panaroo:    pre-annotated GFF3 files (Prokka output) as input.
#   - PPanGGOLiN: pre-annotated GFF3 files (Prokka output) as input, listed
#                 in a tab-separated "organisms" file.
#
# This script annotates once with Prokka (needed for Panaroo/PPanGGOLiN),
# and lets RIBAP do its own annotation internally (recommended, since RIBAP's
# Nextflow pipeline manages Prokka versions/params itself). You choose which
# tool(s) to run with -t / --tool.
#
# ---------------------------------------------------------------------------
# REQUIREMENTS (install ahead of time; see README.md for one-liners)
# ---------------------------------------------------------------------------
#   - prokka       (genome annotation)          conda install -c bioconda prokka
#   - nextflow     (to run RIBAP)                conda install -c bioconda nextflow
#   - RIBAP        (github.com/hoelzer-lab/ribap) - no install needed, nextflow
#                  pulls it directly via -r/--repo below
#   - panaroo      conda install -c bioconda panaroo
#   - ppanggolin   conda install -c bioconda ppanggolin
#
# ---------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------
#   ./run_pangenome.sh -i /path/to/genomes_fasta_dir -o /path/to/output -t ribap
#   ./run_pangenome.sh -i /path/to/genomes_fasta_dir -o /path/to/output -t panaroo
#   ./run_pangenome.sh -i /path/to/genomes_fasta_dir -o /path/to/output -t ppanggolin
#   ./run_pangenome.sh -i /path/to/genomes_fasta_dir -o /path/to/output -t all
#
# Input directory must contain one FASTA file per genome (extensions:
# .fasta, .fa, .fna — one genome/assembly per file, contigs allowed).
#
# ---------------------------------------------------------------------------

set -euo pipefail

# ---------------------------- defaults --------------------------------
INPUT_DIR=""
OUTPUT_DIR=""
TOOL=""
THREADS=8
RIBAP_PROFILE="local"     # switch to "slurm" if you have the HPC profile set up
PROKKA_KINGDOM="Bacteria"
SIMILARITY=0.95            # Panaroo/Roary-style default; loosen for genus-level
usage() {
  cat <<EOF
Usage: $0 -i <genome_dir> -o <output_dir> -t <ribap|panaroo|ppanggolin|all> [options]

Required:
  -i, --input       Directory of genome FASTA files (.fasta/.fa/.fna)
  -o, --output       Output directory (will be created)
  -t, --tool         Which pipeline to run: ribap | panaroo | ppanggolin | all

Optional:
  -c, --threads       Number of CPU threads (default: $THREADS)
  -s, --similarity    Sequence similarity threshold for Panaroo/PPanGGOLiN
                       (default: $SIMILARITY). Lower this (e.g. 0.7-0.8) for
                       genus-level / cross-species comparisons.
  -p, --profile       Nextflow profile for RIBAP: local | slurm (default: local)
  -h, --help          Show this help

Examples:
  $0 -i genomes/ -o results/ -t ribap
  $0 -i genomes/ -o results/ -t all -c 16 -s 0.8
EOF
  exit 1
}

# ---------------------------- arg parsing --------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)      INPUT_DIR="$2"; shift 2 ;;
    -o|--output)     OUTPUT_DIR="$2"; shift 2 ;;
    -t|--tool)       TOOL="$2"; shift 2 ;;
    -c|--threads)    THREADS="$2"; shift 2 ;;
    -s|--similarity) SIMILARITY="$2"; shift 2 ;;
    -p|--profile)    RIBAP_PROFILE="$2"; shift 2 ;;
    -h|--help)       usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

[[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$TOOL" ]] && usage
[[ ! -d "$INPUT_DIR" ]] && { echo "ERROR: input dir '$INPUT_DIR' not found"; exit 1; }

case "$TOOL" in
  ribap|panaroo|ppanggolin|all) ;;
  *) echo "ERROR: --tool must be one of: ribap, panaroo, ppanggolin, all"; usage ;;
esac

mkdir -p "$OUTPUT_DIR"
LOGDIR="$OUTPUT_DIR/logs"
mkdir -p "$LOGDIR"

GENOMES=$(find "$INPUT_DIR" -maxdepth 1 \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) | sort)
N_GENOMES=$(echo "$GENOMES" | grep -c . || true)

if [[ "$N_GENOMES" -eq 0 ]]; then
  echo "ERROR: no .fasta/.fa/.fna files found in $INPUT_DIR"
  exit 1
fi

echo "=================================================================="
echo " Combined pangenome pipeline"
echo " Input genomes : $N_GENOMES (from $INPUT_DIR)"
echo " Output dir    : $OUTPUT_DIR"
echo " Tool(s)       : $TOOL"
echo " Threads       : $THREADS"
echo "=================================================================="

# ---------------------------------------------------------------------------
# STEP 1: Prokka annotation (needed for Panaroo / PPanGGOLiN only)
# ---------------------------------------------------------------------------
run_prokka() {
  local PROKKA_DIR="$OUTPUT_DIR/prokka"
  mkdir -p "$PROKKA_DIR"

  echo "[1/2] Running Prokka annotation on $N_GENOMES genomes..."
  while read -r genome; do
    [[ -z "$genome" ]] && continue
    base=$(basename "$genome")
    sample="${base%.*}"
    if [[ -f "$PROKKA_DIR/$sample/$sample.gff" ]]; then
      echo "  - $sample already annotated, skipping"
      continue
    fi
    echo "  - annotating $sample"
    prokka --outdir "$PROKKA_DIR/$sample" \
           --prefix "$sample" \
           --kingdom "$PROKKA_KINGDOM" \
           --cpus "$THREADS" \
           --force \
           "$genome" \
           > "$LOGDIR/prokka_${sample}.log" 2>&1
  done <<< "$GENOMES"

  echo "$PROKKA_DIR"
}

# ---------------------------------------------------------------------------
# STEP 2a: RIBAP (raw genomes in, handles Prokka + Roary + ILP internally)
# ---------------------------------------------------------------------------
run_ribap() {
  echo "[RIBAP] Launching Nextflow pipeline..."
  local RIBAP_DIR="$OUTPUT_DIR/ribap"
  mkdir -p "$RIBAP_DIR"

  # RIBAP expects a directory of genome fastas as --fasta_path (glob) or a
  # samplesheet depending on version; consult `nextflow run hoelzer-lab/ribap --help`
  # for your installed version's exact flags. This is the common invocation:
  nextflow run hoelzer-lab/ribap \
    -r master \
    -profile "$RIBAP_PROFILE" \
    --cores "$THREADS" \
    --fasta "$INPUT_DIR/*.{fasta,fa,fna}" \
    --outdir "$RIBAP_DIR" \
    -resume \
    2>&1 | tee "$LOGDIR/ribap.log"

  echo "[RIBAP] Done. Core gene groups (RIBAP groups) written to: $RIBAP_DIR"
}

# ---------------------------------------------------------------------------
# STEP 2b: Panaroo (needs Prokka GFF3s)
# ---------------------------------------------------------------------------
run_panaroo() {
  local PROKKA_DIR="$1"
  echo "[Panaroo] Collecting GFF3 files..."
  local GFFS
  GFFS=$(find "$PROKKA_DIR" -name "*.gff" | sort)
  local PANAROO_DIR="$OUTPUT_DIR/panaroo"
  mkdir -p "$PANAROO_DIR"

  echo "[Panaroo] Running (mode=strict, threshold=$SIMILARITY)..."
  # mode: strict (single species) | moderate | sensitive (more diverse input)
  panaroo -i $GFFS \
          -o "$PANAROO_DIR" \
          --clean-mode strict \
          -t "$THREADS" \
          --threshold "$SIMILARITY" \
          2>&1 | tee "$LOGDIR/panaroo.log"

  echo "[Panaroo] Done. Core/accessory gene tables in: $PANAROO_DIR"
}

# ---------------------------------------------------------------------------
# STEP 2c: PPanGGOLiN (needs Prokka GFF3s listed in an organisms.tsv)
# ---------------------------------------------------------------------------
run_ppanggolin() {
  local PROKKA_DIR="$1"
  echo "[PPanGGOLiN] Building organisms.tsv..."
  local PPANG_DIR="$OUTPUT_DIR/ppanggolin"
  mkdir -p "$PPANG_DIR"
  local ORG_TSV="$PPANG_DIR/organisms.tsv"
  > "$ORG_TSV"

  while read -r genome; do
    [[ -z "$genome" ]] && continue
    base=$(basename "$genome")
    sample="${base%.*}"
    gff="$PROKKA_DIR/$sample/$sample.gff"
    if [[ -f "$gff" ]]; then
      echo -e "${sample}\t${gff}" >> "$ORG_TSV"
    fi
  done <<< "$GENOMES"

  echo "[PPanGGOLiN] Running 'ppanggolin all' (identity=$SIMILARITY)..."
  ppanggolin all \
    --anno "$ORG_TSV" \
    --output "$PPANG_DIR/run" \
    --identity "$SIMILARITY" \
    --cpu "$THREADS" \
    -f \
    2>&1 | tee "$LOGDIR/ppanggolin.log"

  echo "[PPanGGOLiN] Done. Partitions (persistent/shell/cloud) in: $PPANG_DIR/run"
}

# ---------------------------------------------------------------------------
# DISPATCH
# ---------------------------------------------------------------------------
PROKKA_DIR=""
if [[ "$TOOL" == "panaroo" || "$TOOL" == "ppanggolin" || "$TOOL" == "all" ]]; then
  PROKKA_DIR=$(run_prokka)
fi

case "$TOOL" in
  ribap)
    run_ribap
    ;;
  panaroo)
    run_panaroo "$PROKKA_DIR"
    ;;
  ppanggolin)
    run_ppanggolin "$PROKKA_DIR"
    ;;
  all)
    run_ribap
    run_panaroo "$PROKKA_DIR"
    run_ppanggolin "$PROKKA_DIR"
    echo ""
    echo "=================================================================="
    echo " All three tools finished. Compare core genome sizes with:"
    echo "   RIBAP:      wc -l $OUTPUT_DIR/ribap/**/ribap_groups.csv (approx. path — check RIBAP outdir)"
    echo "   Panaroo:    grep -c 'Yes' $OUTPUT_DIR/panaroo/gene_presence_absence.csv (or inspect summary_statistics.txt)"
    echo "   PPanGGOLiN: cat $OUTPUT_DIR/ppanggolin/run/genomes_statistics.tsv"
    echo "=================================================================="
    ;;
esac

echo ""
echo "Pipeline finished. Logs in: $LOGDIR"
