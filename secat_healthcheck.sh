#!/bin/bash
# =============================================================================
# SeCAT health check  (multi-stage, system-agnostic)
#
# USAGE:  ./secat_healthcheck.sh [pre|post1|post2|all]
#   pre    Readiness check to run BEFORE launching: manifest, per-study input
#          files, params, and roster are all present and consistent.
#   post1  After phase 1 (analysis -> verdicts -> report).
#   post2  After --step standardize (merged dataset, verification, comparison).
#   all    pre + post1 + post2.
#   (default) auto: always run 'pre'; add post1/post2 if their outputs exist.
#
# System-agnostic: all paths are read from params.yaml (override the params file
# with SECAT_PARAMS=/path/to/params.yaml). No hard-coded locations. Run it from
# the repo root, or from anywhere — it locates itself.
# =============================================================================
set -uo pipefail
MODE="${1:-auto}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
PARAMS="${SECAT_PARAMS:-params.yaml}"

C_RED='\033[0;31m'; C_GREEN='\033[0;32m'; C_YELLOW='\033[0;33m'; C_BLUE='\033[0;34m'; C_NC='\033[0m'
PASS=0; FAILN=0; WARNN=0
hdr()  { echo -e "\n${C_BLUE}=== $1 ===${C_NC}"; }
ok()   { echo -e "  ${C_GREEN}OK:${C_NC} $1"; PASS=$((PASS+1)); }
bad()  { echo -e "  ${C_RED}FAIL:${C_NC} $1"; FAILN=$((FAILN+1)); }
warn() { echo -e "  ${C_YELLOW}WARN:${C_NC} $1"; WARNN=$((WARNN+1)); }

# tiny reader for a top-level "key: value" scalar in params.yaml
yget() { grep -E "^[[:space:]]*$1[[:space:]]*:" "$PARAMS" 2>/dev/null | head -1 \
         | sed -E 's/^[^:]*:[[:space:]]*//; s/[[:space:]]*#.*$//; s/[[:space:]]+$//' \
         | sed "s/[\"']//g"; }

[ -f "$PARAMS" ] || { echo -e "${C_RED}Cannot find params file: $PARAMS${C_NC}"; exit 2; }
MANIFEST="$(yget manifest)";   MANIFEST="${MANIFEST:-secat_manifest.tsv}"
OUTDIR="$(yget outdir)";       OUTDIR="${OUTDIR:-output}"
FINAL="$(yget final_outputs)"; FINAL="${FINAL:-$OUTDIR/final_outputs}"
SEL_MODE="$(yget selection_mode)"; SEL_MODE="${SEL_MODE:-auto}"
SEL_FILE="$(yget selection_file)"
RPT=1_report; DS=2_dataset; VD=3_verdicts; VF=4_verification; CMP=5_comparison

check_pre() {
  hdr "PRE-RUN READINESS"
  if [ -f "$MANIFEST" ]; then
    ok "manifest found ($MANIFEST)"
    local h; h="$(head -1 "$MANIFEST")"
    for col in study_name asv_fasta_path asv_counts_path taxonomy_path metadata_path; do
      printf '%s' "$h" | tr '\t' '\n' | grep -qx "$col" && ok "column '$col' present" || bad "manifest missing column '$col'"
    done
    local miss; miss="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}{split("asv_fasta_path asv_counts_path taxonomy_path metadata_path",P," ");for(k in P){col=P[k];if(col in c){p=$(c[col]);if(p!="" && system("[ -f \"" p "\" ]")!=0)print "    missing "col" for "$(c["study_name"])": "p}}}' "$MANIFEST")"
    if [ -z "$miss" ]; then ok "all per-study input files exist"; else bad "some manifest input files are missing:"; echo -e "${C_RED}$miss${C_NC}"; fi
  else
    bad "manifest NOT found ($MANIFEST)"
  fi
  if [ "$SEL_MODE" = "roster" ] || [ "$SEL_MODE" = "file" ]; then
    if [ -n "${SEL_FILE:-}" ] && [ -f "$SEL_FILE" ]; then
      ok "selection_mode=$SEL_MODE; selection_file found ($SEL_FILE)"
      if [ -f "$MANIFEST" ]; then
        local studies; studies="$(grep -vE '^[[:space:]]*#|^[[:space:]]*$' "$SEL_FILE" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
        local known; known="$(cut -f1 "$MANIFEST" | tail -n +2)"
        local notfound=""
        while IFS= read -r s; do [ -z "$s" ] && continue; grep -qxF "$s" <<<"$known" || notfound+="    not in manifest: $s"$'\n'; done <<<"$studies"
        if [ -z "$notfound" ]; then ok "all rostered studies exist in manifest"; else warn "rostered studies not in manifest (will be skipped):"; echo -e "${C_YELLOW}$notfound${C_NC}"; fi
      fi
    else
      bad "selection_mode=$SEL_MODE but selection_file missing (${SEL_FILE:-unset})"
    fi
  else
    ok "selection_mode=$SEL_MODE (no roster file required)"
  fi
}

check_post1() {
  hdr "POST PHASE 1 (verdicts + report)"
  [ -f "$FINAL/$RPT/SeCAT_Master_Summary_Report.pdf" ] && ok "master report ($RPT/SeCAT_Master_Summary_Report.pdf)" || bad "master report missing"
  local n; n=$(find "$FINAL/$RPT/per_study" -name '*.pdf' 2>/dev/null | wc -l); [ "$n" -gt 0 ] && ok "$n per-study report(s)" || bad "no per-study reports"
  [ -f "$FINAL/$VD/master_verdict_table.csv" ] && ok "master_verdict_table.csv" || bad "master_verdict_table.csv missing"
}

check_post2() {
  hdr "POST STANDARDIZE (merged dataset + verification)"
  for f in combined_feature_table.tsv combined_taxonomy.tsv combined_sequences.fasta; do
    [ -f "$FINAL/$DS/$f" ] && ok "$DS/$f" || bad "$DS/$f missing"
  done
  local n; n=$(find "$FINAL/$VF" -name 'permanova_results.csv' 2>/dev/null | wc -l); [ "$n" -gt 0 ] && ok "verification stats present ($n level(s))" || warn "no verification stats (run_validation off?)"
  local c; c=$(find "$FINAL/$CMP" -name '*.tsv' 2>/dev/null | wc -l); [ "$c" -gt 0 ] && ok "$CMP has $c table(s)" || bad "$CMP empty"
}

echo -e "${C_BLUE}SeCAT health check${C_NC}  ($(date))  params=$PARAMS  final=$FINAL"
case "$MODE" in
  pre)   check_pre ;;
  post1) check_post1 ;;
  post2) check_post2 ;;
  all)   check_pre; check_post1; check_post2 ;;
  auto|*) check_pre
          [ -d "$FINAL/$RPT" ] && check_post1
          [ -d "$FINAL/$DS" ]  && check_post2 ;;
esac
echo -e "\n${C_BLUE}=== SUMMARY ===${C_NC}  ${C_GREEN}$PASS ok${C_NC}, ${C_YELLOW}$WARNN warn${C_NC}, ${C_RED}$FAILN fail${C_NC}"
[ "$FAILN" -eq 0 ]
