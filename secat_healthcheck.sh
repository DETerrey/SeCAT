#!/bin/bash
# ===================================================================
# MESAP Comprehensive Pipeline Health Check (V2.0)
#
# USAGE: Run from the MESAP project root after a pipeline run.
#        ./secat_healthcheck.sh
# ===================================================================

# --- Color Definitions ---
C_RED='\033[0;31m'
C_GREEN='\033[0;32m'
C_YELLOW='\033[0;33m'
C_BLUE='\033[0;34m'
C_NC='\033[0m' # No Color

# --- Helper Functions ---
print_header() { echo -e "\n${C_BLUE}=== $1 ===${C_NC}"; }
print_success() { echo -e "  ${C_GREEN}✓ SUCCESS:${C_NC} $1"; }
print_failure() { echo -e "  ${C_RED}✗ FAILURE:${C_NC} $1"; }
print_warning() { echo -e "  ${C_YELLOW}⚠ WARNING:${C_NC} $1"; }

check_file_exists() {
    if [ -f "$2" ]; then print_success "$1 found ($2)"; return 0;
    else print_failure "$1 NOT FOUND ($2)"; return 1; fi
}

check_dir_populated() {
    if [ -d "$2" ] && [ -n "$(find "$2" -name "$3" -print -quit)" ]; then
        local count=$(find "$2" -name "$3" | wc -l)
        print_success "$1 contains $count files ($2)"
        return 0
    else
        print_failure "$1 is empty or missing ($2)"
        return 1
    fi
}

# ===============================================================
# --- SCRIPT START ---
# ===============================================================
print_header "MESAP Pipeline Health Check"
echo "Running diagnostics on $(date)"

# --- Step 1: Input Validation ---
print_header "Step 1: Input Manifest"
check_file_exists "Manifest file" "secat_manifest.tsv"
if [ -f "secat_manifest.tsv" ]; then
    echo "  -> Validating paths in manifest..."
    # Quick check using awk to verify files exist
    awk -F'\t' 'NR>1 {if(system("[ -f " $4 " ]") != 0) print "Missing ASV file: " $4}' secat_manifest.tsv
fi

# --- Step 2: Intermediate Files ---
print_header "Step 2: Intermediate Artifacts"
# Phase 1/2 (Mode Dependent)
if [ -f "output/intermediate/study_alignment_coords.csv" ]; then
    print_success "Study Mode coordinates found"
    check_dir_populated "Aligned FASTAs" "output/intermediate/aligned_fastas" "*_aligned.fasta"
elif [ -f "output/intermediate/primer_coords_phase1_output.csv" ]; then
    print_success "Primer Mode coordinates found"
    check_dir_populated "Primer Databases" "output/primer_databases" "*.fasta"
else
    print_warning "No coordinate files found (neither Study nor Primer mode artifacts)"
fi

check_file_exists "Simulation Task List" "output/intermediate/simulation_tasks.csv"

# --- Step 3: Analysis Results ---
print_header "Step 3: Analysis Results"
check_dir_populated "Simulation Results (RDS)" "output/simulation_results" "*results.rds"
check_dir_populated "Real Data Results (RDS)" "output/real_data_results" "*_results.rds"

# --- Step 4: Aggregation ---
print_header "Step 4: Aggregated Data"
check_file_exists "Baseline Statistics" "output/aggregated_data/simulation_baseline_statistics.csv"
check_file_exists "Retention Curves" "output/aggregated_data/simulation_retention_curves.csv"
check_file_exists "Master Verdict Table" "output/aggregated_data/master_verdict_table.csv"

# --- Step 5: Final Report ---
print_header "Step 5: Final Report & Plots"
check_file_exists "Master Summary PDF" "output/master_summary_report.pdf"
check_dir_populated "Individual Study Reports" "output/final_plots" "*.pdf"

# --- Step 6: Log Analysis ---
print_header "Step 6: Log Analysis"
if [ -d "logs" ]; then
    # Check for FATAL errors
    errors=$(grep -l -E "FATAL|Error:|Execution halted|Killed" logs/*.err 2>/dev/null)
    if [ -n "$errors" ]; then
        print_failure "Critical errors found in logs:"
        for log in $errors; do
            echo -e "    - ${C_RED}$log${C_NC}"
            # Extract the actual error message line
            grep -E "FATAL|Error:|Execution halted|Killed" "$log" | tail -n 1 | sed 's/^/      /'
        done
    else
        print_success "No critical errors found in logs."
    fi

    # Check for Warnings (Non-fatal)
    warns=$(grep -l "WARNING" logs/*.err 2>/dev/null)
    if [ -n "$warns" ]; then
        echo -e "  ${C_YELLOW}Note: Warnings found in $(echo "$warns" | wc -l) log files (check for 'WARN' or 'WARNING').${C_NC}"
    fi
else
    print_warning "Logs directory not found."
fi

echo -e "\n${C_BLUE}--- Health Check Complete ---${C_NC}"
