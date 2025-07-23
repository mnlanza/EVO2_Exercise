#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Running EVO2 Analysis ===${NC}\n"

# Function to check if a command was successful
check_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ Analysis completed successfully${NC}\n"
    else
        echo -e "\n❌ Analysis failed"
        exit 1
    fi
}

# Create output directories if they don't exist
mkdir -p plots/ex1_plots plots/ex1_shifting_plots plots/rc_Ex1

# Standard variant analysis
echo "Running standard variant analysis..."
Rscript scripts/Ex1_main.R \
    --fasta raw_data/gene_variants.fasta \
    --variants 83_S1,83_L,83_S2,83_Stop \
    --input_logits_dir raw_data \
    --output_dir plots/ex1_plots \
    --only_single

check_status

# Shifting positions analysis
echo -e "\nRunning shifting positions analysis..."
Rscript scripts/Ex1_main.R \
    --fasta raw_data/shifted_variants.fasta \
    --variants 83_S1,83_S2,83_S1_shifted1r,83_S1_shifted2r \
    --input_logits_dir raw_data \
    --output_dir plots/ex1_shifting_plots \
    --reverse_comp \
    --only_single

check_status

# Reverse complement analysis
############ TODO: add reverse complement analysis (get logits for rc strands)
echo -e "\nRunning reverse complement analysis..."
Rscript scripts/Ex1_main.R \
    --fasta raw_data/reverse_complements_83.fasta \
    --variants 83_AGC_rc,83_ATC_rc,83_AGT_rc,83_TGA_rc \
    --input_logits_dir raw_data \
    --output_dir plots/rc_Ex1 \
    --reverse_comp

check_status

echo -e "${GREEN}All analyses completed successfully!${NC}"
echo "Results can be found in:"
echo "  - plots/ex1_plots"
echo "  - plots/ex1_shifting_plots"
echo "  - plots/rc_Ex1" 