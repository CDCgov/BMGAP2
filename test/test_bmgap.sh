#!/usr/bin/env bats
# Test BMGAP results for a specific SRR which is exported in the environment
# assume that conda environment or other environment is already activated

# Check expected environment variables
set -u
echo "SRR: $SRR"
echo "indir: $indir"
echo "species_analysis_json: $species_analysis_json"
echo "molecular_data_json: $molecular_data_json"
echo "scheme_counts_json: $scheme_counts_json"
echo "amr_data_json: $amr_data_json"

# this directory is the directory that this script is in
export thisDir=$(dirname "${BATS_TEST_FILENAME}")
# The data for this SRR test is in $SRR.exp in the test folder
export testDir="$thisDir/$SRR.exp"

# test that jq is installed
@test "jq is installed" {
  run command -v jq
  [ "$status" -eq 0 ]
}

# test that species is Neisseria meningitidis from the json
@test "species_analysis" {
  # get species analysis json. Remove mash_hash key since it can vary slightly.
  obs=$(jq --sort-keys . $species_analysis_json | sed 's/mash_hash//g')
  exp=$(jq --sort-keys . ${testDir}/species_analysis.json | sed 's/mash_hash//g')
  echo -e "obs: $obs" | sed 's/^/# /' >&3
  echo -e "exp: $exp" | sed 's/^/# /' >&3
  [[ "$obs" == "$exp" ]]
}
