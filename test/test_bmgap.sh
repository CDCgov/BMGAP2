#!/usr/bin/env bats
# Test BMGAP results for a specific SRR which is exported in the environment
# assume that conda environment or other environment is already activated

# Check expected environment variables
set -u

# this directory is the directory that this script is in
export thisDir=$(dirname "${BATS_TEST_FILENAME}")
# The data for this SRR test is in $SRR.exp in the test folder
export testDir="$thisDir/$SRR.exp"

# test that jq is installed
@test "jq is installed and env variables are found" {
  echo "# SRR: $SRR" >&3
  echo "# indir: $indir" >&3
  echo "# species_analysis_json: $species_analysis_json" >&3
  echo "# molecular_data_json: $molecular_data_json" >&3
  echo "# scheme_counts_json: $scheme_counts_json" >&3
  echo "# amr_data_json: $amr_data_json" >&3

  run command -v jq
  [ "$status" -eq 0 ]
}

# test that species is Neisseria meningitidis from the json
@test "species_analysis" {
  # get species analysis json. Remove mash_hash key since it can vary slightly.
  # I'm sure it can be done more correctly with jq instead of sed but this is okay for now.
  obs=$(jq --sort-keys . $species_analysis_json | sed -e '/mash_hash/d' -e '/score/d')
  exp=$(jq --sort-keys . ${testDir}/species_analysis.json | sed -e '/mash_hash/d' -e '/score/d')
  echo -e "obs: $obs" | sed 's/^/# /' >&3
  echo -e "exp: $exp" | sed 's/^/# /' >&3
  [[ "$obs" == "$exp" ]]
}

@test "molecular_data" {
  jq_del='del(.Filename, .Analysis_User, .Analysis_Time, .Lab_ID, .Analysis_Version)'
  sed_del="/Not found\|null\|Not applicable\|Error\|None/d"

  mlst_keys=""
  for locus in abcZ adk aroE fumC gdh pdhC pgm ST cc; do
    mlst_keys="$mlst_keys .Nm_MLST_$locus"
  done
  for locus in adk atpG frdB fucK mdh pgi recA ST; do
    mlst_keys="$mlst_keys .Hi_MLST_$locus"
  done

  for key in $mlst_keys; do
    obs=$(jq --sort-keys $key $molecular_data_json | sed -e 's/,$//' -e "$sed_del")
    exp=$(jq --sort-keys $key ${testDir}/molecular_data.json | sed -e 's/,$//' -e "$sed_del")
    echo -e "# checking key $key" >&3
    echo -e "obs: $obs" | sed 's/^/# /' >&3
    echo -e "exp: $exp" | sed 's/^/# /' >&3
    [[ "$obs" == "$exp" ]]
  done
  
}

@test "AMR" {
  obs=$(jq --sort-keys .antimicrobics $amr_data_json)
  exp=$(jq --sort-keys .antimicrobics ${testDir}/amr_data.json)
  [[ "$obs" == "$exp" ]]
}