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

  run command -v csvtk
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

@test "serogroup predictions gene-by-gene" {
  genesIdx=$(head -n 1 $serogroup_tab)
}

@test "serogroup predictions exactly the same" {
  header=$(head -n1 $serogroup_tab | tr '\t' ',')
  # use csvtk to reorder the columns
  csvtk -t cut -f "$header" $serogroup_tab > $BATS_TEST_TMPDIR/obs.tab
  csvtk -t cut -f "$header" ${testDir}/serogroup_predictions.tab > $BATS_TEST_TMPDIR/exp.tab

  # compare the Genes_Present column's values.
  genesIdx=$(head -n 1 $serogroup_tab | tr '\t' '\n' | grep -n '^Genes_Present$' | cut -d: -f1)
  obs_genes=$(cut -f$genesIdx $BATS_TEST_TMPDIR/obs.tab | tail -n +2 | tr ',' '\n' | sort)
  exp_genes=$(cut -f$genesIdx $BATS_TEST_TMPDIR/exp.tab | tail -n +2 | tr ',' '\n' | sort)
  #echo -e "obs_genes:\n$obs_genes" >&3
  #echo -e "exp_genes:\n$exp_genes" >&3
  # go gene by gene to see which are different and err on the first one that's different
  # print any diffs to >&3
  run bash -c '
    paste <(echo "$obs_genes") <(echo "$exp_genes") | awk -F"\t" "
      NR==1 { next }
      {
        if (\$1 != \$2) {
          # print to >&3
          printf \"# Gene %s differs: %s vs %s\n\", NR, \$1, \$2 > "/dev/stderr"
          diff_found=1
        }
      }
      END { exit diff_found }
    "
  '
  [ "$status" -eq 0 ]

  run bash -c '
    paste exp.tab obs.tab | awk -F"\t" "
      NR==1 {
        n = NF/2
        for (i=1; i<=n; i++) colname[i]=\$i
        next
      }
      {
        for (i=1; i<=n; i++) {
          if (\$i != \$(i+n)) {
            printf \"Row %d, Column %s differs: %s vs %s\n\", NR, colname[i], \$i, \$(i+n)
            diff_found=1
          }
        }
      }
      END { exit diff_found }
    "
  '
  [ "$status" -eq 0 ]
}