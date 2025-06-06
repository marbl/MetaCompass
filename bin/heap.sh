#!/bin/bash

function find_best_ref {
  local all_reads_fasta=$1
  local subset_reads_fasta=$2
  local fasta_dir=$3
  local selected_refs=$4

  echo "all_reads_fasta ${all_reads_fasta}"
  echo "subset_reads_fasta ${subset_reads_fasta}"
  echo "fasta_dir ${fasta_dir}"
  echo "selected_refs ${selected_refs}"

  [ ! -f $selected_refs ] && touch $selected_refs

  local subset_reads_dir=$(realpath $(dirname "$subset_reads_fasta"))
  local tmpdir=$(mktemp -d -p "$subset_reads_dir" "kmc-XXXXXX")
  mkdir ${tmpdir}/subset_reads
  echo "Temp directory is ${tmpdir}"

  dir_name_1=$(realpath $(dirname "$all_reads_fasta"))

  local all_reads="${dir_name_1}/kmc_last_used.kmc"
  if [ ! -f ${all_reads}.kmc_pre ]; then
    all_reads=$1
  fi

  if [ -s $subset_reads_fasta ]; then
    kmc -k28 -t8  -hp -ci1  -fq $subset_reads_fasta "${tmpdir}/subset_reads_kmc_tmp" ${tmpdir}/subset_reads
    kmc_tools -t8 simple  ${all_reads} ${tmpdir}/subset_reads_kmc_tmp counters_subtract ${tmpdir}/difference_kmc_tmp
    kmc_file_to_use=${tmpdir}/difference_kmc_tmp
  else
    echo "Subset file is empty or doesn't exist. Skipping subtraction."
    kmc_file_to_use=${all_reads}
  fi

  echo "Selected kmc file ${kmc_file_to_use}"
  best_ref_count=0
  best_ref=""

  for ref_kmc_file in $(find $fasta_dir -name "cluster*.kmc_pre*"); do
    if [[ ! -f $ref_kmc_file ]]; then
      echo "Warning: $ref_kmc_file does not exist or cannot be read."
      continue
    fi

    ref_basename=$(basename "$ref_kmc_file" .kmc_pre)
    dir_name=$(realpath $(dirname "$ref_kmc_file"))
    echo ${ref_basename}
    if ! grep -qF "$ref_basename" $selected_refs; then
       echo kmc_tools -t8 -hp simple $kmc_file_to_use -ci1 ${dir_name}/${ref_basename} -ci1 intersect ${tmpdir}/${ref_basename}.intersect
       kmc_tools -t8 -hp simple $kmc_file_to_use -ci1 ${dir_name}/${ref_basename}  -ci1 intersect ${tmpdir}/${ref_basename}.intersect
       total_kmers=$(kmc_tools info "${tmpdir}/${ref_basename}.intersect" | grep "total k-mers" | awk -F':' '{print $2}' | tr -d '[:space:]')

      if [[ ! "$total_kmers" =~ ^[0-9]+$ ]]; then
        echo "Warning: Could not obtain a valid k-mer count for $ref_kmc_file."
        continue
      fi

      if [ "$total_kmers" -gt "$best_ref_count" ]; then
        best_ref_count=$total_kmers
        best_ref=$ref_basename
      fi
    fi
  done



  if [ -n "$best_ref" ]; then
    echo "$best_ref" >> $selected_refs
  fi
  echo $kmc_file_to_use
  if [ -f ${kmc_file_to_use}.kmc_pre ]; then
    echo $kmc_file_to_use
    mv ${kmc_file_to_use}.kmc_pre "${dir_name_1}/kmc_last_used.kmc.kmc_pre"
    mv  ${kmc_file_to_use}.kmc_suf "${dir_name_1}/kmc_last_used.kmc.kmc_suf"

  fi
  rm -r ${tmpdir}


  echo "Selected Reference: $best_ref, Intersections: $best_ref_count"
}

find_best_ref "$@"
