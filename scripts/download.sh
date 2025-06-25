#!/bin/bash
# Fail on any error
set -euo pipefail
LOG=$4 # Assuming $4 is like your !{LOG}
to_download=$6   # this is a directory location
download_start=$SECONDS
if [ "$1" == "T" ]; then        # Assuming $1 is like your !{REFS}
  echo "Debug: reference_db=$2" # $2 is like your !{params.reference_db}
  out="$2/collect_refs"
  mkdir -p "$out"
  cp $5 $out
  cd "$out"
  declare -A downloaded
  downloadedList=$3 # $3 is like your !{params.downloadedList}
  script_dir="$(cd "$(dirname "$0")" && pwd)"
  file_access_locker="$script_dir/utils/file_access_locker.sh"
  download_lock_file="$out/download_lock_file"

  bash "$file_access_locker" "$download_lock_file" "lock" "$LOG"

  if [ ! -f "$downloadedList" ]; then
    touch "$downloadedList"
  fi
  while read -r line; do
    downloaded[$line]=1
  done <"$downloadedList"

  echo "DEBUG: Checking current directory"
  pwd

  # Read accession list and download only the ones not downloaded yet
  if [ -e "${to_download}/accessions_to_download.txt" ]; then
    rm ${to_download}/accessions_to_download.txt
  fi
  touch ${to_download}/accessions_to_download.txt
  while read -r acc || [[ -n "$acc" ]]; do
    acc=$(echo $acc | tr -d '[:space:]')                    # Remove any extra spaces
    [ -z "$acc" ] && echo "Error: acc is empty" && continue # Print an error message and skip to the next iteration if acc is empty
    echo acc $acc
    # debug: print out the key and the corresponding value from the array
    if [[ -z ${downloaded[$acc]+x} ]]; then
      echo "Accession $acc is not downloaded. Downloading now."
      echo $acc >>${to_download}/accessions_to_download.txt
      echo $acc >>"$downloadedList"
    fi
  done <"$5" # Assuming $5 is like your !{params.ref_candidates}
  # The actual downloading part
  if [ -s ${to_download}/accessions_to_download.txt ]; then
    attempt=1
    while [ $attempt -le 3 ]; do
        datasets download genome accession --inputfile ${to_download}/accessions_to_download.txt --filename refs.zip 2>&1 | tee -a "$LOG"
        if [ ! -s refs.zip ]; then
	    echo "Download failed or refs.zip is empty. Retrying..."
  
            continue
        fi
	unzip -n refs.zip
        unzip_exit_status=$?
        if [ $unzip_exit_status -eq 0 ]; then
            rm -f refs.zip
            rm -f README.md
            cat ${to_download}/accessions_to_download.txt >>"$downloadedList"
            break
        else
            echo "Attempt $attempt to unzip failed. Retrying..."
            rm -f refs.zip
            ((attempt++))
        fi
    done

    if [ $attempt -gt 3 ]; then
        echo "process collect_refs failed after 3 attempts"
    fi
  else
    echo "No accessions to download"
  fi
  cd -
  REF_CANDS="$out/ncbi_dataset/data/"
  # Change permissions of downloaded files
  chmod -R 777 "$REF_CANDS" || true
  bash "$file_access_locker" "$download_lock_file" "unlock" "$LOG"
  echo "Unlocked"
  printf %s "NOTE: All references may not be available for download from NCBI.\n" >>"$LOG"
else
  printf %s "-----Refs found! Skipping collect_refs-----\n$(date)\n\n" >>"$LOG"
fi

download_end=$SECONDS
download_time=$((download_end - download_start))
echo "Download time: $download_time" >>"$LOG"
