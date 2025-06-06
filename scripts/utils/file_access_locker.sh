#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <file> <lock|unlock>"
    exit 1
fi

# Define the file to lock/unlock
file="$1"
lock_option="$2"
LOG="$3"

# Create a lock file name based on the provided file
lock_file="${file}.lock"

# Function to lock the file
lock_file() {
    while [ -e "$lock_file" ]; do
        echo "Lock file $lock_file already exists. File is locked. Waiting for unlock..." >>"$LOG"
        sleep 1
    done

    touch "$lock_file"
    echo "File locked: $file" >>"$LOG"
}

# Function to unlock the file
unlock_file() {
    if [ ! -e "$lock_file" ]; then
        echo "Lock file $lock_file does not exist. File is not locked." >>"$LOG"
        exit 1
    fi

    rm -f "$lock_file"
    echo "File unlocked: $file" >>"$LOG"
}

# Check the lock option and perform the corresponding action
case "$lock_option" in
    "lock")
        lock_file
        ;;
    "unlock")
        unlock_file
        ;;
    *)
        echo "Invalid option. Use 'lock' or 'unlock'." >>"$LOG"
        exit 1
        ;;
esac

exit 0
