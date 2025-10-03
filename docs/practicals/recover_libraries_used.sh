#!/bin/bash
# find_libraries.sh
# This script prints all lines starting with "library(" in .qmd files,
# and saves a unique sorted list of libraries to a file called "all_libraries".

# Exit if no .qmd files exist
shopt -s nullglob
files=( *.qmd )
if [ ${#files[@]} -eq 0 ]; then
  echo "No .qmd files found."
  exit 1
fi

# Temporary file to collect library names
tmpfile=$(mktemp)

echo "=== Searching for libraries in .qmd files ==="
echo

# Loop through all .qmd files
for file in "${files[@]}"; do
  echo "=== $file ==="
  # Print lines that start with library(
  grep '^library(' "$file" || echo "(none)"
  
  # Extract just the library names (handle quotes and spaces)
  grep '^library(' "$file" \
    | sed -E 's/^library\(([^)]+)\).*/\1/' \
    | tr -d '"' \
    | tr -d "'" \
    | tr -d ' ' >> "$tmpfile"
  echo
done

# Sort, remove duplicates, and overwrite output file
sort -u "$tmpfile" > all_libraries

echo "✅ Unique libraries saved to: all_libraries"
echo
echo "=== Unique Libraries Used ==="
cat all_libraries

# Clean up
rm -f "$tmpfile"
