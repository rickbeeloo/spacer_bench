prefixes=("6b1148fc" "01da6dbd" "f25c897e" "6ef5b781" "4b7b4f05" "adb3760b" "6951d55c" "936dac80" "a238f23b" "493dc49e" "6a19eb7c" "59e8522e" "6b1148fc" "578a5025" "37d6023d" "9be1d7f4" "93cc9b21")

# function to return true if  a file is  relevant  (.txt, .tsv, .sam, .tab)
is_relevant_file() {
  local file="$1"
  if [[ "$file" == *.txt || "$file" == *.tsv || "$file" == *.sam || "$file" == *.tab || "$file" == *.fa ]]; then
    return 0
  else
    return 1
  fi
}

process_file() {
  local file="$1"
  local found_prefixes=()
  
  for prefix in "${prefixes[@]}"; do
    if grep -q -F -m 1 "${prefix}_" "$file"; then
      found_prefixes+=("$prefix")
    fi
  done
  
  if [ ${#found_prefixes[@]} -gt 0 ]; then
    echo "Found prefixes in $file: ${found_prefixes[*]}"
    # Escape each prefix for sed and join with |
    local regex_pattern=""
    for prefix in "${found_prefixes[@]}"; do
      if [ -n "$regex_pattern" ]; then
        regex_pattern="${regex_pattern}|"
      fi
      # Escape special characters for sed
      prefix_escaped=$(echo "$prefix" | sed 's/[]\/$*.^[]/\\&/g')
      regex_pattern="${regex_pattern}${prefix_escaped}"
    done
    echo "Using pattern: ${regex_pattern}"
    echo "Before processing $file:"
    head -n 1 "$file"
    sed -i "s/${regex_pattern}_//g" "$file"
    echo "After processing $file:"
    head -n 1 "$file"
  fi
}

for run_dir in ./results/simulated/run_t_*/; do
  echo "run_dir $run_dir"
  for file in "$run_dir"/*; do
    if [ -f "$file" ]; then
      echo "$file"
      if is_relevant_file "$file"; then
        process_file "$file"
      fi
    fi
    
    if [ -d "$file" ]; then
      echo "$file"
      for subfile in "$file"/*; do
        if [ -f "$subfile" ]; then
          if is_relevant_file "$subfile"; then
            process_file "$subfile"
          fi
        fi
      done
    fi
  done
done 