# Use the current working directory as the search directory
search_dir="$PWD"

# Output file to save the list of .SAO filenames
output_file="sao_files_list.txt"

# Find all .SAO files and save only the filenames to the output file
find "$search_dir" -type f -iname "*.sao" -exec basename {} \; > "$output_file"

# Optional: Display the content of the output file
cat "$output_file"