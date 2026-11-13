#!/bin/bash

# Check if the filename is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <filename.tex>"
  exit 1
fi

# Remove the file extension
filename="${1%.*}"

# Run xelatex, bibtex, and xelatex two more times
echo "Running xelatex on $filename.tex..."
xelatex "$filename.tex"

echo "Running bibtex on $filename.aux..."
bibtex "$filename"

echo "Running xelatex again on $filename.tex..."
xelatex "$filename.tex"

echo "Running xelatex one last time on $filename.tex..."
xelatex "$filename.tex"

echo "Compilation finished. Check for any errors in the output above."

# Optional: Open the PDF file after compiling
# open "$filename.pdf"  # Uncomment for macOS
# xdg-open "$filename.pdf"  # Uncomment for Linux
