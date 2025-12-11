#!/bin/bash
# Create README.md file with <r^3> values for each element.

echo "| Element | <r^3> [Bohr^3] |"
echo "|---------|----------------|"

for dir in c h n o s; do
    if [ -d "$dir" ]; then
        cd "$dir"
        r3_value=$(python read_horton.py 2>/dev/null)
        cd ..
        printf "| %-7s | %-14s |\n" "${dir^^}" "$r3_value"
    fi
done
