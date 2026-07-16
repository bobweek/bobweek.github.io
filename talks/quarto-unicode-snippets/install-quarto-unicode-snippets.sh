#!/usr/bin/env bash
set -euo pipefail

SRC_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SNIP="quarto-unicode-math.code-snippets"

targets=()

if command -v codium >/dev/null 2>&1 || [ -d "$HOME/.config/VSCodium" ]; then
  targets+=("$HOME/.config/VSCodium/User/snippets")
fi

if command -v code >/dev/null 2>&1 || [ -d "$HOME/.config/Code" ]; then
  targets+=("$HOME/.config/Code/User/snippets")
fi

if [ "${#targets[@]}" -eq 0 ]; then
  targets+=("$HOME/.config/VSCodium/User/snippets")
fi

for dir in "${targets[@]}"; do
  mkdir -p "$dir"
  cp "$SRC_DIR/$SNIP" "$dir/$SNIP"
  echo "installed: $dir/$SNIP"
done

echo
echo "Reload VS Code/VSCodium after installation."
echo "Use by typing e.g. alpha<Tab>, \\alpha<Tab>, rho<Tab>, \\rho<Tab>."
