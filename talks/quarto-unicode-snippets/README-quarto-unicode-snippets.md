# Quarto / Markdown Unicode math snippets

This installs a separate VS Code / VSCodium snippet file:

```text
quarto-unicode-math.code-snippets
```

It does not overwrite your existing Quarto layout snippets.

## Install

From this folder:

```bash
./install-quarto-unicode-snippets.sh
```

Then reload VS Code / VSCodium.

## Usage

Type a prefix and accept the snippet completion.

Examples:

```text
alpha<Tab>   -> α
\alpha<Tab>  -> α
rho<Tab>     -> ρ
\rho<Tab>    -> ρ
Delta<Tab>   -> Δ
\Delta<Tab>  -> Δ
partial<Tab> -> ∂
leq<Tab>     -> ≤
approx<Tab>  -> ≈
infty<Tab>   -> ∞
```

If `Tab` does not expand snippets, use `Ctrl+Space` and pick the completion, or set:

```json
"editor.tabCompletion": "on"
```

## Notes

These snippets are scoped to:

```text
markdown, quarto
```

That keeps them out of normal code files.
