# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

[... existing content remains unchanged ...]

## Git and GitHub Tips

- When using gh cli, first check remotes, only keep original - remove others

## Style and Language Guidelines

- Any human writing in this project should use British English - docs/code annotations etc

## Testing Resources

### Benchmark Test Files
- For testing: configuration file `p_config = Path("test/benchmark1/benchmark1.yml")` 
- For testing: forcing data file `p_forcing = Path("test/benchmark1/forcing/Kc1_2011_data_5.txt")`

## Documentation Guidelines

- Remember the yaml rst files are generated - so modify the `generate_datamodel_rst.py` script rather than the rst files if edits are needed