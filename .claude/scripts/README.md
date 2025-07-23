# CHANGELOG Management Scripts

This directory contains scripts for maintaining the CHANGELOG.md file.

## Scripts

### `changelog_restructure.py`
Main script for restructuring and maintaining CHANGELOG.md:
- Parses CHANGELOG to structured JSON format
- Removes duplicate entries
- Sorts entries in reverse chronological order
- Groups entries by category within each date
- Creates backup before making changes

**Usage:**
```bash
cd /path/to/SUEWS
python3 .claude/scripts/changelog_restructure.py
```

The script will:
1. Parse the current CHANGELOG.md
2. Clean any duplicate entries
3. Sort all entries by date (newest first)
4. Group entries by category within each date
5. Create a backup file (CHANGELOG.md.backup)
6. Save the restructured CHANGELOG.md
7. Generate a debug JSON file (changelog_debug.json) for inspection

## CHANGELOG Categories
- `[feature]`: New user-facing functionality
- `[bugfix]`: Bug fixes (create GitHub issue)
- `[change]`: User-facing behaviour changes
- `[maintenance]`: Codebase maintenance (including Claude Code/dev tooling)
- `[doc]`: Documentation updates

## Best Practices
1. Run `changelog_restructure.py` after adding new entries to ensure proper sorting
2. For analyzing git commits and generating new entries, use the `/log-changes` slash command
3. The script automatically creates backups before making changes
4. Clean up temporary files (*.backup, *.json) after verifying changes