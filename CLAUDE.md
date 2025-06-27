# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Worktree Context Management

### Branch-Specific Plans
When working in a git worktree or on a specific feature branch, check for branch-specific context and plans:

1. **First, identify current branch:**
   ```bash
   git branch --show-current
   ```

2. **Then load the corresponding plan if it exists:**
   - For branch `feature/hide-internal-options` → See `.claude-plans/feature-hide-internal-options.md`
   - For branch `feature/rsl-physics-fixes` → See `.claude-plans/feature-rsl-physics-fixes.md`
   - For branch `feature/fast-dev-build` → See `.claude-plans/feature-fast-dev-build.md`
   - For other feature branches → Check `.claude-plans/feature-{branch-name}.md`

3. **If no plan exists**, proceed with standard development practices.

### Plan Lifecycle Management

**Creating a new plan:**
- When starting complex multi-session work, create `.claude-plans/feature-{branch-name}.md`
- Include: current context, progress tracking, key decisions, implementation steps
- Commit to master/main branch so it's available in all worktrees

**Maintaining plans:**
- Update progress status in the plan as work proceeds
- Add new findings or decisions that affect implementation
- Keep plans focused and actionable

**Cleaning up completed plans:**
- When a feature branch is merged, remove its plan file
- Archive important decisions to main documentation if needed
- Clean up with: `git rm .claude-plans/feature-{branch-name}.md`

### Example Plan Structure
```markdown
# Feature: [Feature Name]

## Context
Brief description of what this feature/fix addresses

## Progress Tracking
- [ ] Task 1
- [x] Task 2 (completed)
- [ ] Task 3

## Key Decisions
- Decision 1: Reasoning
- Decision 2: Reasoning

## Implementation Notes
Specific implementation details, gotchas, or important context

## Files to Modify
- `path/to/file1.py`
- `path/to/file2.py`
```


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