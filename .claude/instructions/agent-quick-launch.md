# Quick Agent Launch Commands

## ⚠️ CRITICAL: Environment Isolation Required

**EVERY worktree MUST have its own conda environment!**
- Root directory: uses `suews-dev`
- Each worktree: uses `suews-dev-{feature-name}`
- NEVER use `suews-dev` in a worktree!

## Copy-Paste Commands for Each Agent

### Agent 1: Core Runtime Fixes (HIGH PRIORITY)
```bash
mamba create -n suews-core-fixes --clone suews-dev -y
mamba activate suews-core-fixes
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/core-bugs
claude --no-chat
```
**First prompt**: "Read ../../.claude/worktree-plans/feature-core-runtime-fixes.md and start working on the first uncompleted task. Run make dev first, then begin with issue #391."

### Agent 2: Adjust Default Values (HIGH PRIORITY)
```bash
mamba create -n suews-defaults --clone suews-dev -y
mamba activate suews-defaults
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/default-values
claude --no-chat
```
**First prompt**: "Read ../../.claude/worktree-plans/feature-adjust-default-values.md and start working on issue #428. Run make dev first, then audit default values in the data model."

### Agent 3: Hide Internal Options (CONTINUE WORK)
```bash
mamba create -n suews-hide-internal --clone suews-dev -y
mamba activate suews-hide-internal
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/hide-internal-options
claude --no-chat
```
**First prompt**: "Check git status, read ../../.claude/worktree-plans/feature-hide-internal-options.md, run make dev, then continue from task 4 (adding metadata markers)."

### Agent 4: RSL Physics Fixes
```bash
mamba create -n suews-rsl-fixes --clone suews-dev -y
mamba activate suews-rsl-fixes
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/rsl-fixes
claude --no-chat
```
**First prompt**: "Read ../../.claude/worktree-plans/feature-rsl-physics-fixes.md and start with issue #419. Run make dev first, then investigate MOST implementation for short buildings."

### Agent 5: SuPy Data Processing
```bash
mamba create -n suews-data-proc --clone suews-dev -y
mamba activate suews-data-proc
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/supy-processing
claude --no-chat
```
**First prompt**: "Read ../../.claude/worktree-plans/feature-supy-data-processing.md and start with issue #408. Run make dev first, then investigate save_supy output issues."

## Monitor All Agents

In a separate terminal:
```bash
# Watch for plan updates
watch -n 30 'cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS && git status .claude/worktree-plans/'

# Check active environments
mamba env list | grep suews-

# See recent commits
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS
git log --oneline --all --since="1 hour ago"
```

## Emergency Stop All Agents
```bash
# List all suews environments
mamba env list | grep suews- | grep -v suews-dev

# Remove all agent environments (keep suews-dev)
for env in suews-core-fixes suews-defaults suews-hide-internal suews-rsl-fixes suews-data-proc; do
    mamba env remove -n $env -y 2>/dev/null
done
```

## Key Reminders for Agents
1. Always run `make dev` before starting work
2. Run `make test` after each task
3. Update plan progress immediately
4. Commit only when tests pass
5. Push plan updates to master branch