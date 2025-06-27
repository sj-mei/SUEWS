# Claude Code Agent Launch Instructions
## Using Separate Conda Environments for Parallel Development

### Prerequisites
- Ensure you have the base `suews-dev` environment already created
- All worktrees are properly set up under `worktrees/`
- Plans are created in `.claude/worktree-plans/`

### Step-by-Step Launch Instructions

## Agent 1: Core Runtime Fixes
```bash
# Terminal 1
mamba create -n suews-core-fixes --clone suews-dev
mamba activate suews-core-fixes
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/core-bugs

# Launch Claude Code
claude --no-chat

# Agent instructions:
# 1. Run: git branch --show-current
# 2. Read: cat ../../.claude/worktree-plans/feature-core-runtime-fixes.md
# 3. Run: make dev
# 4. Start working on tasks in the plan
# 5. Run: make test (after each task)
# 6. Update plan progress as you work
```

## Agent 2: Adjust Default Values
```bash
# Terminal 2
mamba create -n suews-defaults --clone suews-dev
mamba activate suews-defaults
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/default-values

# Launch Claude Code
claude --no-chat

# Agent instructions:
# 1. Run: git branch --show-current
# 2. Read: cat ../../.claude/worktree-plans/feature-adjust-default-values.md
# 3. Run: make dev
# 4. Start working on removing problematic defaults (#428)
# 5. Run: make test (after each change)
# 6. Update plan progress as you work
```

## Agent 3: Hide Internal Options (Continue Existing Work)
```bash
# Terminal 3
mamba create -n suews-hide-internal --clone suews-dev
mamba activate suews-hide-internal
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/hide-internal-options

# Launch Claude Code
claude --no-chat

# Agent instructions:
# 1. Check uncommitted changes: git status
# 2. Read: cat ../../.claude/worktree-plans/feature-hide-internal-options.md
# 3. Run: make dev
# 4. Continue from task 4 (metadata markers)
# 5. Commit existing work first if needed
# 6. Run: make test (after changes)
```

## Agent 4: RSL Physics Fixes
```bash
# Terminal 4
mamba create -n suews-rsl-fixes --clone suews-dev
mamba activate suews-rsl-fixes
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/rsl-fixes

# Launch Claude Code
claude --no-chat

# Agent instructions:
# 1. Run: git branch --show-current
# 2. Read: cat ../../.claude/worktree-plans/feature-rsl-physics-fixes.md
# 3. Run: make dev
# 4. Focus on #419 (MOST for short buildings)
# 5. Run: make test (especially RSL-related tests)
# 6. Update plan with findings
```

## Agent 5: SuPy Data Processing
```bash
# Terminal 5
mamba create -n suews-data-proc --clone suews-dev
mamba activate suews-data-proc
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS/worktrees/supy-processing

# Launch Claude Code
claude --no-chat

# Agent instructions:
# 1. Run: git branch --show-current
# 2. Read: cat ../../.claude/worktree-plans/feature-supy-data-processing.md
# 3. Run: make dev
# 4. Start with #408 (save_supy output issues)
# 5. Run: make test (focus on I/O tests)
# 6. Update plan progress
```

### Common Instructions for All Agents

```markdown
## Your Development Workflow

1. **Verify Setup**
   - Confirm branch: `git branch --show-current`
   - Read your plan: `cat ../../.claude/worktree-plans/feature-{your-branch}.md`
   - Build: `make dev`

2. **Development Cycle**
   - Pick a task from the plan
   - Make changes
   - Test: `make test` (or specific tests)
   - Update plan progress immediately
   - Commit when tests pass

3. **Testing Commands**
   - Quick: `pytest test/test_supy.py -v`
   - Specific: `pytest test/test_file.py::test_function -v`
   - Full: `make test`
   - Rebuild: `make clean && make dev && make test`

4. **Updating Plans**
   - Edit: `../../.claude/worktree-plans/feature-{your-branch}.md`
   - Mark completed: Change `- [ ]` to `- [x]`
   - Add notes about discoveries/blockers
   - Document any new tasks found

5. **Committing Work**
   ```bash
   # Ensure tests pass first!
   make test
   
   # Stage changes
   git add -A
   
   # Commit with proper message
   git commit -m "fix: resolve validation bug in SuPy init
   
   - Fixed type conversion in validation.py line 234
   - Added test case for edge condition
   - Updated error messages for clarity
   
   Addresses #391"
   
   # Push to remote
   git push origin feature/{your-branch}
   ```

6. **End of Session**
   - Update plan with session summary
   - Note any blocking issues
   - Commit plan updates to master:
     ```bash
     cd ../..  # Go to main repo
     git add .claude/worktree-plans/feature-{your-branch}.md
     git commit -m "docs: update {feature} plan progress"
     git push origin master
     ```

## Important Reminders

- Each agent has its own conda environment - no conflicts!
- Always run tests before committing
- Update plans frequently - other agents might need to know
- Focus on the GitHub issues listed in your plan
- Ask for clarification if a task is unclear
- Document any significant findings in the plan
```

### Monitoring Progress

From the main terminal, you can monitor all agents:

```bash
# Check which agents are active
mamba env list | grep suews-

# View recent commits across all worktrees
cd /Users/tingsun/Dropbox\ \(Personal\)/6.Repos/SUEWS
git worktree list

# Check plan updates
ls -la .claude/worktree-plans/feature-*.md

# Monitor test results
find worktrees -name "*.log" -mtime -1
```

### Cleanup After Work

When an agent completes its work:

```bash
# In the worktree
make clean

# Deactivate and remove environment
mamba deactivate
mamba env remove -n suews-{feature-name}
```

### Troubleshooting

1. **Build Fails**
   - Check gfortran is installed
   - Try: `make clean && make dev`
   - Check error messages for missing dependencies

2. **Tests Fail**
   - Run specific failing test with `-vv` for details
   - Check if you need to rebuild: `make clean && make dev`
   - Look for recent changes that might cause issues

3. **Import Errors**
   - Ensure you're in the right conda environment
   - Verify `make dev` completed successfully
   - Check Python path: `which python`

4. **Git Conflicts**
   - Pull latest master: `git pull origin master`
   - Rebase if needed: `git rebase master`
   - Resolve conflicts carefully

Ready to launch agents! Each agent should work independently without conflicts.