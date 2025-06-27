# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Git Worktrees for Claude Code

This repository uses nested git worktrees to enable parallel development with Claude Code. All worktrees are located under `worktrees/` directory for Claude Code accessibility.

### Worktree Structure
```
SUEWS/
├── worktrees/              # All worktrees nested here (in .gitignore)
│   ├── core-bugs/         # feature/core-runtime-fixes
│   ├── enhancements/      # feature/infrastructure-enhancements
│   ├── fast-dev-build/    # feature/fast-dev-build
│   └── ...
└── .claude/              # Claude Code workspace
    └── worktree-plans/   # Branch-specific context (in master)
        ├── README.md
        ├── feature-core-runtime-fixes.md
        └── ...
```

### Working with Worktrees

**Creating a new worktree with its environment:**
```bash
# Step 1: Create the worktree
git worktree add worktrees/my-feature feature/my-feature-name

# Step 2: Create dedicated environment
mamba create -n suews-dev-my-feature --clone suews-dev

# Step 3: Navigate and activate
cd worktrees/my-feature
mamba activate suews-dev-my-feature

# Step 4: Build in the new environment
make dev
```

**Switching between worktrees:**
```bash
# Step 1: Deactivate current environment (if any)
mamba deactivate

# Step 2: Navigate to worktree
cd worktrees/my-feature

# Step 3: Activate matching environment
mamba activate suews-dev-my-feature

# Step 4: Verify correct environment
echo "Active env: $CONDA_DEFAULT_ENV"  # Should show: suews-dev-my-feature
```

**Listing worktrees:**
```bash
git worktree list
```

**Removing a worktree (including environment):**
```bash
# Step 1: Remove the worktree
git worktree remove worktrees/my-feature

# Step 2: Remove the associated environment
mamba env remove -n suews-dev-my-feature

# Step 3: Verify cleanup
git worktree list
mamba env list | grep suews-dev-my-feature  # Should return nothing
```

### Best Practices
- Always create worktrees under `worktrees/` directory
- Use descriptive names matching the feature
- **Each worktree gets its own `suews-dev-{name}` environment**
- Clean up BOTH worktree AND environment after merging
- Never share environments between worktrees
- Pull master in each worktree to access latest `.claude/worktree-plans/`

### Build System and Testing

**CRITICAL**: Each worktree MUST use a separate conda/mamba environment. Never use the base `suews-dev` environment in worktrees!

#### Environment Setup Rules

1. **Root Directory** (`/SUEWS/`): Uses base `suews-dev` environment
2. **Each Worktree** (`/worktrees/{name}/`): MUST use its own environment

#### Creating Worktree Environments

**For single worktree development:**
```bash
# Create environment for your worktree
mamba create -n suews-dev-{feature-name} --clone suews-dev
cd worktrees/{feature-name}
mamba activate suews-dev-{feature-name}
make dev        # This installs to the worktree-specific environment
make test       # Run tests in isolated environment
```

**Environment Naming Convention:**
- Base environment: `suews-dev` (root directory only)
- Worktree environments: `suews-dev-{feature-name}`
- Examples: `suews-dev-core-fixes`, `suews-dev-rsl-physics`

#### Why Separate Environments Are Required

- `make dev` creates an editable install linked to the current directory
- Only ONE editable install can exist per Python environment
- Multiple worktrees sharing an environment will conflict
- Each worktree may have different code versions

#### Testing Requirements
- **After each task**: Run `make test`
- **Before commits**: Run full test suite
- **For Fortran changes**: `make clean && make dev && make test`
- **Quick tests**: `pytest test/test_supy.py -v`

See:
- `.claude/workspace/worktree-build-analysis.md` for build isolation strategies
- `.claude/workspace/environment-isolation-guide.md` for environment setup

### Current Development Status
- For an overview of all active branches and their associated GitHub issues, see `.claude/worktree-plans/claude-dev-notes.md`
- For agent launch instructions, see `.claude/instructions/`

## Worktree Context Management

### Branch-Specific Plans
When working in a git worktree or on a specific feature branch, check for branch-specific context and plans:

1. **First, identify current branch:**
   ```bash
   git branch --show-current
   ```

2. **Then load the corresponding plan if it exists:**
   - For branch `feature/hide-internal-options` → See `.claude/worktree-plans/feature-hide-internal-options.md`
   - For branch `feature/rsl-physics-fixes` → See `.claude/worktree-plans/feature-rsl-physics-fixes.md`
   - For branch `feature/fast-dev-build` → See `.claude/worktree-plans/feature-fast-dev-build.md`
   - For other feature branches → Check `.claude/worktree-plans/feature-{branch-name}.md`
   
   **Note**: When in a worktree, plans are in the parent directory:
   ```bash
   # From worktree, access plans with:
   cat ../../.claude/worktree-plans/feature-{branch-name}.md
   # Or navigate to parent:
   cd ../.. && cat .claude/worktree-plans/feature-{branch-name}.md
   ```

3. **If no plan exists**, proceed with standard development practices.

### IMPORTANT: Updating Plans During Work

**When working in a worktree, Claude Code MUST:**
1. **Update progress tracking** - Mark tasks as completed immediately when done
2. **Add implementation notes** - Document any important discoveries or decisions
3. **Note blocking issues** - Record any problems encountered
4. **Update file lists** - Add any newly identified files to modify
5. **Record completion status** - Note what was accomplished in each session

**Example update during work:**
```bash
# After completing a task, update the plan:
# Edit .claude/worktree-plans/feature-{branch-name}.md and:
# - Change "- [ ] Fix validation bug" to "- [x] Fix validation bug"
# - Add notes like "Found issue in line 234 of validation.py"
# - Document any new tasks discovered
```

**End of session checklist:**
- [ ] Update all completed tasks in the plan
- [ ] Add notes about any unfinished work
- [ ] Document any blocking issues for next session
- [ ] Commit plan updates to master branch

### Development and Testing Workflow

**Every Claude Code session should follow this workflow:**

1. **Start of Session**
   ```bash
   cd worktrees/{feature-name}
   git branch --show-current  # Verify branch
   cat ../../.claude/worktree-plans/feature-{branch-name}.md  # Read plan
   make dev  # Ensure build is ready
   ```

2. **During Development**
   - After Python changes: `make test`
   - After Fortran changes: `make clean && make dev && make test`
   - Update plan progress immediately when tasks complete

3. **Before Committing**
   ```bash
   make clean && make dev && make test  # Full rebuild and test
   # Only commit if all tests pass!
   ```

4. **Commit Message Format**
   ```bash
   git commit -m "type: brief description

   - Detailed change 1
   - Detailed change 2
   
   Addresses #issue-number"
   ```
   
   Types: feat, fix, docs, test, refactor, chore

### Plan Lifecycle Management

**Creating a new plan:**
- When starting complex multi-session work, create `.claude/worktree-plans/feature-{branch-name}.md`
- Include: current context, progress tracking, key decisions, implementation steps
- Commit to master/main branch so it's available in all worktrees

**Maintaining plans:**
- Update progress status in the plan as work proceeds
- Add new findings or decisions that affect implementation
- Keep plans focused and actionable

**Cleaning up completed plans:**
- When a feature branch is merged, remove its plan file
- Archive important decisions to main documentation if needed
- Clean up with: `git rm .claude/worktree-plans/feature-{branch-name}.md`

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

### Worktree Plan Writing Guide

When creating a plan for Claude Code, follow these guidelines to ensure consistency and clarity:

#### 1. **Context Section**
- Provide a clear, concise description of the feature/fix purpose
- Explain why this work is needed
- Keep it to 2-3 sentences

#### 2. **GitHub Issues Section**
- List all related GitHub issues with their numbers
- Mark the PRIMARY issue if there are multiple
- Include issue labels in parentheses (e.g., "bug", "enhancement")
- Note if issues are CLOSED but need verification

#### 3. **Progress Tracking**
- Use checkbox format for easy visual tracking
- Group related tasks under main headings
- Be specific and actionable (avoid vague tasks)
- Mark completed items immediately when done

#### 4. **Key Decisions**
- Document architectural or design decisions
- Include rationale for each decision
- Note any trade-offs considered
- Keep for future reference

#### 5. **Implementation Notes**
- Technical details that affect implementation
- Known gotchas or edge cases
- Dependencies on other work
- Performance considerations

#### 6. **Files to Modify**
- List specific files that will be changed
- Include brief notes about what changes are needed
- Group by component or subsystem
- Add "(create)" for new files

#### 7. **Additional Sections (as needed)**
- **Testing Strategy**: For complex features
- **Migration Guide**: For breaking changes
- **Performance Goals**: For optimisation work
- **Physics Background**: For scientific features
- **Current Status**: For ongoing work

#### Best Practices
- Keep plans focused and actionable
- Update progress regularly during development
- Remove completed plans after merging
- Reference specific functions/classes when possible
- Include links to relevant documentation
- Note any blocking issues or dependencies


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