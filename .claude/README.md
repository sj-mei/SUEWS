# Claude Code Workspace Directory

This directory contains all Claude Code-specific files, including worktree plans, instructions, and workspace documents.

## Directory Structure

```
.claude/
├── worktree-plans/      # Feature branch plans (committed to git)
│   ├── README.md
│   ├── claude-dev-notes.md
│   └── feature-*.md
├── instructions/        # Agent launch instructions and guides
│   ├── agent-launch-instructions.md
│   └── agent-quick-launch.md
├── workspace/           # Analysis and working documents
│   ├── worktree-analysis-*.md
│   └── worktree-build-analysis.md
├── settings.json        # Claude Code settings (committed)
└── settings.local.json  # Local Claude Code settings (ignored)
```

## Directory Purposes

### worktree-plans/
**Status**: Committed to git  
Contains branch-specific plans that provide context across Claude Code sessions. These plans track:
- GitHub issues to address
- Progress on tasks
- Key decisions and findings
- Implementation notes

### instructions/
**Status**: Committed to git  
Contains detailed instructions for:
- Launching multiple Claude Code agents
- Managing parallel development
- Build system isolation strategies
- Testing requirements

### workspace/
**Status**: Committed to git (except temp-* files)  
Contains analysis documents and workspace information:
- Build system analysis
- Worktree status reports
- Temporary working documents (temp-* are ignored)

### settings files
- `settings.json`: Shared Claude Code configuration (committed)
- `settings.local.json`: Local overrides (ignored)

## Quick Reference

### For Claude Code Agents
1. Check your branch: `git branch --show-current`
2. Read your plan: `cat ../../.claude/worktree-plans/feature-{branch-name}.md`
3. Follow instructions: `cat ../../.claude/instructions/agent-quick-launch.md`

### For Developers
- All worktree plans are in `.claude/worktree-plans/`
- Launch instructions are in `.claude/instructions/`
- Build analysis is in `.claude/workspace/worktree-build-analysis.md`

## Git Policy
- ✅ Commit: worktree-plans/, instructions/, workspace/*-analysis.md, settings.json
- ❌ Ignore: settings.local.json, workspace/temp-*