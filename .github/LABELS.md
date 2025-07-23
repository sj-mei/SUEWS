# SUEWS Issue Labeling System

This document describes the issue labeling system used in the SUEWS project to effectively triage, prioritize, and categorize issues.

## Overview

SUEWS uses a three-tier labeling system:
1. **Type** - What kind of issue is this?
2. **Priority** - How urgent is this?
3. **Status** - What's the current state?

Additionally, we use category labels to identify which part of the codebase is affected.

## Tier 1: Issue Type (1-*)

Identifies the nature of the issue:

- `1-bug` - Something isn't working
- `1-feature` - New functionality request
- `1-question` - User question or support request

**Usage**: Every issue should have exactly ONE type label.

## Tier 2: Priority (2-P*)

Indicates urgency and importance:

- `2-P0` - Critical priority (system breaking, data corruption, security)
- `2-P1` - High priority (major functionality broken, blocking users)
- `2-P2` - Medium priority (minor issues, nice-to-have features)

**Usage**: Assigned during triage based on impact and urgency.

## Tier 3: Status (3-*)

Tracks the current state of the issue:

- `3-triage` - Needs assessment and understanding
- `3-ready` - Understood and ready to be worked on
- `3-in-progress` - Someone is actively working on it
- `3-needs-science` - Requires scientific input or decision
- `3-needs-deps` - Waiting on other issues/PRs to complete first

**Usage**: Updated as the issue progresses through its lifecycle.

### Status Workflow

```
New Issue → 3-triage → 3-ready → 3-in-progress → Closed
                   ↓                      ↓
              3-needs-science        3-needs-deps
                   ↓                      ↓
                3-ready               (resolved)
```

## Category Labels

These identify which part of the codebase is affected:

### Documentation (doc:*)
- `doc:api` - API documentation
- `doc:dev` - Developer documentation
- `doc:user` - User guides and tutorials

### Infrastructure (infra:*)
- `infra:build` - Build system (meson, setup.py)
- `infra:ci` - CI/CD pipelines
- `infra:code-refactor` - Code refactoring
- `infra:data-model` - Pydantic data models
- `infra:fortran-python` - Fortran-Python interface
- `infra:input` - Input data handling
- `infra:logging` - Logging and diagnostics
- `infra:output` - Output formatting
- `infra:packaging` - Package distribution
- `infra:test` - Testing infrastructure
- `infra:type-safety` - Type safety and compatibility
- `infra:utility` - Utility functions

### Physics Modules (module:*)
Each label corresponds to a specific SUEWS physics module (e.g., `module:evap` for evaporation processes).

### Meta (meta:*)
- `meta:governance` - Project governance
- `meta:release` - Release planning

### Special Labels
- `good-first-issue` - Good for newcomers
- `help-wanted` - Extra attention needed

## Labeling Guidelines

1. **New Issues**: Start with `3-triage` until assessed
2. **During Triage**: 
   - Assign appropriate type (1-*)
   - Set priority (2-P*)
   - Identify affected categories
   - Update status to `3-ready` or `3-needs-science`/`3-needs-deps`
3. **When Work Begins**: Update to `3-in-progress`
4. **If Blocked**: Update to appropriate `3-needs-*` label

## Examples

### Example 1: Bug Report
```
Labels: 1-bug, 2-P1, 3-ready, module:evap, infra:fortran-python
```
A high-priority bug in the evaporation module affecting the Fortran-Python interface, ready to be worked on.

### Example 2: Feature Request Needing Science Input
```
Labels: 1-feature, 2-P2, 3-needs-science, module:anthro
```
A medium-priority feature for anthropogenic heat that requires scientific decision-making.

### Example 3: Documentation Question
```
Labels: 1-question, 3-triage, doc:user
```
A user question about documentation that needs initial assessment.

## Label Management

- Labels are managed by maintainers
- Community members can suggest label changes in comments
- Regular label cleanup to remove unused labels
- See `.github/workflows/README.md` for automation details