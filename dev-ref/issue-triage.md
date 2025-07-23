# SUEWS Issue Triage Guide

This guide provides a systematic approach for triaging issues in the SUEWS repository using a MECE (Mutually Exclusive, Collectively Exhaustive) label system.

## Quick Reference

### Label System Overview
Every issue gets labels from these dimensions:
- **1-** prefix = Type (what is it?)
- **2-** prefix = Priority (how urgent?)
- **3-** prefix = Status (current state)
- **No prefix** = Area (what part of system?)

### Required Labels
1. **Type**: `1-bug`, `1-feature`, `1-question`
2. **Priority**: `2-P0`, `2-P1`, `2-P2`
3. **Status**: `3-triage`, `3-ready`, `3-in-progress`, `3-needs-science`, `3-needs-deps`
4. **Area**: `module:*`, `infra:*`, `doc:*`

### Optional Labels
- **Special**: `good-first-issue`, `help-wanted`

## Triage Decision Tree

### Step 1: Issue Type (PICK ONE)
```
Is it reporting a problem?
  └─> 1-bug

Is it requesting new functionality?
  └─> 1-feature

Is it asking a question or unclear about something?
  └─> 1-question
```

### Step 2: Priority (PICK ONE)
```
For Bugs:
├─ Data corruption/crashes? → 2-P0
├─ Core functionality broken? → 2-P1
└─ Minor issue? → 2-P2

For Features/Docs/Questions → 2-P2 (default)
```

### Step 3: Status (PICK ONE)
```
What's the current state?
├─ Needs assessment → 3-triage
├─ Ready for work → 3-ready
├─ Being worked on → 3-in-progress
├─ Needs scientific input → 3-needs-science
└─ Waiting on dependencies → 3-needs-deps
```

### Step 4: Area (PICK ONE)
```
What part of the system?
├─ Project meta? → meta:xxx
│   └─> meta:governance, meta:release
├─ Physics module? → module:xxx
│   └─> module:ohm, module:snow, etc.
├─ Infrastructure? → infra:xxx
│   └─> infra:data-model, infra:ci, infra:build, etc.
└─ Documentation? → doc:xxx
    └─> doc:user, doc:api, doc:dev
```

### Step 5: Special Labels (OPTIONAL)
```
Is this good for newcomers? → good-first-issue
Need community help? → help-wanted
```

## Complete Label Reference

### Type Labels (Pick ONE)
- `1-bug` - Something isn't working
- `1-feature` - New functionality request
- `1-question` - User question/support/unclear

### Priority Labels (Pick ONE)
- `2-P0` - Critical (crashes/data loss)
- `2-P1` - High (major functionality broken)
- `2-P2` - Normal (everything else)

### Status Labels (Pick ONE)
- `3-triage` - Needs assessment
- `3-ready` - Ready for work
- `3-in-progress` - Being worked on
- `3-needs-science` - Requires scientific input/decision
- `3-needs-deps` - Waiting on other issues/PRs

### Area Labels (Pick ONE)

#### Meta/Project
- `meta:governance` - Project governance and process
- `meta:release` - Release planning and management

#### Physics Modules
- `module:anohm` - Analytical OHM
- `module:anthro` - Anthropogenic heat
- `module:atmmoiststab` - Atmospheric stability
- `module:beers` - Building radiation
- `module:biogenco2` - Biogenic CO2
- `module:bluews` - Building water balance
- `module:dailystate` - Daily state updates
- `module:ehc` - Explicit Heat Conduction
- `module:estm` - Element Surface Temperature
- `module:evap` - Evaporation processes
- `module:lumps` - LUMPS energy balance
- `module:narp` - Net All-wave Radiation
- `module:ohm` - Objective Hysteresis Model
- `module:resist` - Resistance calculations
- `module:rslprof` - RSL profiles
- `module:snow` - Snow processes
- `module:solweig` - Solar/longwave radiation
- `module:spartacus` - SPARTACUS radiation
- `module:stebbs` - Surface Temperature Energy Balance
- `module:waterdist` - Water distribution

#### Infrastructure
- `infra:data-model` - Pydantic data models & validation
- `infra:ci` - CI/CD pipelines, GitHub Actions
- `infra:build` - Build system (meson, setup.py, compilation)
- `infra:packaging` - Package distribution, dependencies
- `infra:test` - Testing infrastructure, pytest
- `infra:fortran-python` - Fortran-Python interface
- `infra:logging` - Logging, debugging, and diagnostic output
- `infra:type-safety` - Type safety and compiler compatibility
- `infra:code-refactor` - Code refactoring and cleanup
- `infra:input` - Input data handling and validation
- `infra:output` - Output formatting and export
- `infra:utility` - Utility functions and helpers

#### Documentation
- `doc:user` - User guides, tutorials
- `doc:api` - API documentation
- `doc:dev` - Developer documentation

### Special Labels (Use when helpful)
- `good-first-issue` - Good for newcomers
- `help-wanted` - Need community help

## Triage Examples

### Example 1: Bug Report
**Issue**: "SUEWS crashes when snow module enabled on Windows"

Labels:
- `1-bug`
- `2-P0` (crashes)
- `3-ready`
- `module:snow`

### Example 2: Feature Request
**Issue**: "Add support for Pydantic validation in YAML config"

Labels:
- `1-feature`
- `2-P2`
- `3-ready`
- `infra:data-model`

### Example 3: User Question
**Issue**: "How to configure anthropogenic heat?"

Labels:
- `1-question`
- `2-P2`
- `3-ready`
- `doc:user`

### Example 4: CI Build Failure
**Issue**: "Build fails on macOS ARM64 with Python 3.11"

Labels:
- `1-bug`
- `2-P0` (blocks releases)
- `3-triage`
- `infra:ci`

### Example 5: Science Question
**Issue**: "What albedo values should be used for new green roof type?"

Labels:
- `1-feature`
- `2-P2`
- `3-needs-science`
- `module:ohm`

### Example 6: Dependency Blocked
**Issue**: "Update to Pydantic v3 once numpy compatibility fixed"

Labels:
- `1-feature`
- `2-P2`
- `3-needs-deps`
- `infra:data-model`

## GitHub Projects Integration

Use these filters in GitHub Projects:

```
# View all bugs
label:1-bug

# High priority items ready to work
label:2-P0,2-P1 label:3-ready

# Infrastructure issues
label:infra:pydantic,infra:ci,infra:build

# Physics module bugs
label:1-bug label:module:ohm,module:snow,module:estm

# In progress work
label:3-in-progress

# Good first issues
label:good-first-issue

# Needs triage
label:3-triage

# Needs science input
label:3-needs-science

# Blocked by dependencies
label:3-needs-deps
```

## MECE Label System

### Design Principles
- **Mutually Exclusive**: Each label dimension has clear, non-overlapping categories
- **Collectively Exhaustive**: Every issue can be properly categorized
- **Scalable**: Easy to add new areas as the project grows

### Label Requirements

**Every issue must have labels from these dimensions**:
1. Type: `1-bug`, `1-feature`, or `1-question`
2. Priority: `2-P0`, `2-P1`, or `2-P2`
3. Status: `3-triage`, `3-ready`, `3-in-progress`, `3-needs-science`, or `3-needs-deps`
4. Area: `meta:*`, `module:*`, `infra:*`, or `doc:*`

**Optional labels**:
- `good-first-issue` or `help-wanted` when appropriate

### Note on Closed Issues
- No need for "done" status label - GitHub tracks closed issues automatically
- When an issue is completed, simply close it rather than changing labels

## Best Practices

1. **Label immediately** - Every new issue gets labels from all 4 dimensions
2. **Start with triage** - New issues begin at `3-triage` status
3. **Update status** - Change when work starts/stops
4. **Use comments** - Details go in comments, not labels
5. **Close completed work** - Don't use "done" label, just close the issue

## Quick Triage Checklist

- [ ] Set type label (`1-bug`, `1-feature`, `1-question`)
- [ ] Set priority label (`2-P0`, `2-P1`, `2-P2`)
- [ ] Set status label (`3-triage`, `3-ready`, `3-in-progress`, `3-needs-science`, `3-needs-deps`)
- [ ] Set area label (`meta:*`, `module:*`, `infra:*`, `doc:*`)
- [ ] Add `good-first-issue` if appropriate
- [ ] Assign if someone will work on it
- [ ] Link related issues