# Feature: FAI Validation Warning

## Context
The current Pydantic validation rule `FAI>0.25*(1-PAI)` is too restrictive for real-world cases, particularly for areas with low PAI (Plan Area Index) such as detached residential houses or industrial areas. This rule needs to be changed from a strict validation error to a warning, allowing users to proceed with their data while being informed of potential issues.

## GitHub Issues
- #326: Pydantic rule FAI>0.25*(1-PAI) unreasonable (PRIMARY - bug)

## Progress Tracking

### Analysis Phase
- [ ] Locate the FAI validation rule in the codebase
- [ ] Understand current validation implementation
- [ ] Identify where warnings should be issued

### Implementation Phase
- [ ] Change validation from error to warning for FAI rule
- [ ] Ensure FAI is set to lower limit when rule is violated
- [ ] Add appropriate warning messages
- [ ] Update any related documentation

### Testing Phase
- [ ] Test with the Gothenburg (Kville) test case mentioned in issue
- [ ] Verify warnings are displayed correctly
- [ ] Ensure FAI is properly forced to lower limit
- [ ] Run full test suite

## Key Decisions
- Convert strict validation error to warning as per Sue Grimmond's comment
- Force FAI to lower limit when validation fails, as suggested by biglimp
- Maintain the rule for informational purposes while allowing processing to continue

## Implementation Notes
- The rule originates from issue #302
- Test case from Gothenburg has 15 grids with 6 reporting FAI value errors
- Common in real-world scenarios with low PAI values
- Warning should be clear but not block processing

## Files to Modify
- Pydantic data model files containing FAI validation
- Validation logic that enforces the FAI>0.25*(1-PAI) rule
- Any warning/logging systems that need to capture this warning
- Related test files

## Current Status
- Worktree created: `worktrees/fai-validation`
- Branch: `feature/fai-validation-warning`
- Environment: `suews-dev-fai-validation`
- Ready to begin implementation