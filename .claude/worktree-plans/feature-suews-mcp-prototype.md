# Feature: SUEWS MCP Prototype

## Context
Create a Model Context Protocol (MCP) server for SUEWS that provides intelligent configuration guidance and result interpretation. This will help users create scientifically sound model configurations and understand their simulation outputs.

## GitHub Issues
- No specific issue yet (prototype development)

## Progress Tracking

### Phase 1: Infrastructure Setup
- [x] Create worktree and dedicated environment
- [ ] Set up basic MCP server structure
- [ ] Install MCP Python SDK
- [ ] Create project structure with tools/resources directories

### Phase 2: Core Integration
- [ ] Build bridge to existing SUEWS pydantic models
- [ ] Implement basic configuration validation tool
- [ ] Test validation with benchmark configuration
- [ ] Create error parsing and reporting system

### Phase 3: Configuration Guidance Tools
- [ ] Implement `validate_config` tool
- [ ] Create `suggest_parameters` tool
- [ ] Build `check_physics_compatibility` tool
- [ ] Develop `generate_config_template` tool
- [ ] Add `explain_parameter` tool

### Phase 4: Result Interpretation Tools
- [ ] Implement `diagnose_energy_balance` tool
- [ ] Create `interpret_thermal_comfort` tool
- [ ] Build `analyze_urban_effects` tool
- [ ] Develop `validate_against_observations` tool
- [ ] Add `generate_insights_report` tool

### Phase 5: Knowledge Base
- [ ] Create parameter documentation database
- [ ] Build physics method compatibility matrix
- [ ] Add typical value ranges by climate
- [ ] Implement interpretation rules engine
- [ ] Include literature reference values

### Phase 6: Testing & Documentation
- [ ] Test with Claude Desktop integration
- [ ] Create example workflows
- [ ] Write user documentation
- [ ] Add performance optimisations

## Key Decisions
- Use Python MCP SDK for implementation
- Leverage existing pydantic validation infrastructure
- Create modular tool architecture for extensibility
- Focus on scientific accuracy and user guidance

## Implementation Notes
- Bridge to SUEWS data models via `src/supy/data_model/`
- Use `ValidationController` for conditional validation logic
- Implement smart suggestions based on urban context
- Provide clear scientific explanations for all guidance

## Files to Create
- `worktrees/suews-mcp/mcp-server/` (main directory)
- `worktrees/suews-mcp/mcp-server/src/server.py` (MCP server)
- `worktrees/suews-mcp/mcp-server/src/tools/` (MCP tools)
- `worktrees/suews-mcp/mcp-server/src/utils/suews_bridge.py` (SUEWS integration)
- `worktrees/suews-mcp/mcp-server/src/knowledge/` (scientific knowledge base)
- `worktrees/suews-mcp/mcp-server/pyproject.toml` (project config)
- `worktrees/suews-mcp/mcp-server/README.md` (documentation)

## Current Status
- Worktree created at `worktrees/suews-mcp`
- Environment `suews-dev-mcp` set up
- Ready to implement MCP server structure