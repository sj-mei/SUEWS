/**
 * SUEWS Configuration Builder - UI Module
 * 
 * This module handles UI utilities, event listeners, search, and display functions.
 */

window.configBuilder = window.configBuilder || {};
window.configBuilder.ui = {};

/**
 * Setup event listeners
 */
window.configBuilder.ui.setupEventListeners = function() {
    // Search functionality
    const searchInput = document.getElementById('parameterSearch');
    if (searchInput) {
        searchInput.addEventListener('input', (e) => {
            window.configBuilder.ui.performSearch(e.target.value);
        });
    }
    
    // Theme toggle
    const themeToggle = document.getElementById('themeToggle');
    if (themeToggle) {
        themeToggle.addEventListener('click', () => {
            document.body.classList.toggle('dark-mode');
            const icon = themeToggle.querySelector('i');
            icon.classList.toggle('fa-moon');
            icon.classList.toggle('fa-sun');
        });
    }
    
    // Import/Export buttons
    const importBtn = document.getElementById('importBtn');
    if (importBtn) {
        importBtn.addEventListener('click', window.configBuilder.ui.showImportModal);
    }
    
    const exportBtn = document.getElementById('exportBtn');
    if (exportBtn) {
        exportBtn.addEventListener('click', () => window.configBuilder.ui.exportConfig('yaml'));
    }
    
    const exportJsonBtn = document.getElementById('exportJsonBtn');
    if (exportJsonBtn) {
        exportJsonBtn.addEventListener('click', () => window.configBuilder.ui.exportConfig('json'));
    }
    
    // Export YAML button (in header)
    const exportYamlBtn = document.getElementById('exportYamlBtn');
    if (exportYamlBtn) {
        exportYamlBtn.addEventListener('click', () => window.configBuilder.ui.exportConfig('yaml'));
    }
    
    // Validate button
    const validateBtn = document.getElementById('validateBtn');
    if (validateBtn) {
        validateBtn.addEventListener('click', () => window.configBuilder.ui.validateConfig());
    }
    
    // Reset button
    const resetBtn = document.getElementById('resetBtn');
    if (resetBtn) {
        resetBtn.addEventListener('click', window.configBuilder.ui.resetForm);
    }
    
    // Modal close buttons
    document.querySelectorAll('[data-bs-dismiss="modal"]').forEach(btn => {
        btn.addEventListener('click', () => {
            const modal = btn.closest('.modal');
            if (modal) {
                const bsModal = bootstrap.Modal.getInstance(modal);
                if (bsModal) {
                    bsModal.hide();
                }
            }
        });
    });
};

/**
 * Perform parameter search
 */
window.configBuilder.ui.performSearch = function(searchTerm) {
    const normalizedSearch = searchTerm.toLowerCase().trim();
    const allFields = document.querySelectorAll('.form-field, .card');
    
    if (!normalizedSearch) {
        // Clear search - show all fields
        allFields.forEach(field => {
            field.classList.remove('search-highlight', 'search-no-match');
        });
        window.configBuilder.ui.updateFieldVisibility();
        return;
    }
    
    allFields.forEach(field => {
        const labels = field.querySelectorAll('label, .card-header');
        let hasMatch = false;
        
        labels.forEach(label => {
            const text = label.textContent.toLowerCase();
            if (text.includes(normalizedSearch)) {
                hasMatch = true;
            }
        });
        
        if (hasMatch) {
            field.classList.remove('search-no-match');
            field.classList.add('search-highlight');
            
            // Expand parent sections if needed
            let parent = field.closest('.collapse');
            while (parent) {
                parent.classList.add('show');
                parent = parent.parentElement.closest('.collapse');
            }
        } else {
            field.classList.remove('search-highlight');
            field.classList.add('search-no-match');
        }
    });
};

/**
 * Update field visibility based on mode and search
 */
window.configBuilder.ui.updateFieldVisibility = function() {
    const isAdvancedMode = document.body.classList.contains('advanced-mode');
    const searchActive = document.getElementById('parameterSearch')?.value.trim() !== '';
    
    document.querySelectorAll('.form-field, .card').forEach(field => {
        if (searchActive) {
            // In search mode, hide non-matches
            if (field.classList.contains('search-no-match')) {
                field.style.display = 'none';
            } else {
                field.style.display = '';
            }
        } else {
            // Normal mode - show all fields
            field.style.display = '';
        }
    });
};

/**
 * Display debug message
 */
window.configBuilder.ui.displayDebug = function(message, containerId = 'general-form-container') {
    console.log('Debug:', message);
    const container = document.getElementById(containerId);
    if (container) {
        const debugDiv = document.createElement('div');
        debugDiv.className = 'alert alert-info';
        debugDiv.textContent = `Debug: ${message}`;
        container.appendChild(debugDiv);
        setTimeout(() => debugDiv.remove(), 3000);
    }
};

/**
 * Display error message
 */
window.configBuilder.ui.displayError = function(message, containerId = 'general-form-container') {
    console.error('Error:', message);
    const container = document.getElementById(containerId);
    if (container) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.textContent = `Error: ${message}`;
        container.appendChild(errorDiv);
    }
};

/**
 * Format property name for display
 */
window.configBuilder.ui.formatPropertyName = function(name) {
    // Convert snake_case to Title Case
    return name
        .split('_')
        .map(word => word.charAt(0).toUpperCase() + word.slice(1).toLowerCase())
        .join(' ');
};

/**
 * Format field label with display name and unit
 */
window.configBuilder.ui.formatFieldLabel = function(propKey, propSchema) {
    console.log(`formatFieldLabel for ${propKey}:`, {
        display_name: propSchema.display_name,
        unit: propSchema.unit,
        title: propSchema.title,
        hasAnyOf: !!propSchema.anyOf
    });
    
    // For anyOf schemas, the display_name and unit are at the top level
    let label = propSchema.display_name || propSchema.title || window.configBuilder.ui.formatPropertyName(propKey);
    
    // Don't append units here - they're shown separately now
    
    return label;
};

/**
 * Render unit with KaTeX
 */
window.configBuilder.ui.renderUnit = function(unit, container) {
    if (!unit || unit === 'dimensionless') {
        return;
    }
    
    const unitSpan = document.createElement('span');
    unitSpan.className = 'field-unit text-muted ms-2';
    
    // Try to render as LaTeX if it contains math symbols
    if (unit.includes('^') || unit.includes('_') || unit.includes('\\') || unit.includes(' ')) {
        try {
            // Convert simple notation to LaTeX with proper formatting
            let latexUnit = unit
                .replace(/\^(-?\d+)/g, '^{$1}')
                .replace(/_(\d+)/g, '_{$1}')
                .replace(/m\^2/g, 'm^2')
                .replace(/m\^3/g, 'm^3')
                .replace(/m\^-1/g, 'm^{-1}')
                .replace(/s\^-1/g, 's^{-1}')
                .replace(/W m\^-2/g, 'W\\,m^{-2}')
                .replace(/kg m\^-3/g, 'kg\\,m^{-3}')
                .replace(/m s\^-1/g, 'm\\,s^{-1}')
                .replace(/deg/g, '°');
            
            // Add display mode delimiters
            if (!latexUnit.includes('$')) {
                latexUnit = `\\(${latexUnit}\\)`;
            }
            
            unitSpan.innerHTML = latexUnit;
            
            // Render with KaTeX if available
            if (typeof renderMathInElement !== 'undefined') {
                renderMathInElement(unitSpan, {
                    delimiters: [
                        {left: '\\(', right: '\\)', display: false},
                        {left: '\\[', right: '\\]', display: true}
                    ],
                    throwOnError: false
                });
            }
        } catch (e) {
            // Fallback to plain text
            unitSpan.textContent = `(${unit})`;
        }
    } else {
        unitSpan.textContent = `(${unit})`;
    }
    
    container.appendChild(unitSpan);
};

/**
 * Show import modal
 */
window.configBuilder.ui.showImportModal = function() {
    const modal = new bootstrap.Modal(document.getElementById('importModal'));
    modal.show();
};

/**
 * Import configuration
 */
window.configBuilder.ui.importConfig = function() {
    const fileInput = document.getElementById('importFile');
    const textInput = document.getElementById('importText');
    
    if (fileInput.files.length > 0) {
        const file = fileInput.files[0];
        const reader = new FileReader();
        
        reader.onload = function(e) {
            try {
                const content = e.target.result;
                let data;
                
                if (file.name.endsWith('.json')) {
                    data = JSON.parse(content);
                } else if (file.name.endsWith('.yml') || file.name.endsWith('.yaml')) {
                    data = jsyaml.load(content);
                } else {
                    throw new Error('Unsupported file format. Please use .json or .yaml/.yml files.');
                }
                
                // Update configData
                window.configBuilderState.configData = data;
                
                // Regenerate form
                window.configBuilder.forms.generateForm();
                
                // Close modal
                const modal = bootstrap.Modal.getInstance(document.getElementById('importModal'));
                modal.hide();
                
                // Clear inputs
                fileInput.value = '';
                textInput.value = '';
                
                alert('Configuration imported successfully!');
            } catch (error) {
                alert('Error importing file: ' + error.message);
            }
        };
        
        reader.readAsText(file);
    } else if (textInput.value.trim()) {
        try {
            let data;
            const text = textInput.value.trim();
            
            // Try to parse as JSON first
            try {
                data = JSON.parse(text);
            } catch {
                // If JSON fails, try YAML
                data = jsyaml.load(text);
            }
            
            // Update configData
            window.configBuilderState.configData = data;
            
            // Regenerate form
            window.configBuilder.forms.generateForm();
            
            // Close modal
            const modal = bootstrap.Modal.getInstance(document.getElementById('importModal'));
            modal.hide();
            
            // Clear inputs
            fileInput.value = '';
            textInput.value = '';
            
            alert('Configuration imported successfully!');
        } catch (error) {
            alert('Error parsing configuration: ' + error.message);
        }
    } else {
        alert('Please select a file or paste configuration text.');
    }
};

/**
 * Export configuration
 */
window.configBuilder.ui.exportConfig = function(format) {
    try {
        const configData = window.configBuilderState.configData;
        
        // Clean up FlexibleRefValues for export
        const cleanedData = window.configBuilder.ui.cleanFlexibleRefValuesForExport(
            JSON.parse(JSON.stringify(configData)), 
            window.configBuilderState.schema
        );
        
        let content, filename, mimeType;
        
        if (format === 'json') {
            content = JSON.stringify(cleanedData, null, 2);
            filename = 'suews-config.json';
            mimeType = 'application/json';
        } else {
            content = jsyaml.dump(cleanedData, {
                indent: 2,
                noRefs: true,
                sortKeys: false
            });
            filename = 'suews-config.yml';
            mimeType = 'text/yaml';
        }
        
        // Create blob and download
        const blob = new Blob([content], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    } catch (error) {
        alert('Error exporting configuration: ' + error.message);
    }
};

/**
 * Validate configuration
 */
window.configBuilder.ui.validateConfig = function() {
    try {
        const configData = window.configBuilderState.configData;
        const schema = window.configBuilderState.schema;
        
        if (!configData || !schema) {
            alert('No configuration to validate');
            return;
        }
        
        // Initialize AJV validator with better error messages
        const ajv = new Ajv({ 
            allErrors: true,
            strict: false,
            verbose: true,
            jsonPointers: true,
            validateSchema: false  // Disable schema validation to avoid meta-schema issues
        });
        
        // Add formats
        ajv.addFormat('date-time', true);
        ajv.addFormat('uri', true);
        
        // Remove custom keyword - not needed for validation
        
        // Add the schema with definitions
        const validate = ajv.compile(schema);
        
        // Clean the data for validation
        const cleanedData = window.configBuilder.ui.cleanFlexibleRefValuesForExport(
            JSON.parse(JSON.stringify(configData)), 
            schema
        );
        
        console.log('Validating data:', cleanedData);
        console.log('Validating with schema:', schema);
        
        // Validate
        const valid = validate(cleanedData);
        
        if (valid) {
            // Update validation status
            window.updateValidationStatus('valid', 'Configuration is valid');
            alert('Configuration is valid!');
        } else {
            // Process errors for better readability
            const processedErrors = window.configBuilder.ui.processValidationErrors(validate.errors, cleanedData, schema);
            
            // Format errors for display
            let errorMessage = 'Validation errors:\n\n';
            
            // Group errors by section
            const errorsBySection = new Map();
            
            processedErrors.forEach(error => {
                const section = error.section || 'General';
                if (!errorsBySection.has(section)) {
                    errorsBySection.set(section, []);
                }
                errorsBySection.get(section).push(error);
            });
            
            // Display errors by section
            errorsBySection.forEach((errors, section) => {
                errorMessage += `${section}:\n`;
                errors.forEach(error => {
                    errorMessage += `  - ${error.field}: ${error.message}\n`;
                });
                errorMessage += '\n';
            });
            
            window.updateValidationStatus('invalid', 'Configuration has errors');
            
            // Show in console for debugging
            console.error('Validation errors:', validate.errors);
            console.error('Processed errors:', processedErrors);
            
            alert(errorMessage);
        }
    } catch (error) {
        console.error('Validation error:', error);
        alert('Error validating configuration: ' + error.message);
    }
};

/**
 * Process validation errors for better readability
 */
window.configBuilder.ui.processValidationErrors = function(errors, data, schema) {
    const processedErrors = [];
    const errorPaths = new Set();
    
    // First pass: collect all error paths to filter out parent anyOf errors
    errors.forEach(err => {
        if (err.instancePath) {
            errorPaths.add(err.instancePath);
        }
    });
    
    errors.forEach(err => {
        // Skip generic anyOf errors if we have more specific errors for child paths
        if (err.keyword === 'anyOf') {
            let hasChildErrors = false;
            errorPaths.forEach(path => {
                if (path.startsWith(err.instancePath) && path !== err.instancePath) {
                    hasChildErrors = true;
                }
            });
            if (hasChildErrors) {
                return;
            }
        }
        
        // Extract meaningful path components
        const pathParts = err.instancePath.split('/').filter(p => p);
        let section = 'General';
        let fieldPath = '';
        
        if (pathParts.length > 0) {
            // Determine section based on top-level property
            const topLevel = pathParts[0];
            switch (topLevel) {
                case 'model':
                    section = 'Model Configuration';
                    break;
                case 'sites':
                    section = 'Site Information';
                    if (pathParts[1] === '0') {
                        pathParts[1] = 'Site 1';
                    }
                    break;
                case 'name':
                case 'description':
                    section = 'General Settings';
                    break;
                default:
                    section = topLevel.charAt(0).toUpperCase() + topLevel.slice(1);
            }
            
            // Build human-readable field path
            fieldPath = pathParts.slice(1).map(part => {
                // Convert array indices to readable format
                if (/^\d+$/.test(part)) {
                    return `[${parseInt(part) + 1}]`;
                }
                // Convert snake_case to readable format
                return part.split('_').map(word => 
                    word.charAt(0).toUpperCase() + word.slice(1)
                ).join(' ');
            }).join(' > ');
            
            if (!fieldPath && pathParts.length === 1) {
                fieldPath = pathParts[0].split('_').map(word => 
                    word.charAt(0).toUpperCase() + word.slice(1)
                ).join(' ');
            }
        }
        
        // Create readable error message
        let message = err.message;
        
        if (err.keyword === 'required') {
            const missingField = err.params.missingProperty;
            const readableField = missingField.split('_').map(word => 
                word.charAt(0).toUpperCase() + word.slice(1)
            ).join(' ');
            message = `Missing required field: ${readableField}`;
            fieldPath = fieldPath || 'Configuration';
        } else if (err.keyword === 'type') {
            const expectedType = err.params.type;
            const actualType = err.data === null ? 'null' : 
                             err.data === undefined ? 'undefined' : 
                             Array.isArray(err.data) ? 'array' : 
                             typeof err.data;
            message = `Expected ${expectedType} but got ${actualType}`;
        } else if (err.keyword === 'additionalProperties') {
            message = `Unknown property: ${err.params.additionalProperty}`;
        } else if (err.keyword === 'enum') {
            const allowedValues = err.params.allowedValues;
            message = `Must be one of: ${allowedValues.join(', ')}`;
        } else if (err.keyword === 'minimum' || err.keyword === 'maximum') {
            message = `Value ${err.data} is out of range (${err.message})`;
        } else if (err.keyword === 'minItems' || err.keyword === 'maxItems') {
            message = err.message.charAt(0).toUpperCase() + err.message.slice(1);
        } else if (err.keyword === 'anyOf') {
            // For anyOf errors, try to provide more context
            if (err.schemaPath.includes('sites') && err.schemaPath.includes('properties')) {
                message = 'Invalid value for this field';
            } else {
                message = 'Value does not match expected format';
            }
        }
        
        processedErrors.push({
            section: section,
            field: fieldPath || 'Root',
            message: message,
            path: err.instancePath,
            keyword: err.keyword
        });
    });
    
    // Sort errors by section and field
    processedErrors.sort((a, b) => {
        if (a.section !== b.section) {
            return a.section.localeCompare(b.section);
        }
        return a.field.localeCompare(b.field);
    });
    
    // Remove duplicates
    const uniqueErrors = [];
    const seenErrors = new Set();
    
    processedErrors.forEach(error => {
        const key = `${error.section}-${error.field}-${error.message}`;
        if (!seenErrors.has(key)) {
            seenErrors.add(key);
            uniqueErrors.push(error);
        }
    });
    
    return uniqueErrors;
};

/**
 * Clean FlexibleRefValues for export
 */
window.configBuilder.ui.cleanFlexibleRefValuesForExport = function(data, schemaObj) {
    if (!schemaObj || !schemaObj.properties) {
        return data;
    }
    
    for (const [key, value] of Object.entries(data)) {
        if (value === null || value === undefined) {
            continue;
        }
        
        const propSchema = schemaObj.properties[key];
        if (!propSchema) {
            continue;
        }
        
        // Check if this is a FlexibleRefValue
        if (propSchema.anyOf && propSchema.anyOf.length === 2) {
            const hasRefValue = propSchema.anyOf.some(option => 
                option.$ref && option.$ref.includes('RefValue'));
            const hasRawType = propSchema.anyOf.some(option => 
                option.type && ['number', 'integer', 'string', 'boolean'].includes(option.type));
            
            if (hasRefValue && hasRawType) {
                // This is a FlexibleRefValue
                if (typeof value === 'object' && value.value !== undefined) {
                    // Convert RefValue to raw value
                    data[key] = value.value;
                }
                // If it's already a raw value, keep it as is
            }
        } else if (typeof value === 'object' && !Array.isArray(value)) {
            // Recursively clean nested objects
            let nestedSchema = propSchema;
            
            // Handle $ref
            if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
                const refPath = propSchema.$ref.replace('#/$defs/', '');
                const schema = window.configBuilderState.schema;
                if (schema.$defs && schema.$defs[refPath]) {
                    nestedSchema = schema.$defs[refPath];
                }
            }
            
            data[key] = window.configBuilder.ui.cleanFlexibleRefValuesForExport(value, nestedSchema);
        } else if (Array.isArray(value) && propSchema.items) {
            // Handle arrays
            let itemsSchema = propSchema.items;
            
            // Handle $ref in items
            if (itemsSchema.$ref && itemsSchema.$ref.startsWith('#/$defs/')) {
                const refPath = itemsSchema.$ref.replace('#/$defs/', '');
                const schema = window.configBuilderState.schema;
                if (schema.$defs && schema.$defs[refPath]) {
                    itemsSchema = schema.$defs[refPath];
                }
            }
            
            // Clean each array item
            data[key] = value.map(item => {
                if (typeof item === 'object' && !Array.isArray(item)) {
                    return window.configBuilder.ui.cleanFlexibleRefValuesForExport(item, itemsSchema);
                }
                return item;
            });
        }
    }
    
    // Remove any __is_copied markers
    if (data.__is_copied) {
        delete data.__is_copied;
    }
    
    return data;
};

/**
 * Reset form
 */
window.configBuilder.ui.resetForm = function() {
    if (confirm('Are you sure you want to reset the form? All changes will be lost.')) {
        // Reset to empty config
        window.configBuilder.schema.initializeEmptyConfig();
        
        // Regenerate form
        window.configBuilder.forms.generateForm();
        
        alert('Form has been reset.');
    }
};

/**
 * Format display value
 */
window.configBuilder.ui.formatDisplayValue = function(value, schema) {
    if (value === null || value === undefined) {
        return 'Not set';
    }
    
    // For FlexibleRefValue, show the actual value
    if (typeof value === 'object' && value.value !== undefined) {
        const displayValue = value.value;
        if (value.ref) {
            return `${displayValue} (Ref: ${value.ref})`;
        }
        return String(displayValue);
    }
    
    // For arrays
    if (Array.isArray(value)) {
        return `Array[${value.length}]`;
    }
    
    // For objects
    if (typeof value === 'object') {
        return 'Object';
    }
    
    return String(value);
};