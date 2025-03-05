document.addEventListener('DOMContentLoaded', function () {
    console.log('DOM loaded, initializing application...');

    // Check for required elements
    const requiredElements = [
        'importModal', 'validationModal', 'preview-container',
        'importBtn', 'exportJsonBtn', 'exportYamlBtn', 'newBtn', 'validateBtn',
        'confirmImport', 'jsonPreview', 'yamlPreview'
    ];

    const missingElements = [];
    requiredElements.forEach(id => {
        if (!document.getElementById(id)) {
            missingElements.push(id);
        }
    });

    if (missingElements.length > 0) {
        console.error('Missing required elements:', missingElements);
        // Create an error message at the top of the page
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.innerHTML = `<strong>Error:</strong> Missing required HTML elements: ${missingElements.join(', ')}`;
        document.body.insertBefore(errorDiv, document.body.firstChild);
    }

    // Initialize Bootstrap modals if they exist
    const importModalEl = document.getElementById('importModal');
    const validationModalEl = document.getElementById('validationModal');

    if (importModalEl) {
        importModal = new bootstrap.Modal(importModalEl);
    } else {
        console.error('Import modal element not found');
    }

    if (validationModalEl) {
        validationModal = new bootstrap.Modal(validationModalEl);
    } else {
        console.error('Validation modal element not found');
    }

    // Initialize Ajv validator if available
    try {
        if (typeof Ajv !== 'undefined') {
            ajv = new Ajv({
                allErrors: true,
                verbose: true,
                jsonPointers: true
            });
            console.log('Ajv initialized successfully');
        } else {
            console.error('Ajv library not loaded');
        }
    } catch (error) {
        console.error('Error initializing Ajv:', error);
    }

    // Set up event listeners for buttons if they exist
    setupEventListeners();

    // Load the schema
    loadSchema();

    // Manually activate the first tab
    const firstTab = document.querySelector('.nav-tabs .nav-link');
    if (firstTab) {
        displayDebug('Activating first tab: ' + firstTab.textContent);
        firstTab.click();
    } else {
        displayDebug('No tabs found to activate');
    }

    // Add toggle layout functionality
    const toggleLayoutBtn = document.getElementById('toggleLayoutBtn');
    if (toggleLayoutBtn) {
        toggleLayoutBtn.addEventListener('click', function() {
            const mainContainer = document.querySelector('.main-container');
            if (mainContainer) {
                mainContainer.classList.toggle('split-layout');

                // Update button icon and text
                const icon = this.querySelector('i');
                if (icon) {
                    if (mainContainer.classList.contains('split-layout')) {
                        icon.className = 'fas fa-columns';
                        this.innerHTML = '<i class="fas fa-columns"></i> Horizontal Layout';
                    } else {
                        icon.className = 'fas fa-stream';
                        this.innerHTML = '<i class="fas fa-stream"></i> Vertical Layout';
                    }
                }
            }
        });

        // Set initial button text
        toggleLayoutBtn.innerHTML = '<i class="fas fa-stream"></i> Vertical Layout';
    }
});

// Add a new function to set up event listeners
function setupEventListeners() {
    // Import button
    const importBtn = document.getElementById('importBtn');
    if (importBtn) {
        importBtn.addEventListener('click', showImportModal);
    }

    // Export buttons
    const exportJsonBtn = document.getElementById('exportJsonBtn');
    if (exportJsonBtn) {
        exportJsonBtn.addEventListener('click', () => exportConfig('json'));
    }

    const exportYamlBtn = document.getElementById('exportYamlBtn');
    if (exportYamlBtn) {
        exportYamlBtn.addEventListener('click', () => exportConfig('yaml'));
    }

    // New button
    const newBtn = document.getElementById('newBtn');
    if (newBtn) {
        newBtn.addEventListener('click', resetForm);
    }

    // Validate button
    const validateBtn = document.getElementById('validateBtn');
    if (validateBtn) {
        validateBtn.addEventListener('click', validateConfig);
    }

    // Import confirm button
    const confirmImport = document.getElementById('confirmImport');
    if (confirmImport) {
        confirmImport.addEventListener('click', importConfig);
    }

    // Preview format radio buttons
    const jsonPreview = document.getElementById('jsonPreview');
    if (jsonPreview) {
        jsonPreview.addEventListener('change', updatePreview);
    }

    const yamlPreview = document.getElementById('yamlPreview');
    if (yamlPreview) {
        yamlPreview.addEventListener('change', updatePreview);
    }
}

// Add this function to display debug information on the page
function displayDebug(message, containerId = 'general-form-container') {
    console.log(message);
    const container = document.getElementById(containerId);
    if (container) {
        const debugDiv = document.createElement('div');
        debugDiv.className = 'alert alert-info';
        debugDiv.textContent = message;
        container.appendChild(debugDiv);
    }
}

// Update the loadSchema function
function loadSchema() {
    displayDebug('Loading schema...');
    showLoading();
    fetch('suews-config-schema.json')
        .then(response => {
            displayDebug(`Schema fetch response: ${response.status}`);
            if (!response.ok) {
                throw new Error(`Failed to load schema: ${response.status} ${response.statusText}`);
            }
            return response.json();
        })
        .then(schemaData => {
            displayDebug('Schema loaded successfully');
            schema = schemaData;
            displayDebug(`Schema top-level properties: ${Object.keys(schema.properties).join(', ')}`);

            // Initialize empty config
            initializeEmptyConfig();
            displayDebug('Empty config initialized');

            // Generate form
            generateForm();
            displayDebug('Form generation attempted');

            // Update preview
            updatePreview();
            hideLoading();
        })
        .catch(error => {
            console.error('Error loading schema:', error);
            displayDebug(`Error loading schema: ${error.message}`);
            hideLoading();
        });
}

// Initialize empty config with required properties
function initializeEmptyConfig() {
    configData = {
        name: "New SUEWS Configuration",
        description: "Description of this configuration"
    };

    // Add required properties from schema
    if (schema && schema.properties) {
        Object.keys(schema.properties).forEach(key => {
            if (key !== 'name' && key !== 'description') {
                if (schema.properties[key].type === 'array') {
                    configData[key] = [];
                    // Add one empty item for arrays
                    if (schema.properties[key].items && schema.properties[key].items.type === 'object') {
                        const emptyItem = createEmptyObject(schema.properties[key].items);
                        configData[key].push(emptyItem);
                    }
                } else if (schema.properties[key].type === 'object') {
                    configData[key] = createEmptyObject(schema.properties[key]);
                }
            }
        });
    }
}

// Create an empty object based on schema
function createEmptyObject(schemaObj) {
    const obj = {};

    // Handle $ref in schema
    let resolvedSchema = schemaObj;
    if (schemaObj.$ref && schemaObj.$ref.startsWith('#/$defs/')) {
        const refPath = schemaObj.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            resolvedSchema = schema.$defs[refPath];
        } else {
            console.error('Could not resolve schema reference', schemaObj.$ref);
            return obj;
        }
    }

    if (resolvedSchema.properties) {
        Object.keys(resolvedSchema.properties).forEach(propKey => {
            const propSchema = resolvedSchema.properties[propKey];

            // Handle $ref in property schema
            let resolvedPropSchema = propSchema;
            if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
                const refPath = propSchema.$ref.replace('#/$defs/', '');
                if (schema.$defs && schema.$defs[refPath]) {
                    resolvedPropSchema = schema.$defs[refPath];
                } else {
                    console.error('Could not resolve property schema reference', propSchema.$ref);
                }
            }

            if (resolvedPropSchema.type === 'object') {
                obj[propKey] = createEmptyObject(resolvedPropSchema);
            } else if (resolvedPropSchema.type === 'array') {
                obj[propKey] = [];
                // Add one empty item for arrays
                if (resolvedPropSchema.items) {
                    let resolvedItems = resolvedPropSchema.items;
                    if (resolvedPropSchema.items.$ref && resolvedPropSchema.items.$ref.startsWith('#/$defs/')) {
                        const refPath = resolvedPropSchema.items.$ref.replace('#/$defs/', '');
                        if (schema.$defs && schema.$defs[refPath]) {
                            resolvedItems = schema.$defs[refPath];
                        }
                    }

                    if (resolvedItems.type === 'object') {
                        const emptyItem = createEmptyObject(resolvedItems);
                        obj[propKey].push(emptyItem);
                    }
                }
            } else {
                // For ValueWithDOI type
                if (resolvedPropSchema.properties &&
                    resolvedPropSchema.properties.value &&
                    resolvedPropSchema.properties.ref) {
                    obj[propKey] = { value: null, ref: null };
                } else {
                    // For primitive types, use default or null
                    obj[propKey] = resolvedPropSchema.default !== undefined ? resolvedPropSchema.default : null;
                }
            }
        });
    }

    return obj;
}

// Update the generateForm function to ensure model content is generated
function generateForm() {
    console.log('Starting form generation');

    // Generate general form fields
    generateGeneralFields();

    // Generate model form fields - make sure this is working
    generateModelFields();

    // Generate site form fields
    generateSiteFields();

    // Log completion
    console.log('Form generation completed');
}

// Update the generateModelFields function to use the schema
function generateModelFields() {
    console.log('Generating model fields from schema...');

    // Get the model container
    const modelContainer = document.getElementById('model-form-container');
    if (!modelContainer) {
        console.error('Model form container not found');
        return;
    }

    // Clear the container
    modelContainer.innerHTML = '';

    // Check if schema is available
    if (!schema || !schema.properties || !schema.properties.model) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.textContent = 'Schema for model not found';
        modelContainer.appendChild(errorDiv);
        console.error('Schema for model not found', schema);
        return;
    }

    // Get the model schema - handle $ref if present
    let modelSchema = schema.properties.model;

    // If the model schema is a reference, resolve it
    if (modelSchema.$ref && modelSchema.$ref.startsWith('#/$defs/')) {
        const refPath = modelSchema.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            modelSchema = schema.$defs[refPath];
            console.log('Resolved model schema reference:', refPath);
        } else {
            const errorDiv = document.createElement('div');
            errorDiv.className = 'alert alert-danger';
            errorDiv.textContent = `Could not resolve schema reference: ${modelSchema.$ref}`;
            modelContainer.appendChild(errorDiv);
            console.error('Could not resolve schema reference', modelSchema.$ref);
            return;
        }
    }

    // Initialize model object if it doesn't exist
    if (!configData.model) {
        configData.model = createEmptyObject(modelSchema);
    }

    // Check if model schema has properties
    if (!modelSchema.properties) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.textContent = 'Model schema does not have properties defined';
        modelContainer.appendChild(errorDiv);
        console.error('Model schema does not have properties', modelSchema);
        return;
    }

    // Get model properties from schema
    const modelProperties = modelSchema.properties;

    // Create tabs for model sections
    const modelTabs = document.createElement('ul');
    modelTabs.className = 'nav nav-tabs';
    modelTabs.id = 'modelTabs';
    modelTabs.setAttribute('role', 'tablist');

    const modelTabContent = document.createElement('div');
    modelTabContent.className = 'tab-content';
    modelTabContent.id = 'modelTabContent';

    // Track if we've set an active tab
    let activeTabSet = false;

    // Create tabs for each top-level property in the model schema
    Object.keys(modelProperties).forEach((propKey, index) => {
        const propSchema = modelProperties[propKey];
        const isActive = !activeTabSet;

        if (isActive) {
            activeTabSet = true;
        }

        // Create tab item
        const tabItem = document.createElement('li');
        tabItem.className = 'nav-item';
        tabItem.setAttribute('role', 'presentation');

        // Create tab button
        const tabButton = document.createElement('button');
        tabButton.className = `nav-link ${isActive ? 'active' : ''}`;
        tabButton.id = `model-${propKey}-tab`;
        tabButton.setAttribute('data-bs-toggle', 'tab');
        tabButton.setAttribute('data-bs-target', `#model-${propKey}`);
        tabButton.setAttribute('type', 'button');
        tabButton.setAttribute('role', 'tab');
        tabButton.textContent = propSchema.title || formatPropertyName(propKey);

        tabItem.appendChild(tabButton);
        modelTabs.appendChild(tabItem);

        // Create tab content
        const tabContent = document.createElement('div');
        tabContent.className = `tab-pane fade ${isActive ? 'show active' : ''}`;
        tabContent.id = `model-${propKey}`;
        tabContent.setAttribute('role', 'tabpanel');

        // Generate fields for this property
        if (propSchema.type === 'object' || (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/'))) {
            // If it's a reference, resolve it first
            let resolvedSchema = propSchema;
            if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
                const refPath = propSchema.$ref.replace('#/$defs/', '');
                if (schema.$defs && schema.$defs[refPath]) {
                    resolvedSchema = schema.$defs[refPath];
                }
            }

            // Initialize the property if it doesn't exist
            if (!configData.model[propKey]) {
                configData.model[propKey] = createEmptyObject(resolvedSchema);
            }

            generateObjectFields(resolvedSchema, configData.model[propKey], tabContent, `model.${propKey}`);
        } else {
            // Handle non-object properties (unlikely in this context)
            generatePrimitiveField(propSchema, configData.model[propKey], tabContent, `model.${propKey}`, propKey);
        }

        modelTabContent.appendChild(tabContent);
    });

    // Add tabs and content to the model container
    modelContainer.appendChild(modelTabs);
    modelContainer.appendChild(modelTabContent);

    console.log('Model fields generation completed');
}

// Helper function to format property names for display
function formatPropertyName(name) {
    return name
        .replace(/([A-Z])/g, ' $1')
        .replace(/^./, str => str.toUpperCase())
        .replace(/_/g, ' ');
}

// Generate fields for an object
function generateObjectFields(objSchema, objData, container, path) {
    console.log(`Generating object fields for ${path}...`);

    // Handle $ref in schema
    let resolvedSchema = objSchema;
    if (objSchema.$ref && objSchema.$ref.startsWith('#/$defs/')) {
        const refPath = objSchema.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            resolvedSchema = schema.$defs[refPath];
            console.log(`Resolved schema reference for ${path}: ${refPath}`);
        } else {
            console.error(`Could not resolve schema reference for ${path}`, objSchema.$ref);
            const errorDiv = document.createElement('div');
            errorDiv.className = 'alert alert-danger';
            errorDiv.textContent = `Could not resolve schema reference: ${objSchema.$ref}`;
            container.appendChild(errorDiv);
            return;
        }
    }

    if (!resolvedSchema.properties) {
        console.error(`No properties found in schema for ${path}`);
        return;
    }

    // Create sections for different property groups
    const sections = {};
    const defaultSection = document.createElement('div');
    defaultSection.className = 'form-section';
    container.appendChild(defaultSection);

    // Process each property in the schema
    Object.keys(resolvedSchema.properties).forEach(propKey => {
        let propSchema = resolvedSchema.properties[propKey];

        // Handle $ref in property schema
        if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
            const refPath = propSchema.$ref.replace('#/$defs/', '');
            if (schema.$defs && schema.$defs[refPath]) {
                propSchema = schema.$defs[refPath];
                console.log(`Resolved property schema reference for ${path}.${propKey}: ${refPath}`);
            } else {
                console.error(`Could not resolve property schema reference for ${path}.${propKey}`, propSchema.$ref);
            }
        }

        // Initialize property data if it doesn't exist
        if (objData[propKey] === undefined) {
            if (propSchema.type === 'object') {
                objData[propKey] = createEmptyObject(propSchema);
            } else if (propSchema.type === 'array') {
                objData[propKey] = [];
            } else {
                objData[propKey] = propSchema.default !== undefined ? propSchema.default : null;
            }
        }

        // Get the section for this property
        let section = defaultSection;
        if (propSchema.section) {
            if (!sections[propSchema.section]) {
                sections[propSchema.section] = document.createElement('div');
                sections[propSchema.section].className = 'form-section';
                const sectionTitle = document.createElement('h4');
                sectionTitle.textContent = propSchema.section;
                sections[propSchema.section].appendChild(sectionTitle);
                container.appendChild(sections[propSchema.section]);
            }
            section = sections[propSchema.section];
        }

        // Generate field based on property type
        const propPath = `${path}.${propKey}`;

        // Special handling for ValueWithDOI type
        if (propSchema.properties &&
            propSchema.properties.value &&
            propSchema.properties.ref) {
            generateValueWithDOIField(propSchema, objData[propKey], section, propPath, propKey);
        } else if (propSchema.type === 'object') {
            // Create a collapsible card for nested objects
            const card = document.createElement('div');
            card.className = 'card mb-3';

            const cardHeader = document.createElement('div');
            cardHeader.className = 'card-header';

            const cardTitle = document.createElement('h5');
            cardTitle.className = 'mb-0';

            const collapseButton = document.createElement('button');
            collapseButton.className = 'btn btn-link';
            collapseButton.setAttribute('data-bs-toggle', 'collapse');
            collapseButton.setAttribute('data-bs-target', `#collapse-${propPath.replace(/\./g, '-')}`);
            collapseButton.textContent = propSchema.title || formatPropertyName(propKey);

            cardTitle.appendChild(collapseButton);
            cardHeader.appendChild(cardTitle);
            card.appendChild(cardHeader);

            const collapseDiv = document.createElement('div');
            collapseDiv.id = `collapse-${propPath.replace(/\./g, '-')}`;
            collapseDiv.className = 'collapse';

            const cardBody = document.createElement('div');
            cardBody.className = 'card-body';

            // Generate fields for nested object
            generateObjectFields(propSchema, objData[propKey], cardBody, propPath);

            collapseDiv.appendChild(cardBody);
            card.appendChild(collapseDiv);
            section.appendChild(card);
        } else if (propSchema.type === 'array') {
            // Create a collapsible card for arrays
            const card = document.createElement('div');
            card.className = 'card mb-3';

            const cardHeader = document.createElement('div');
            cardHeader.className = 'card-header';

            const cardTitle = document.createElement('h5');
            cardTitle.className = 'mb-0';

            const collapseButton = document.createElement('button');
            collapseButton.className = 'btn btn-link';
            collapseButton.setAttribute('data-bs-toggle', 'collapse');
            collapseButton.setAttribute('data-bs-target', `#collapse-${propPath.replace(/\./g, '-')}`);
            collapseButton.textContent = propSchema.title || formatPropertyName(propKey);

            cardTitle.appendChild(collapseButton);
            cardHeader.appendChild(cardTitle);
            card.appendChild(cardHeader);

            const collapseDiv = document.createElement('div');
            collapseDiv.id = `collapse-${propPath.replace(/\./g, '-')}`;
            collapseDiv.className = 'collapse';

            const cardBody = document.createElement('div');
            cardBody.className = 'card-body';

            // Generate fields for array
            generateArrayFields(propSchema, objData[propKey], cardBody, propPath);

            collapseDiv.appendChild(cardBody);
            card.appendChild(collapseDiv);
            section.appendChild(card);
        } else {
            // Generate field for primitive type
            generatePrimitiveField(propSchema, objData[propKey], section, propPath, propKey);
        }
    });
}

// Generate field for ValueWithDOI type
function generateValueWithDOIField(propSchema, propData, container, path, propKey) {
    // Handle $ref in schema
    let resolvedSchema = propSchema;
    if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
        const refPath = propSchema.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            resolvedSchema = schema.$defs[refPath];
            console.log(`Resolved ValueWithDOI schema reference for ${path}: ${refPath}`);
        } else {
            console.error(`Could not resolve ValueWithDOI schema reference for ${path}`, propSchema.$ref);
            return;
        }
    }

    // Initialize data if needed
    if (!propData || typeof propData !== 'object') {
        propData = { value: null, ref: null };
        setNestedProperty(configData, path, propData);
    }

    // Get value schema
    let valueSchema = resolvedSchema.properties && resolvedSchema.properties.value
        ? resolvedSchema.properties.value
        : { type: 'string' };

    // Handle $ref in value schema
    if (valueSchema.$ref && valueSchema.$ref.startsWith('#/$defs/')) {
        const refPath = valueSchema.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            valueSchema = schema.$defs[refPath];
        }
    }

    // Get field type for value
    const fieldType = getInputType(valueSchema);

    // Get field options for enum
    const options = getEnumOptions(valueSchema);

    // Get min/max for number fields
    const min = valueSchema.minimum !== undefined ? valueSchema.minimum : null;
    const max = valueSchema.maximum !== undefined ? valueSchema.maximum : null;

    // Format the label to include both the wrapper type and the actual field name
    const displayLabel = resolvedSchema.title || `ValueWithDOI[${formatPropertyName(propKey)}]`;

    // Get the description from the property schema first, then from the resolved schema
    const description = propSchema.description || resolvedSchema.description || null;

    console.log(`ValueWithDOI field: ${path}, description: ${description}`);

    // Create form field with DOI
    createValueWithDOIField(
        container,
        path.replace(/\./g, '-'),
        displayLabel,
        fieldType,
        propData.value,
        propData.ref,
        description,
        (value) => {
            // Convert value to appropriate type
            const convertedValue = convertValueToType(value, valueSchema.type);

            // Update value
            propData.value = convertedValue;

            // Update preview
            updatePreview();
        },
        (ref) => {
            // Update ref
            propData.ref = ref;

            // Update preview
            updatePreview();
        },
        options,
        min,
        max
    );
}

// Helper function to determine input type from schema
function getInputType(propSchema) {
    if (propSchema.enum) {
        return 'select';
    }

    switch (propSchema.type) {
        case 'string':
            return propSchema.format === 'textarea' ? 'textarea' : 'text';
        case 'number':
        case 'integer':
            return 'number';
        case 'boolean':
            return 'checkbox';
        default:
            return 'text';
    }
}

// Helper function to get enum options
function getEnumOptions(propSchema) {
    if (!propSchema.enum) return null;

    return propSchema.enum.map((value, index) => {
        let text = value.toString();
        if (propSchema.enumNames && propSchema.enumNames[index]) {
            text = propSchema.enumNames[index];
        }
        return { value, text };
    });
}

// Helper function to convert value to the correct type
function convertValueToType(value, type) {
    if (value === '' || value === null || value === undefined) {
        return null;
    }

    switch (type) {
        case 'number':
        case 'integer':
            return parseFloat(value);
        case 'boolean':
            return Boolean(value);
        default:
            return value;
    }
}

// Update the preview
function updatePreview() {
    try {
        const previewContainer = document.getElementById('preview-container');

        // Check if preview container exists
        if (!previewContainer) {
            console.error('Preview container not found in the DOM');
            return; // Exit the function if container doesn't exist
        }

        const previewFormat = document.querySelector('input[name="previewFormat"]:checked')?.value || 'json';

        try {
            if (previewFormat === 'json') {
                previewContainer.innerHTML = `<pre>${JSON.stringify(configData, null, 2)}</pre>`;
            } else {
                // Check if jsyaml is defined
                if (typeof jsyaml === 'undefined') {
                    previewContainer.innerHTML = `<div class="alert alert-warning">YAML library not loaded. Showing JSON instead.</div>
                                                <pre>${JSON.stringify(configData, null, 2)}</pre>`;
                } else {
                    previewContainer.innerHTML = `<pre>${jsyaml.dump(configData)}</pre>`;
                }
            }
        } catch (error) {
            console.error('Error generating preview:', error);
            previewContainer.innerHTML = `<div class="alert alert-danger">Error generating preview: ${error.message}</div>`;
        }
    } catch (outerError) {
        // Catch any errors in the outer function
        console.error('Fatal error in updatePreview:', outerError);
    }
}

// Show import modal
function showImportModal() {
    document.getElementById('importData').value = '';
    importModal.show();
}

// Import configuration
function importConfig() {
    const importData = document.getElementById('importData').value;
    const importFormat = document.querySelector('input[name="importFormat"]:checked').value;

    try {
        let importedConfig;

        if (importFormat === 'json') {
            importedConfig = JSON.parse(importData);
        } else {
            importedConfig = jsyaml.load(importData);
        }

        // Validate imported config against schema
        if (schema) {
            const validate = ajv.compile(schema);
            const valid = validate(importedConfig);

            if (!valid) {
                showValidationErrors(validate.errors);
                return;
            }
        }

        // Update config data
        configData = importedConfig;

        // Regenerate form
        generateForm();

        // Update preview
        updatePreview();

        // Close modal
        importModal.hide();
    } catch (error) {
        alert(`Error importing configuration: ${error.message}`);
    }
}

// Export configuration
function exportConfig(format) {
    try {
        let exportData;
        let mimeType;
        let filename;

        if (format === 'json') {
            exportData = JSON.stringify(configData, null, 2);
            mimeType = 'application/json';
            filename = `${configData.name || 'suews-config'}.json`;
        } else {
            exportData = jsyaml.dump(configData);
            mimeType = 'text/yaml';
            filename = `${configData.name || 'suews-config'}.yaml`;
        }

        // Create download link
        const blob = new Blob([exportData], { type: mimeType });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    } catch (error) {
        alert(`Error exporting configuration: ${error.message}`);
    }
}

// Reset form
function resetForm() {
    if (confirm('Are you sure you want to create a new configuration? All current data will be lost.')) {
        initializeEmptyConfig();
        generateForm();
        updatePreview();
    }
}

// Validate configuration
function validateConfig() {
    if (!schema) {
        alert('Schema not loaded. Cannot validate configuration.');
        return;
    }

    const validate = ajv.compile(schema);
    const valid = validate(configData);

    if (valid) {
        document.getElementById('validationResults').innerHTML =
            '<div class="alert alert-success">Configuration is valid!</div>';
    } else {
        showValidationErrors(validate.errors);
    }

    validationModal.show();
}

// Show validation errors
function showValidationErrors(errors) {
    const resultsContainer = document.getElementById('validationResults');
    resultsContainer.innerHTML = '<div class="alert alert-danger">Configuration has errors:</div>';

    const errorList = document.createElement('ul');
    errorList.className = 'validation-error-list';

    errors.forEach(error => {
        const errorItem = document.createElement('li');
        errorItem.className = 'validation-error';

        // Format the error message
        let message = `${error.instancePath || '/'}: ${error.message}`;

        // Add additional details for specific error types
        if (error.keyword === 'enum') {
            message += `. Allowed values: ${error.params.allowedValues.join(', ')}`;
        } else if (error.keyword === 'required') {
            message += `: ${error.params.missingProperty}`;
        }

        errorItem.textContent = message;

        // Add button to navigate to the field if possible
        const fieldId = error.instancePath.substring(1).replace(/\//g, '.').replace(/\[(\d+)\]/g, '[$1]');
        const field = document.getElementById(fieldId);

        if (field) {
            const navButton = document.createElement('button');
            navButton.className = 'btn btn-sm btn-outline-primary ms-2';
            navButton.textContent = 'Go to field';
            navButton.addEventListener('click', () => {
                validationModal.hide();

                // Navigate to the tab containing the field
                const pathParts = fieldId.split('.');
                let tabId;

                if (pathParts[0] === 'model') {
                    // For model properties
                    if (pathParts[1] === 'control') {
                        tabId = 'model-control-tab';
                    } else if (pathParts[1] === 'physics') {
                        tabId = 'model-physics-tab';
                    } else {
                        tabId = 'model-tab';
                    }
                } else if (pathParts[0] === 'site') {
                    // For site properties
                    tabId = 'site-tab';

                    // Extract site index and property
                    const siteMatch = fieldId.match(/site\[(\d+)\]/);
                    if (siteMatch) {
                        const siteIndex = parseInt(siteMatch[1]);
                        const siteTabId = `site-${siteIndex}-tab`;

                        // Click the site tab
                        const siteTab = document.getElementById(siteTabId);
                        if (siteTab) {
                            siteTab.click();
                        }
                    }
                } else {
                    tabId = 'general-tab';
                }

                // Click the tab
                const tab = document.getElementById(tabId);
                if (tab) {
                    tab.click();
                }

                // Scroll to and highlight the field
                setTimeout(() => {
                    field.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    field.classList.add('highlight-error');

                    // Remove highlight after a delay
                    setTimeout(() => {
                        field.classList.remove('highlight-error');
                    }, 3000);
                }, 500);
            });

            errorItem.appendChild(navButton);
        }

        errorList.appendChild(errorItem);
    });

    resultsContainer.appendChild(errorList);
}

// Show loading overlay
function showLoading() {
    const loadingOverlay = document.createElement('div');
    loadingOverlay.id = 'loadingOverlay';
    loadingOverlay.className = 'loading-overlay';
    loadingOverlay.innerHTML = `
        <div class="spinner-border text-light" role="status">
            <span class="visually-hidden">Loading...</span>
        </div>
        <div class="loading-text">Loading...</div>
    `;

    document.body.appendChild(loadingOverlay);
}

// Hide loading overlay
function hideLoading() {
    const loadingOverlay = document.getElementById('loadingOverlay');
    if (loadingOverlay) {
        loadingOverlay.remove();
    }
}

// Add a new site
function addNewSite() {
    // Create new empty site
    const newSite = createEmptyObject(schema.properties.site.items);

    // Set default name
    newSite.name = `Site ${configData.site.length + 1}`;

    // Add to config data
    configData.site.push(newSite);

    // Regenerate site fields
    const siteContainer = document.getElementById('site-form-container');
    generateArrayFields(schema.properties.site, configData.site, siteContainer, 'site');

    // Update preview
    updatePreview();

    // Activate the new site tab
    const newSiteIndex = configData.site.length - 1;
    const newSiteTab = document.getElementById(`site-${newSiteIndex}-tab`);
    if (newSiteTab) {
        newSiteTab.click();
    }
}

// Remove a site
function removeSite(index) {
    if (configData.site.length <= 1) {
        alert('Cannot remove the last site. At least one site is required.');
        return;
    }

    if (confirm(`Are you sure you want to remove site "${configData.site[index].name || `Site ${index + 1}`}"?`)) {
        // Remove from config data
        configData.site.splice(index, 1);

        // Regenerate site fields
        const siteContainer = document.getElementById('site-form-container');
        generateArrayFields(schema.properties.site, configData.site, siteContainer, 'site');

        // Update preview
        updatePreview();

        // Activate the first site tab
        const firstSiteTab = document.getElementById('site-0-tab');
        if (firstSiteTab) {
            firstSiteTab.click();
        }
    }
}

// Helper function to get a nested property from an object using a path string
function getNestedProperty(obj, path) {
    const pathParts = path.split('.');
    let current = obj;

    for (const part of pathParts) {
        // Handle array indices
        const match = part.match(/(\w+)\[(\d+)\]/);
        if (match) {
            const arrayName = match[1];
            const index = parseInt(match[2]);

            if (!current[arrayName] || !current[arrayName][index]) {
                return undefined;
            }

            current = current[arrayName][index];
        } else {
            if (current[part] === undefined) {
                return undefined;
            }

            current = current[part];
        }
    }

    return current;
}

// Helper function to set a nested property in an object using a path string
function setNestedProperty(obj, path, value) {
    const pathParts = path.split('.');
    let current = obj;

    for (let i = 0; i < pathParts.length - 1; i++) {
        const part = pathParts[i];

        // Handle array indices
        const match = part.match(/(\w+)\[(\d+)\]/);
        if (match) {
            const arrayName = match[1];
            const index = parseInt(match[2]);

            if (!current[arrayName]) {
                current[arrayName] = [];
            }

            if (!current[arrayName][index]) {
                current[arrayName][index] = {};
            }

            current = current[arrayName][index];
        } else {
            if (!current[part]) {
                current[part] = {};
            }

            current = current[part];
        }
    }

    const lastPart = pathParts[pathParts.length - 1];

    // Handle array index in the last part
    const match = lastPart.match(/(\w+)\[(\d+)\]/);
    if (match) {
        const arrayName = match[1];
        const index = parseInt(match[2]);

        if (!current[arrayName]) {
            current[arrayName] = [];
        }

        current[arrayName][index] = value;
    } else {
        current[lastPart] = value;
    }
}

// Helper function to get schema for a specific path
function getSchemaForPath(path) {
    if (!schema) return null;

    const pathParts = path.split('.');
    let currentSchema = schema;

    for (const part of pathParts) {
        // Handle array indices
        const match = part.match(/(\w+)\[(\d+)\]/);
        if (match) {
            const arrayName = match[1];

            if (!currentSchema.properties || !currentSchema.properties[arrayName]) {
                return null;
            }

            currentSchema = currentSchema.properties[arrayName].items;
        } else {
            if (!currentSchema.properties || !currentSchema.properties[part]) {
                return null;
            }

            currentSchema = currentSchema.properties[part];
        }
    }

    return currentSchema;
}

// Initialize tabs when the DOM is loaded
document.addEventListener('DOMContentLoaded', function () {
    // Set up event listeners for main tabs
    const mainTabs = document.querySelectorAll('.nav-tabs .nav-link');
    mainTabs.forEach(tab => {
        tab.addEventListener('click', function () {
            const tabId = this.getAttribute('data-bs-target');
            const tabContent = document.querySelector(tabId);

            // Hide all tab contents
            document.querySelectorAll('.tab-pane').forEach(content => {
                content.classList.remove('show', 'active');
            });

            // Show selected tab content
            tabContent.classList.add('show', 'active');

            // Update active tab
            mainTabs.forEach(t => {
                t.classList.remove('active');
            });
            this.classList.add('active');
        });
    });

    // Activate the first tab by default
    if (mainTabs.length > 0) {
        mainTabs[0].click();
    }
});

// Add the missing generateGeneralFields function
function generateGeneralFields() {
    console.log('Generating general fields...');

    const container = document.getElementById('general-form-container');
    if (!container) {
        console.error('General form container not found');
        return;
    }

    // Clear the container
    container.innerHTML = '';

    // Check if schema is available
    if (!schema || !schema.properties) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.textContent = 'Schema not found or invalid';
        container.appendChild(errorDiv);
        console.error('Schema not found or invalid', schema);
        return;
    }

    // Create name field
    if (schema.properties.name) {
        createFormField(
            container,
            'name',
            'Name',
            'text',
            configData.name || '',
            schema.properties.name.description || 'Name of the SUEWS configuration',
            value => {
                configData.name = value;
                updatePreview();
            }
        );
    }

    // Create description field
    if (schema.properties.description) {
        createFormField(
            container,
            'description',
            'Description',
            'textarea',
            configData.description || '',
            schema.properties.description.description || 'Description of this SUEWS configuration',
            value => {
                configData.description = value;
                updatePreview();
            }
        );
    }

    console.log('General fields generated');
}

// Generate field for primitive value
function generatePrimitiveField(propSchema, propData, container, path, propKey) {
    // Handle $ref in schema
    let resolvedSchema = propSchema;
    if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
        const refPath = propSchema.$ref.replace('#/$defs/', '');
        if (schema.$defs && schema.$defs[refPath]) {
            resolvedSchema = schema.$defs[refPath];
            console.log(`Resolved primitive schema reference for ${path}: ${refPath}`);
        } else {
            console.error(`Could not resolve primitive schema reference for ${path}`, propSchema.$ref);
            return;
        }
    }

    // Get field type
    const fieldType = getInputType(resolvedSchema);

    // Get field options for enum
    const options = getEnumOptions(resolvedSchema);

    // Get min/max for number fields
    const min = resolvedSchema.minimum !== undefined ? resolvedSchema.minimum : null;
    const max = resolvedSchema.maximum !== undefined ? resolvedSchema.maximum : null;

    // Create form field
    createFormField(
        container,
        path.replace(/\./g, '-'),
        resolvedSchema.title || formatPropertyName(propKey),
        fieldType,
        propData,
        resolvedSchema.description,
        (value) => {
            // Convert value to appropriate type
            const convertedValue = convertValueToType(value, resolvedSchema.type);

            // Update data
            setNestedProperty(configData, path, convertedValue);

            // Update preview
            updatePreview();
        },
        options,
        min,
        max
    );
}

// Add the missing generateArrayFields function
function generateArrayFields(arraySchema, arrayData, container, path) {
    console.log(`Generating array fields for ${path}...`);

    // Create array container
    const arrayContainer = document.createElement('div');
    arrayContainer.className = 'array-container';
    container.appendChild(arrayContainer);

    // Add description if available
    if (arraySchema.description) {
        const descriptionDiv = document.createElement('div');
        descriptionDiv.className = 'field-description mb-2';
        descriptionDiv.textContent = arraySchema.description;
        arrayContainer.appendChild(descriptionDiv);
    }

    // Create items container
    const itemsContainer = document.createElement('div');
    itemsContainer.className = 'array-items';
    itemsContainer.id = `${path.replace(/\./g, '-')}-items`;
    arrayContainer.appendChild(itemsContainer);

    // Function to add a new item
    const addItem = () => {
        const itemIndex = arrayData.length;
        const itemPath = `${path}[${itemIndex}]`;

        // Create item container
        const itemDiv = document.createElement('div');
        itemDiv.className = 'array-item card mb-3';
        itemDiv.dataset.index = itemIndex;

        // Create item header
        const itemHeader = document.createElement('div');
        itemHeader.className = 'card-header d-flex justify-content-between align-items-center';

        const itemTitle = document.createElement('h5');
        itemTitle.className = 'mb-0';
        itemTitle.textContent = `Item ${itemIndex + 1}`;

        const removeButton = document.createElement('button');
        removeButton.className = 'btn btn-sm btn-danger';
        removeButton.innerHTML = '<i class="fas fa-times"></i> Remove';
        removeButton.addEventListener('click', () => {
            // Remove item from array
            arrayData.splice(itemIndex, 1);

            // Remove item from DOM
            itemDiv.remove();

            // Update indices for remaining items
            const items = itemsContainer.querySelectorAll('.array-item');
            items.forEach((item, idx) => {
                item.dataset.index = idx;
                item.querySelector('h5').textContent = `Item ${idx + 1}`;
            });

            // Update preview
            updatePreview();
        });

        itemHeader.appendChild(itemTitle);
        itemHeader.appendChild(removeButton);
        itemDiv.appendChild(itemHeader);

        // Create item body
        const itemBody = document.createElement('div');
        itemBody.className = 'card-body';
        itemDiv.appendChild(itemBody);

        // Create new item data
        let newItemData;

        // Handle different item types
        if (arraySchema.items) {
            // Handle $ref in items schema
            let itemsSchema = arraySchema.items;
            if (itemsSchema.$ref && itemsSchema.$ref.startsWith('#/$defs/')) {
                const refPath = itemsSchema.$ref.replace('#/$defs/', '');
                if (schema.$defs && schema.$defs[refPath]) {
                    itemsSchema = schema.$defs[refPath];
                    console.log(`Resolved items schema reference for ${path}: ${refPath}`);
                } else {
                    console.error(`Could not resolve items schema reference for ${path}`, itemsSchema.$ref);
                }
            }

            if (itemsSchema.type === 'object') {
                newItemData = createEmptyObject(itemsSchema);
                arrayData.push(newItemData);
                generateObjectFields(itemsSchema, newItemData, itemBody, itemPath);
            } else if (itemsSchema.type === 'array') {
                newItemData = [];
                arrayData.push(newItemData);
                generateArrayFields(itemsSchema, newItemData, itemBody, itemPath);
            } else {
                // For primitive types
                newItemData = itemsSchema.default !== undefined ? itemsSchema.default : null;
                arrayData.push(newItemData);
                generatePrimitiveField(itemsSchema, newItemData, itemBody, itemPath, 'value');
            }
        } else {
            // Default to empty object if items schema not specified
            newItemData = {};
            arrayData.push(newItemData);
        }

        // Add item to container
        itemsContainer.appendChild(itemDiv);

        // Update preview
        updatePreview();
    };

    // Add button to add new items
    const addButton = document.createElement('button');
    addButton.className = 'btn btn-primary mt-2';
    addButton.innerHTML = '<i class="fas fa-plus"></i> Add Item';
    addButton.addEventListener('click', addItem);
    arrayContainer.appendChild(addButton);

    // Generate fields for existing items
    if (Array.isArray(arrayData)) {
        arrayData.forEach((itemData, itemIndex) => {
            const itemPath = `${path}[${itemIndex}]`;

            // Create item container
            const itemDiv = document.createElement('div');
            itemDiv.className = 'array-item card mb-3';
            itemDiv.dataset.index = itemIndex;

            // Create item header
            const itemHeader = document.createElement('div');
            itemHeader.className = 'card-header d-flex justify-content-between align-items-center';

            const itemTitle = document.createElement('h5');
            itemTitle.className = 'mb-0';
            itemTitle.textContent = `Item ${itemIndex + 1}`;

            const removeButton = document.createElement('button');
            removeButton.className = 'btn btn-sm btn-danger';
            removeButton.innerHTML = '<i class="fas fa-times"></i> Remove';
            removeButton.addEventListener('click', () => {
                // Remove item from array
                arrayData.splice(itemIndex, 1);

                // Remove item from DOM
                itemDiv.remove();

                // Update indices for remaining items
                const items = itemsContainer.querySelectorAll('.array-item');
                items.forEach((item, idx) => {
                    item.dataset.index = idx;
                    item.querySelector('h5').textContent = `Item ${idx + 1}`;
                });

                // Update preview
                updatePreview();
            });

            itemHeader.appendChild(itemTitle);
            itemHeader.appendChild(removeButton);
            itemDiv.appendChild(itemHeader);

            // Create item body
            const itemBody = document.createElement('div');
            itemBody.className = 'card-body';
            itemDiv.appendChild(itemBody);

            // Handle different item types
            if (arraySchema.items) {
                // Handle $ref in items schema
                let itemsSchema = arraySchema.items;
                if (itemsSchema.$ref && itemsSchema.$ref.startsWith('#/$defs/')) {
                    const refPath = itemsSchema.$ref.replace('#/$defs/', '');
                    if (schema.$defs && schema.$defs[refPath]) {
                        itemsSchema = schema.$defs[refPath];
                    }
                }

                if (itemsSchema.type === 'object') {
                    generateObjectFields(itemsSchema, itemData, itemBody, itemPath);
                } else if (itemsSchema.type === 'array') {
                    generateArrayFields(itemsSchema, itemData, itemBody, itemPath);
                } else {
                    // For primitive types
                    generatePrimitiveField(itemsSchema, itemData, itemBody, itemPath, 'value');
                }
            }

            // Add item to container
            itemsContainer.appendChild(itemDiv);
        });
    }

    // If array is empty, add one item by default
    if (!arrayData.length) {
        addItem();
    }
}

// Add the missing createFormField function
function createFormField(container, id, label, type, value, description, onChange, options = null, min = null, max = null) {
    console.log(`Creating form field: ${id}, type: ${type}`);

    const fieldDiv = document.createElement('div');
    fieldDiv.className = 'form-field';

    // Create label with description tooltip
    const labelDiv = document.createElement('div');
    labelDiv.className = 'field-label';

    const labelEl = document.createElement('label');
    labelEl.setAttribute('for', id);
    labelEl.textContent = label;

    labelDiv.appendChild(labelEl);

    if (description) {
        const icon = document.createElement('i');
        icon.className = 'fas fa-info-circle description-icon';
        icon.setAttribute('title', description);
        icon.setAttribute('data-bs-toggle', 'tooltip');
        icon.setAttribute('data-bs-placement', 'top');
        labelDiv.appendChild(icon);

        // Initialize tooltip if Bootstrap is available
        if (typeof bootstrap !== 'undefined' && bootstrap.Tooltip) {
            new bootstrap.Tooltip(icon);
        }

        console.log(`Added info icon for field: ${id}, description: ${description}`);
    } else {
        console.log(`No description available for field: ${id}`);
    }

    fieldDiv.appendChild(labelDiv);

    // Create input based on type
    let input;

    switch (type) {
        case 'text':
            input = document.createElement('input');
            input.type = 'text';
            input.className = 'form-control';
            input.id = id;
            input.value = value !== undefined && value !== null ? value : '';
            break;

        case 'number':
            input = document.createElement('input');
            input.type = 'number';
            input.className = 'form-control';
            input.id = id;
            input.value = value !== undefined && value !== null ? value : '';
            if (min !== null) input.min = min;
            if (max !== null) input.max = max;
            break;

        case 'textarea':
            input = document.createElement('textarea');
            input.className = 'form-control';
            input.id = id;
            input.rows = 3;
            input.value = value !== undefined && value !== null ? value : '';
            break;

        case 'select':
            input = document.createElement('select');
            input.className = 'form-select';
            input.id = id;

            if (options) {
                options.forEach(option => {
                    const optEl = document.createElement('option');
                    optEl.value = option.value;
                    optEl.textContent = option.text;
                    if (value === option.value) {
                        optEl.selected = true;
                    }
                    input.appendChild(optEl);
                });
            }
            break;

        case 'checkbox':
            input = document.createElement('div');
            input.className = 'form-check';

            const checkbox = document.createElement('input');
            checkbox.type = 'checkbox';
            checkbox.className = 'form-check-input';
            checkbox.id = id;
            checkbox.checked = value || false;

            const checkLabel = document.createElement('label');
            checkLabel.className = 'form-check-label';
            checkLabel.setAttribute('for', id);
            checkLabel.textContent = 'Enabled';

            input.appendChild(checkbox);
            input.appendChild(checkLabel);

            // Add change event
            checkbox.addEventListener('change', () => {
                onChange(checkbox.checked);
            });
            break;
    }

    // Add change event for non-checkbox inputs
    if (type !== 'checkbox') {
        input.addEventListener('change', () => {
            onChange(input.value);
        });
    }

    fieldDiv.appendChild(input);
    container.appendChild(fieldDiv);

    return input;
}

// Add the missing createValueWithDOIField function
function createValueWithDOIField(container, id, label, type, value, ref, description, onValueChange, onRefChange, options = null, min = null, max = null) {
    console.log(`Creating ValueWithDOI field: ${id}, type: ${type}`);

    const fieldDiv = document.createElement('div');
    fieldDiv.className = 'form-field value-with-doi';

    // Create label with description tooltip
    const labelDiv = document.createElement('div');
    labelDiv.className = 'field-label';

    const labelEl = document.createElement('label');
    labelEl.setAttribute('for', id);

    // Extract the actual field name from the path if it's in the format ValueWithDOI[FieldName]
    const bracketMatch = label.match(/ValueWithDOI\[(.*?)\]/);
    if (bracketMatch && bracketMatch[1]) {
        labelEl.textContent = bracketMatch[1];
    } else {
        labelEl.textContent = label;
    }

    labelDiv.appendChild(labelEl);

    if (description) {
        const icon = document.createElement('i');
        icon.className = 'fas fa-info-circle description-icon';
        icon.setAttribute('title', description);
        icon.setAttribute('data-bs-toggle', 'tooltip');
        icon.setAttribute('data-bs-placement', 'top');
        labelDiv.appendChild(icon);

        // Initialize tooltip if Bootstrap is available
        if (typeof bootstrap !== 'undefined' && bootstrap.Tooltip) {
            new bootstrap.Tooltip(icon);
        }

        console.log(`Added info icon for field: ${id}, description: ${description}`);
    } else {
        console.log(`No description available for field: ${id}`);
    }

    fieldDiv.appendChild(labelDiv);

    // Create value input
    const valueDiv = document.createElement('div');
    valueDiv.className = 'value-input';

    let input;

    switch (type) {
        case 'text':
            input = document.createElement('input');
            input.type = 'text';
            input.className = 'form-control';
            input.id = id;
            input.value = value !== undefined && value !== null ? value : '';
            break;

        case 'number':
            input = document.createElement('input');
            input.type = 'number';
            input.className = 'form-control';
            input.id = id;
            input.value = value !== undefined && value !== null ? value : '';
            if (min !== null) input.min = min;
            if (max !== null) input.max = max;
            break;

        case 'select':
            input = document.createElement('select');
            input.className = 'form-select';
            input.id = id;

            if (options) {
                options.forEach(option => {
                    const optEl = document.createElement('option');
                    optEl.value = option.value;
                    optEl.textContent = option.text;
                    if (value === option.value) {
                        optEl.selected = true;
                    }
                    input.appendChild(optEl);
                });
            }
            break;
    }

    // Add change event
    input.addEventListener('change', () => {
        let newValue = input.value;
        if (type === 'number' || type === 'select') {
            newValue = parseFloat(newValue);
            if (isNaN(newValue)) {
                newValue = null;
            }
        }
        onValueChange(newValue);
    });

    valueDiv.appendChild(input);
    fieldDiv.appendChild(valueDiv);

    // Create reference toggle
    const refToggle = document.createElement('div');
    refToggle.className = 'reference-toggle';
    refToggle.innerHTML = '<i class="fas fa-book"></i> ' + (ref ? 'Edit Reference' : 'Add Reference');
    refToggle.addEventListener('click', () => {
        refFields.classList.toggle('show');
        refToggle.innerHTML = '<i class="fas fa-book"></i> ' + (refFields.classList.contains('show') ? 'Hide Reference' : (ref ? 'Edit Reference' : 'Add Reference'));
    });

    fieldDiv.appendChild(refToggle);

    // Create reference fields
    const refFields = document.createElement('div');
    refFields.className = 'reference-fields';
    if (ref) {
        refFields.classList.add('show');
    }

    const refInput = document.createElement('input');
    refInput.type = 'text';
    refInput.className = 'form-control';
    refInput.id = `${id}-ref`;
    refInput.value = ref || '';
    refInput.placeholder = 'DOI or reference URL';

    refInput.addEventListener('change', () => {
        onRefChange(refInput.value || null);
    });

    refFields.appendChild(refInput);
    fieldDiv.appendChild(refFields);

    container.appendChild(fieldDiv);

    return { valueInput: input, refInput: refInput };
}

// Add the missing generateSiteFields function
function generateSiteFields() {
    console.log('Generating site fields...');

    const container = document.getElementById('site-form-container');
    if (!container) {
        console.error('Site form container not found');
        return;
    }

    // Clear the container
    container.innerHTML = '';

    // Check if schema is available
    if (!schema || !schema.properties || !schema.properties.site) {
        const errorDiv = document.createElement('div');
        errorDiv.className = 'alert alert-danger';
        errorDiv.textContent = 'Schema for site not found';
        container.appendChild(errorDiv);
        console.error('Schema for site not found', schema);
        return;
    }

    const siteSchema = schema.properties.site;

    // Initialize site array if it doesn't exist
    if (!configData.site) {
        configData.site = [createEmptyObject(siteSchema.items)];
    }

    // Generate array fields for site
    generateArrayFields(siteSchema, configData.site, container, 'site');

    console.log('Site fields generation completed');
}