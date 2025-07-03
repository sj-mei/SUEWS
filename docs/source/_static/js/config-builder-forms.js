/**
 * SUEWS Configuration Builder - Forms Module
 * 
 * This module handles form generation and field rendering.
 */

window.configBuilder = window.configBuilder || {};
window.configBuilder.forms = {};

/**
 * Generate the main form
 */
window.configBuilder.forms.generateForm = function() {
    const schema = window.configBuilderState.schema;
    const configData = window.configBuilderState.configData;
    
    if (!schema || !configData) {
        console.error('Schema or config data not available');
        return;
    }
    
    // Clear existing form
    const formContainer = document.getElementById('general-form-container');
    if (formContainer) {
        formContainer.innerHTML = '';
    }
    
    // Generate Model Configuration fields
    window.configBuilder.forms.generateModelFields();
    
    // Generate Site Information fields if site exists
    if (schema.properties && schema.properties.site) {
        const siteContainer = document.getElementById('site-form-container');
        if (siteContainer) {
            siteContainer.innerHTML = '';
            
            // Ensure site data exists
            if (!configData.site) {
                configData.site = window.configBuilder.schema.createEmptyObject(schema.properties.site);
            }
            
            window.configBuilder.forms.generateObjectFields(
                schema.properties.site,
                configData.site,
                siteContainer,
                'site'
            );
        }
    }
    
    // Initial preview update
    window.configBuilder.preview.updatePreview();
};

/**
 * Generate model configuration fields
 */
window.configBuilder.forms.generateModelFields = function() {
    const schema = window.configBuilderState.schema;
    const configData = window.configBuilderState.configData;
    const modelContainer = document.getElementById('model-form-container');
    
    if (!modelContainer || !schema || !schema.properties || !schema.properties.model || !schema.properties.model.properties) {
        console.warn('Model container, schema or model properties not available');
        return;
    }
    
    modelContainer.innerHTML = '';
    
    // Ensure model data exists
    if (!configData.model) {
        configData.model = window.configBuilder.schema.createEmptyObject(schema.properties.model);
    }
    
    // Create tabs for model sections
    const tabContainer = document.createElement('div');
    tabContainer.className = 'model-tabs mb-4';
    
    const tabList = document.createElement('ul');
    tabList.className = 'nav nav-tabs';
    tabList.setAttribute('role', 'tablist');
    
    const tabContent = document.createElement('div');
    tabContent.className = 'tab-content mt-3';
    
    // Process each model property
    Object.keys(schema.properties.model.properties).forEach((propKey, index) => {
        const propSchema = schema.properties.model.properties[propKey];
        
        // Create tab item
        const tabItem = document.createElement('li');
        tabItem.className = 'nav-item';
        tabItem.setAttribute('role', 'presentation');
        
        const tabLink = document.createElement('a');
        tabLink.className = `nav-link ${index === 0 ? 'active' : ''}`;
        tabLink.id = `${propKey}-tab`;
        tabLink.href = `#${propKey}-content`;
        tabLink.setAttribute('data-bs-toggle', 'tab');
        tabLink.setAttribute('role', 'tab');
        tabLink.textContent = window.configBuilder.ui.formatFieldLabel(propKey, propSchema);
        
        tabItem.appendChild(tabLink);
        tabList.appendChild(tabItem);
        
        // Create tab pane
        const tabPane = document.createElement('div');
        tabPane.className = `tab-pane fade ${index === 0 ? 'show active' : ''}`;
        tabPane.id = `${propKey}-content`;
        tabPane.setAttribute('role', 'tabpanel');
        
        // Ensure model property data exists
        if (!configData.model[propKey]) {
            configData.model[propKey] = window.configBuilder.schema.createEmptyObject(propSchema);
        }
        
        // Handle $ref in schema
        let resolvedSchema = propSchema;
        if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
            const refPath = propSchema.$ref.replace('#/$defs/', '');
            if (schema.$defs && schema.$defs[refPath]) {
                resolvedSchema = schema.$defs[refPath];
            }
        }
        
        // Generate fields for this model section
        window.configBuilder.forms.generateObjectFields(resolvedSchema, configData.model[propKey], tabPane, `model.${propKey}`);
        
        tabContent.appendChild(tabPane);
    });
    
    tabContainer.appendChild(tabList);
    modelContainer.appendChild(tabContainer);
    modelContainer.appendChild(tabContent);
};

/**
 * Generate object fields
 */
window.configBuilder.forms.generateObjectFields = function(objSchema, objData, container, path) {
    console.log(`Generating object fields for ${path}...`);
    
    // Handle $ref in schema
    let resolvedSchema = objSchema;
    if (objSchema.$ref && objSchema.$ref.startsWith('#/$defs/')) {
        const refPath = objSchema.$ref.replace('#/$defs/', '');
        const schema = window.configBuilderState.schema;
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
    
    // Special handling for VerticalLayers to ensure arrays are synchronized
    if (path.includes('vertical_layers') && resolvedSchema.title === 'VerticalLayers') {
        window.configBuilder.arrays.ensureVerticalLayersSync(objData);
    }
    
    // Create sections for different property groups
    const sections = {};
    const defaultSection = document.createElement('div');
    defaultSection.className = 'form-section';
    container.appendChild(defaultSection);
    
    // Process each property in the schema
    Object.keys(resolvedSchema.properties).forEach(propKey => {
        const originalPropSchema = resolvedSchema.properties[propKey];
        let propSchema = originalPropSchema;
        
        // Handle $ref in property schema
        if (propSchema.$ref && propSchema.$ref.startsWith('#/$defs/')) {
            const refPath = propSchema.$ref.replace('#/$defs/', '');
            const schema = window.configBuilderState.schema;
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
                objData[propKey] = window.configBuilder.schema.createEmptyObject(propSchema);
                // Special handling for vertical_layers object creation
                if (propKey === 'vertical_layers' && objData[propKey]) {
                    window.configBuilder.arrays.ensureVerticalLayersSync(objData[propKey]);
                }
            } else if (propSchema.type === 'array') {
                objData[propKey] = [];
            } else {
                objData[propKey] = propSchema.default !== undefined ? propSchema.default : null;
            }
        }
        
        // Get the section for this property
        let section = defaultSection;
        if (originalPropSchema.section) {
            if (!sections[originalPropSchema.section]) {
                sections[originalPropSchema.section] = document.createElement('div');
                sections[originalPropSchema.section].className = 'form-section';
                const sectionTitle = document.createElement('h4');
                sectionTitle.textContent = originalPropSchema.section;
                sections[originalPropSchema.section].appendChild(sectionTitle);
                container.appendChild(sections[originalPropSchema.section]);
            }
            section = sections[originalPropSchema.section];
        }
        
        // Generate field based on property type
        const propPath = `${path}.${propKey}`;
        
        // Debug logging for roofs/walls
        if (propKey === 'roofs' || propKey === 'walls') {
            console.log(`Processing ${propKey} field:`, {
                propPath,
                propSchema,
                originalPropSchema,
                dataValue: objData[propKey],
                schemaType: propSchema.type
            });
        }
        
        // Check if this is a FlexibleRefValue (anyOf with RefValue and raw type)
        let isFlexibleRefValue = false;
        let hasArrayOption = false;
        if (originalPropSchema.anyOf && originalPropSchema.anyOf.length === 2) {
            // Check if one option is a RefValue and the other is a raw type
            const hasRefValue = originalPropSchema.anyOf.some(option => 
                option.$ref && option.$ref.includes('RefValue'));
            const hasRawType = originalPropSchema.anyOf.some(option => 
                option.type && ['number', 'integer', 'string', 'boolean'].includes(option.type));
            const hasEnum = originalPropSchema.anyOf.some(option => 
                option.$ref && window.configBuilderState.schema.$defs && 
                window.configBuilderState.schema.$defs[option.$ref.replace('#/$defs/', '')] &&
                window.configBuilderState.schema.$defs[option.$ref.replace('#/$defs/', '')].enum);
            const arrayOption = originalPropSchema.anyOf.find(option => 
                option.type === 'array');
            
            if (hasRefValue && hasRawType) {
                isFlexibleRefValue = true;
            } else if (hasRefValue && arrayOption) {
                // This is a FlexibleRefValue with array option - handle as array
                hasArrayOption = true;
            }
        }
        
        // Special handling for ValueWithDOI type or FlexibleRefValue
        if ((propSchema.properties &&
            propSchema.properties.value &&
            propSchema.properties.ref) || isFlexibleRefValue) {
            window.configBuilder.forms.generateValueWithDOIField(originalPropSchema, objData[propKey], section, propPath, propKey);
        } else if (hasArrayOption || (originalPropSchema.anyOf && originalPropSchema.anyOf.some(opt => opt.type === 'array'))) {
            // Handle anyOf with array option as a regular array
            const arrayOption = originalPropSchema.anyOf.find(option => option.type === 'array');
            if (arrayOption) {
                console.log(`Handling ${propKey} as array from anyOf`);
                // Initialize array if null
                if (objData[propKey] === null || objData[propKey] === undefined) {
                    objData[propKey] = arrayOption.default || [];
                }
                
                // Check if this is an inline array that needs a label
                const isVerticalLayerInlineArray = path.includes('vertical_layers') && 
                    (propKey === 'height' || propKey === 'veg_frac' || propKey === 'veg_scale' || 
                     propKey === 'building_frac' || propKey === 'building_scale');
                
                // Note: roofs and walls are arrays of objects, not inline primitive arrays
                // They should use the standard collapsible card format
                const needsInlineLabel = isVerticalLayerInlineArray;
                
                if (needsInlineLabel) {
                    // For inline arrays, create a simple container with label
                    const fieldContainer = document.createElement('div');
                    fieldContainer.className = 'form-field mb-3';
                    
                    // Add field label
                    const label = document.createElement('label');
                    label.className = 'form-label';
                    label.textContent = window.configBuilder.ui.formatFieldLabel(propKey, originalPropSchema);
                    fieldContainer.appendChild(label);
                    
                    // Add unit if available
                    if (originalPropSchema.unit) {
                        window.configBuilder.ui.renderUnit(originalPropSchema.unit, label);
                    }
                    
                    // Generate the array fields
                    window.configBuilder.forms.generateArrayFields(arrayOption, objData[propKey], fieldContainer, propPath);
                    
                    section.appendChild(fieldContainer);
                } else {
                    // For regular arrays (like roofs/walls), create a collapsible card
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
                    collapseButton.textContent = window.configBuilder.ui.formatFieldLabel(propKey, originalPropSchema);

                    cardTitle.appendChild(collapseButton);
                    cardHeader.appendChild(cardTitle);
                    card.appendChild(cardHeader);

                    const collapseDiv = document.createElement('div');
                    collapseDiv.id = `collapse-${propPath.replace(/\./g, '-')}`;
                    collapseDiv.className = 'collapse';

                    const cardBody = document.createElement('div');
                    cardBody.className = 'card-body';

                    // Generate fields for array
                    window.configBuilder.forms.generateArrayFields(arrayOption, objData[propKey], cardBody, propPath);

                    collapseDiv.appendChild(cardBody);
                    card.appendChild(collapseDiv);
                    section.appendChild(card);
                }
            } else {
                // Fallback to primitive field
                window.configBuilder.forms.generatePrimitiveField(originalPropSchema, objData[propKey], section, propPath, propKey);
            }
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
            collapseButton.textContent = window.configBuilder.ui.formatFieldLabel(propKey, originalPropSchema);
            
            cardTitle.appendChild(collapseButton);
            cardHeader.appendChild(cardTitle);
            card.appendChild(cardHeader);
            
            const collapseDiv = document.createElement('div');
            collapseDiv.id = `collapse-${propPath.replace(/\./g, '-')}`;
            collapseDiv.className = 'collapse';
            
            const cardBody = document.createElement('div');
            cardBody.className = 'card-body';
            
            // Generate fields for nested object
            window.configBuilder.forms.generateObjectFields(propSchema, objData[propKey], cardBody, propPath);
            
            collapseDiv.appendChild(cardBody);
            card.appendChild(collapseDiv);
            section.appendChild(card);
        } else if (propSchema.type === 'array') {
            // Ensure the array exists in data so that fields render even when initially missing (e.g., roofs / walls)
            if (objData[propKey] === undefined || objData[propKey] === null) {
                objData[propKey] = Array.isArray(propSchema.default) ? JSON.parse(JSON.stringify(propSchema.default)) : [];
            }
            
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
            collapseButton.textContent = window.configBuilder.ui.formatFieldLabel(propKey, originalPropSchema);
            
            cardTitle.appendChild(collapseButton);
            cardHeader.appendChild(cardTitle);
            card.appendChild(cardHeader);
            
            const collapseDiv = document.createElement('div');
            collapseDiv.id = `collapse-${propPath.replace(/\./g, '-')}`;
            collapseDiv.className = 'collapse';
            
            const cardBody = document.createElement('div');
            cardBody.className = 'card-body';
            
            // Generate fields for array
            window.configBuilder.forms.generateArrayFields(propSchema, objData[propKey], cardBody, propPath);
            
            collapseDiv.appendChild(cardBody);
            card.appendChild(collapseDiv);
            section.appendChild(card);
        } else {
            // Generate field for primitive type
            console.log(`Generating primitive field for ${propKey}, type: ${propSchema.type}, value:`, objData[propKey]);
            window.configBuilder.forms.generatePrimitiveField(originalPropSchema, objData[propKey], section, propPath, propKey);
        }
    });
};

/**
 * Generate field for ValueWithDOI type
 */
window.configBuilder.forms.generateValueWithDOIField = function(propSchema, propData, container, path, propKey) {
    // Handle $ref in schema
    let resolvedSchema = propSchema;
    const schema = window.configBuilderState.schema;
    
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
    
    // Handle FlexibleRefValue (anyOf case)
    let valueSchema = null;
    let isFlexibleRefValue = false;
    
    if (propSchema.anyOf && propSchema.anyOf.length === 2) {
        // This is a FlexibleRefValue - find the RefValue option
        const refValueOption = propSchema.anyOf.find(option => 
            option.$ref && option.$ref.includes('RefValue'));
        const enumOption = propSchema.anyOf.find(option => 
            option.$ref && !option.$ref.includes('RefValue'));
        
        if (refValueOption && enumOption) {
            isFlexibleRefValue = true;
            // Resolve the RefValue schema
            const refValuePath = refValueOption.$ref.replace('#/$defs/', '');
            if (schema.$defs && schema.$defs[refValuePath]) {
                resolvedSchema = schema.$defs[refValuePath];
            }
            
            // Get the enum schema for the value
            const enumPath = enumOption.$ref.replace('#/$defs/', '');
            if (schema.$defs && schema.$defs[enumPath]) {
                valueSchema = schema.$defs[enumPath];
            }
        }
    }
    
    // Initialize data if needed
    // For FlexibleRefValue, propData might be either a raw value or a RefValue object
    if (isFlexibleRefValue && propData !== null && typeof propData !== 'object') {
        // Convert raw value to RefValue object
        propData = { value: propData, ref: null };
        window.configBuilder.forms.setNestedProperty(window.configBuilderState.configData, path, propData);
    } else if (!propData || typeof propData !== 'object') {
        propData = { value: null, ref: null };
        window.configBuilder.forms.setNestedProperty(window.configBuilderState.configData, path, propData);
    }
    
    // If not already set by FlexibleRefValue handling, get value schema normally
    if (!valueSchema) {
        valueSchema = resolvedSchema.properties && resolvedSchema.properties.value 
            ? resolvedSchema.properties.value 
            : { type: 'string' };
        
        // Handle $ref in value schema
        if (valueSchema.$ref && valueSchema.$ref.startsWith('#/$defs/')) {
            const refPath = valueSchema.$ref.replace('#/$defs/', '');
            if (schema.$defs && schema.$defs[refPath]) {
                valueSchema = schema.$defs[refPath];
            }
        }
    }
    
    // Get field type for value
    const fieldType = window.configBuilder.forms.getInputType(valueSchema);
    
    // Get field options for enum
    const options = window.configBuilder.schema.getEnumOptions(valueSchema);
    
    // Get min/max for number fields
    const min = valueSchema.minimum !== undefined ? valueSchema.minimum : null;
    const max = valueSchema.maximum !== undefined ? valueSchema.maximum : null;
    
    // Use the helper function to format the field label with units
    const displayLabel = window.configBuilder.ui.formatFieldLabel(propKey, propSchema);
    
    // Get the description from the property schema first, then from the resolved schema
    const description = propSchema.description || resolvedSchema.description || null;
    
    console.log(`ValueWithDOI field: ${path}, label: ${displayLabel}, unit: ${propSchema.unit}, display_name: ${propSchema.display_name}`);
    
    // Create form field with DOI
    window.configBuilder.forms.createValueWithDOIField(
        container,
        path.replace(/\./g, '-'),
        displayLabel,
        fieldType,
        propData.value,
        propData.ref,
        description,
        (value) => {
            // Convert value to appropriate type
            const convertedValue = window.configBuilder.forms.convertValueToType(value, valueSchema.type);
            
            // Update value
            propData.value = convertedValue;
            
            // Special handling for nlayer field in vertical_layers
            console.log(`Field changed at path: ${path}, value: ${convertedValue}`);
            if (path.includes('vertical_layers') && path.includes('nlayer')) {
                console.log('Detected nlayer change, synchronizing arrays...');
                window.configBuilder.arrays.synchronizeVerticalLayerArrays(path, convertedValue);
                // Also synchronize initial state arrays that depend on nlayer
                window.configBuilder.arrays.synchronizeInitialStateArrays(convertedValue);
            }
            
            // Update preview
            window.configBuilder.preview.updatePreview();
        },
        (ref) => {
            // Update ref
            propData.ref = ref;
            
            // Update preview
            window.configBuilder.preview.updatePreview();
        },
        options,
        min,
        max,
        propSchema.unit  // Pass the unit
    );
};

/**
 * Generate primitive field
 */
window.configBuilder.forms.generatePrimitiveField = function(propSchema, propData, container, path, propKey) {
    const fieldType = window.configBuilder.forms.getInputType(propSchema);
    const options = window.configBuilder.schema.getEnumOptions(propSchema);
    
    // Get min/max for number fields
    const min = propSchema.minimum !== undefined ? propSchema.minimum : null;
    const max = propSchema.maximum !== undefined ? propSchema.maximum : null;
    
    // Use the helper function to format the field label
    const displayLabel = window.configBuilder.ui.formatFieldLabel(propKey, propSchema);
    
    // Create form field
    window.configBuilder.forms.createFormField(
        container,
        path.replace(/\./g, '-'),
        displayLabel,
        fieldType,
        propData,
        propSchema.description,
        (value) => {
            // Convert value to appropriate type
            const convertedValue = window.configBuilder.forms.convertValueToType(value, propSchema.type);
            
            // Update data
            window.configBuilder.forms.setNestedProperty(window.configBuilderState.configData, path, convertedValue);
            
            // Update preview
            window.configBuilder.preview.updatePreview();
        },
        options,
        min,
        max,
        propSchema.unit  // Pass the unit
    );
};

/**
 * Get input type from schema
 */
window.configBuilder.forms.getInputType = function(propSchema) {
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
};

/**
 * Convert value to appropriate type
 */
window.configBuilder.forms.convertValueToType = function(value, type) {
    switch (type) {
        case 'number':
            return parseFloat(value);
        case 'integer':
            return parseInt(value);
        case 'boolean':
            return value === true || value === 'true';
        default:
            return value;
    }
};

/**
 * Set nested property in object
 */
window.configBuilder.forms.setNestedProperty = function(obj, path, value) {
    const parts = path.split('.');
    let current = obj;
    
    for (let i = 0; i < parts.length - 1; i++) {
        const part = parts[i];
        if (!current[part]) {
            current[part] = {};
        }
        current = current[part];
    }
    
    current[parts[parts.length - 1]] = value;
};

/**
 * Generate array fields
 */
window.configBuilder.forms.generateArrayFields = function(arraySchema, arrayData, container, path) {
    console.log(`Generating array fields for ${path}...`);
    
    const schema = window.configBuilderState.schema;
    
    // Check if this is a vertical layer primitive array
    const isVerticalLayerArray = path.includes('vertical_layers') && 
        (path.includes('height') || path.includes('veg_frac') || path.includes('veg_scale') || 
         path.includes('building_frac') || path.includes('building_scale'));
    
    // Only vertical layer arrays should be inline
    const isInlineArray = isVerticalLayerArray;
    
    const isPrimitiveArray = arraySchema.items && 
        (arraySchema.items.type === 'number' || arraySchema.items.type === 'integer' || 
         arraySchema.items.type === 'string' || arraySchema.items.type === 'boolean');
    
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
    
    // For inline primitive arrays, use a special inline display
    if (isInlineArray && isPrimitiveArray) {
        window.configBuilder.forms.generateInlinePrimitiveArray(arraySchema, arrayData, arrayContainer, path);
        return;
    }
    
    // Create items container for regular arrays
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
        const isCopied = arrayData[itemIndex] && arrayData[itemIndex].__is_copied;
        itemTitle.textContent = `Item ${itemIndex + 1}${isCopied ? ' (copied)' : ''}`;
        
        // Create button container
        const buttonContainer = document.createElement('div');
        buttonContainer.className = 'd-flex gap-2';
        
        // Create copy button
        const copyButton = document.createElement('button');
        copyButton.className = 'btn btn-sm btn-secondary';
        copyButton.innerHTML = '<i class="fas fa-copy"></i> Copy';
        
        // Check if this is a vertical layer array or initial state array
        const isVerticalLayerArrayForCopy = path.includes('vertical_layers') && 
            (path.includes('height') || path.includes('veg_frac') || path.includes('veg_scale') || 
             path.includes('building_frac') || path.includes('building_scale') || 
             path.includes('roofs') || path.includes('walls'));
        
        const isInitialStateArrayForCopy = path.includes('initial_state') && 
            (path.includes('roofs') || path.includes('walls'));
        
        if (isVerticalLayerArrayForCopy || isInitialStateArrayForCopy) {
            copyButton.style.display = 'none';
        }
        
        copyButton.addEventListener('click', () => {
            // Deep copy the item data
            const copiedData = JSON.parse(JSON.stringify(arrayData[itemIndex]));
            
            // Mark as copied
            copiedData.__is_copied = true;
            
            // Add the copied item to the array
            arrayData.push(copiedData);
            
            // Add the new item to the UI
            addItem();
            
            // Update preview
            window.configBuilder.preview.updatePreview();
        });
        
        const removeButton = document.createElement('button');
        removeButton.className = 'btn btn-sm btn-danger';
        removeButton.innerHTML = '<i class="fas fa-times"></i> Remove';
        
        // Check if this is an nlayer-controlled array
        const isNlayerControlled = window.configBuilder.arrays.isNlayerControlledArray(path);
        
        if (isNlayerControlled) {
            removeButton.style.display = 'none';
        }
        
        removeButton.addEventListener('click', () => {
            if (!isNlayerControlled) {
                // Remove item from array
                arrayData.splice(itemIndex, 1);
                
                // Remove item from DOM
                itemDiv.remove();
                
                // Update indices for remaining items
                const items = itemsContainer.querySelectorAll('.array-item');
                items.forEach((item, idx) => {
                    item.dataset.index = idx;
                    const isCopied = arrayData[idx] && arrayData[idx].__is_copied;
                    item.querySelector('h5').textContent = `Item ${idx + 1}${isCopied ? ' (copied)' : ''}`;
                });
                
                // Update preview
                window.configBuilder.preview.updatePreview();
            }
        });
        
        buttonContainer.appendChild(copyButton);
        buttonContainer.appendChild(removeButton);
        
        itemHeader.appendChild(itemTitle);
        itemHeader.appendChild(buttonContainer);
        itemDiv.appendChild(itemHeader);
        
        // Create collapsible wrapper for item body
        const collapseDiv = document.createElement('div');
        collapseDiv.id = `collapse-${path.replace(/\./g, '-')}-${itemIndex}`;
        collapseDiv.className = 'collapse show'; // Show by default
        
        // Make the header clickable to toggle collapse
        itemHeader.style.cursor = 'pointer';
        itemHeader.setAttribute('data-bs-toggle', 'collapse');
        itemHeader.setAttribute('data-bs-target', `#${collapseDiv.id}`);
        
        // Create item body
        const itemBody = document.createElement('div');
        itemBody.className = 'card-body';
        
        collapseDiv.appendChild(itemBody);
        itemDiv.appendChild(collapseDiv);
        
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
                newItemData = window.configBuilder.schema.createEmptyObject(itemsSchema);
                arrayData.push(newItemData);
                window.configBuilder.forms.generateObjectFields(itemsSchema, newItemData, itemBody, itemPath);
            } else if (itemsSchema.type === 'array') {
                newItemData = [];
                arrayData.push(newItemData);
                window.configBuilder.forms.generateArrayFields(itemsSchema, newItemData, itemBody, itemPath);
            } else {
                // For primitive types
                newItemData = itemsSchema.default !== undefined ? itemsSchema.default : null;
                arrayData.push(newItemData);
                window.configBuilder.forms.generatePrimitiveField(itemsSchema, newItemData, itemBody, itemPath, 'value');
            }
        } else {
            // Default to empty object if items schema not specified
            newItemData = {};
            arrayData.push(newItemData);
        }
        
        // Add item to container
        itemsContainer.appendChild(itemDiv);
        
        // Update preview
        window.configBuilder.preview.updatePreview();
    };
    
    // Add button to add new items
    const addButton = document.createElement('button');
    addButton.className = 'btn btn-primary mt-2';
    addButton.innerHTML = '<i class="fas fa-plus"></i> Add Item';
    
    // Check if this is an nlayer-controlled array
    const isNlayerControlled = window.configBuilder.arrays.isNlayerControlledArray(path);
    
    if (isNlayerControlled) {
        addButton.disabled = true;
        addButton.title = 'Array length is controlled by nlayer';
        addButton.innerHTML = '<i class="fas fa-lock"></i> Array length controlled by nlayer';
    }
    
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
            const isCopied = itemData && itemData.__is_copied;
            itemTitle.textContent = `Item ${itemIndex + 1}${isCopied ? ' (copied)' : ''}`;
            
            // Create button container
            const buttonContainer = document.createElement('div');
            buttonContainer.className = 'd-flex gap-2';
            
            // Create copy button
            const copyButton = document.createElement('button');
            copyButton.className = 'btn btn-sm btn-secondary';
            copyButton.innerHTML = '<i class="fas fa-copy"></i> Copy';
            
            // Check if this is an nlayer-controlled array
            if (isNlayerControlled) {
                copyButton.style.display = 'none';
            }
            
            copyButton.addEventListener('click', () => {
                // Deep copy the item data
                const copiedData = JSON.parse(JSON.stringify(itemData));
                
                // Mark as copied
                copiedData.__is_copied = true;
                
                // Add the copied item to the array
                arrayData.push(copiedData);
                
                // Regenerate the array to show the new item
                container.innerHTML = '';
                window.configBuilder.forms.generateArrayFields(arraySchema, arrayData, container, path);
                
                // Update preview
                window.configBuilder.preview.updatePreview();
            });
            
            const removeButton = document.createElement('button');
            removeButton.className = 'btn btn-sm btn-danger';
            removeButton.innerHTML = '<i class="fas fa-times"></i> Remove';
            
            if (isNlayerControlled) {
                removeButton.style.display = 'none';
            }
            
            removeButton.addEventListener('click', () => {
                if (!isNlayerControlled) {
                    // Remove item from array
                    arrayData.splice(itemIndex, 1);
                    
                    // Remove item from DOM
                    itemDiv.remove();
                    
                    // Update indices for remaining items
                    const items = itemsContainer.querySelectorAll('.array-item');
                    items.forEach((item, idx) => {
                        item.dataset.index = idx;
                        const isCopied = arrayData[idx] && arrayData[idx].__is_copied;
                        item.querySelector('h5').textContent = `Item ${idx + 1}${isCopied ? ' (copied)' : ''}`;
                    });
                    
                    // Update preview
                    window.configBuilder.preview.updatePreview();
                }
            });
            
            buttonContainer.appendChild(copyButton);
            buttonContainer.appendChild(removeButton);
            
            itemHeader.appendChild(itemTitle);
            itemHeader.appendChild(buttonContainer);
            itemDiv.appendChild(itemHeader);
            
            // Create collapsible wrapper for item body
            const collapseDiv = document.createElement('div');
            collapseDiv.id = `collapse-${path.replace(/\./g, '-')}-${itemIndex}`;
            collapseDiv.className = 'collapse show'; // Show by default
            
            // Make the header clickable to toggle collapse
            itemHeader.style.cursor = 'pointer';
            itemHeader.setAttribute('data-bs-toggle', 'collapse');
            itemHeader.setAttribute('data-bs-target', `#${collapseDiv.id}`);
            
            // Create item body
            const itemBody = document.createElement('div');
            itemBody.className = 'card-body';
            
            collapseDiv.appendChild(itemBody);
            itemDiv.appendChild(collapseDiv);
            
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
                    window.configBuilder.forms.generateObjectFields(itemsSchema, itemData, itemBody, itemPath);
                } else if (itemsSchema.type === 'array') {
                    window.configBuilder.forms.generateArrayFields(itemsSchema, itemData, itemBody, itemPath);
                } else {
                    // For primitive types
                    window.configBuilder.forms.generatePrimitiveField(itemsSchema, itemData, itemBody, itemPath, 'value');
                }
            }
            
            // Add item to container
            itemsContainer.appendChild(itemDiv);
        });
    }
    
    // Special handling for nlayer-controlled arrays
    if (isNlayerControlled) {
        // Hide the Add button
        addButton.style.display = 'none';
        
        // Also hide all remove buttons
        const allRemoveButtons = itemsContainer.querySelectorAll('.btn-danger');
        allRemoveButtons.forEach(btn => {
            btn.style.display = 'none';
        });
    }
    
    // If array is empty and not nlayer-controlled, add one item by default
    if (!arrayData.length && !isNlayerControlled) {
        addItem();
    }
};

/**
 * Generate inline primitive array
 */
window.configBuilder.forms.generateInlinePrimitiveArray = function(arraySchema, arrayData, container, path) {
    console.log(`Generating inline primitive array for ${path}`);
    
    // Create a container for the inline inputs
    const inlineContainer = document.createElement('div');
    inlineContainer.className = 'inline-array-container';
    inlineContainer.style.display = 'flex';
    inlineContainer.style.flexWrap = 'wrap';
    inlineContainer.style.gap = '10px';
    inlineContainer.style.marginTop = '10px';
    
    // Render each array element as an individual input
    arrayData.forEach((value, index) => {
        const itemDiv = document.createElement('div');
        itemDiv.className = 'inline-array-item';
        itemDiv.style.flex = '1';
        itemDiv.style.minWidth = '80px';
        itemDiv.style.maxWidth = '150px';
        
        // Create label
        const label = document.createElement('label');
        label.className = 'form-label small text-muted';
        label.textContent = `Layer ${index + 1}`;
        label.style.fontSize = '0.8rem';
        label.style.marginBottom = '2px';
        itemDiv.appendChild(label);
        
        // Create input
        const input = document.createElement('input');
        input.type = arraySchema.items.type === 'number' || arraySchema.items.type === 'integer' ? 'number' : 'text';
        input.className = 'form-control form-control-sm';
        input.value = value !== null && value !== undefined ? value : '';
        
        // Add min/max if specified
        if (arraySchema.items.minimum !== undefined) {
            input.min = arraySchema.items.minimum;
        }
        if (arraySchema.items.maximum !== undefined) {
            input.max = arraySchema.items.maximum;
        }
        
        // Add change handler
        input.addEventListener('change', (e) => {
            let newValue = e.target.value;
            
            // Convert to appropriate type
            if (arraySchema.items.type === 'number') {
                newValue = parseFloat(newValue);
                if (isNaN(newValue)) newValue = 0;
            } else if (arraySchema.items.type === 'integer') {
                newValue = parseInt(newValue);
                if (isNaN(newValue)) newValue = 0;
            }
            
            // Update array data
            arrayData[index] = newValue;
            
            // Update preview
            window.configBuilder.preview.updatePreview();
        });
        
        itemDiv.appendChild(input);
        inlineContainer.appendChild(itemDiv);
    });
    
    container.appendChild(inlineContainer);
    
    // Add note about nlayer control
    const noteDiv = document.createElement('div');
    noteDiv.className = 'text-muted small mt-2';
    noteDiv.innerHTML = '<i class="fas fa-info-circle"></i> Number of layers is controlled by nlayer parameter';
    container.appendChild(noteDiv);
};

/**
 * Create form field
 */
window.configBuilder.forms.createFormField = function(container, id, label, type, value, description, onChange, options = null, min = null, max = null, unit = null) {
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
    
    // Add unit display if available
    window.configBuilder.ui.renderUnit(unit, labelDiv);
    
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
            // Don't try to display arrays or objects as text
            if (Array.isArray(value) || (typeof value === 'object' && value !== null)) {
                console.warn(`Attempting to display ${typeof value} as text field:`, value);
                input.value = '';
                input.placeholder = 'Complex data - use array editor';
            } else {
                input.value = value !== undefined && value !== null ? value : '';
            }
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
            // Don't try to display arrays or objects as text
            if (Array.isArray(value) || (typeof value === 'object' && value !== null)) {
                console.warn(`Attempting to display ${typeof value} as textarea:`, value);
                input.value = '';
                input.placeholder = 'Complex data - use array editor';
            } else {
                input.value = value !== undefined && value !== null ? value : '';
            }
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
};

/**
 * Create value with DOI field
 */
window.configBuilder.forms.createValueWithDOIField = function(container, id, label, type, value, ref, description, onValueChange, onRefChange, options = null, min = null, max = null, unit = null) {
    console.log(`Creating ValueWithDOI field: ${id}, label: ${label}, type: ${type}`);
    
    const fieldDiv = document.createElement('div');
    fieldDiv.className = 'form-field value-with-doi';
    
    // Create label with description tooltip
    const labelDiv = document.createElement('div');
    labelDiv.className = 'field-label';
    
    const labelEl = document.createElement('label');
    labelEl.setAttribute('for', id);
    labelEl.textContent = label;
    
    labelDiv.appendChild(labelEl);
    
    // Add unit display if available
    window.configBuilder.ui.renderUnit(unit, labelDiv);
    
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
};