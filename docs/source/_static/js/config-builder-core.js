/**
 * SUEWS Configuration Builder - Core Module
 * 
 * This module contains the core functionality, global state, and initialization
 * for the SUEWS Configuration Builder.
 */

// Global state
window.configBuilderState = {
    schema: null,
    configData: {},
    isLoading: false,
    validationErrors: []
};

// Export functions to global scope for use by other modules
window.configBuilder = window.configBuilder || {};

/**
 * Initialize the configuration builder
 */
window.configBuilder.init = async function() {
    console.log('Initializing SUEWS Configuration Builder...');
    
    try {
        // Load the schema
        await window.configBuilder.schema.loadSchema();
        
        // Initialize empty configuration
        window.configBuilder.schema.initializeEmptyConfig();
        
        // Generate the form
        window.configBuilder.forms.generateForm();
        
        // Setup event listeners
        window.configBuilder.ui.setupEventListeners();
        
        // Initial preview update
        window.configBuilder.preview.updatePreview();
        
        console.log('Configuration builder initialized successfully');
    } catch (error) {
        console.error('Failed to initialize configuration builder:', error);
        window.configBuilder.ui.displayError('Failed to initialize: ' + error.message);
    }
};

/**
 * Get the current configuration data
 */
window.configBuilder.getConfigData = function() {
    return window.configBuilderState.configData;
};

/**
 * Set the configuration data
 */
window.configBuilder.setConfigData = function(data) {
    window.configBuilderState.configData = data;
};

/**
 * Get the schema
 */
window.configBuilder.getSchema = function() {
    return window.configBuilderState.schema;
};

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
    // Check if all required modules are loaded
    const requiredModules = ['schema', 'forms', 'arrays', 'ui', 'preview'];
    const missingModules = requiredModules.filter(module => !window.configBuilder[module]);
    
    if (missingModules.length > 0) {
        console.error('Missing required modules:', missingModules);
        return;
    }
    
    // Initialize the builder
    window.configBuilder.init();
});