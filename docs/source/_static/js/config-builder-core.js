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
        console.log('Step 1: Loading schema...');
        // Load the schema
        await window.configBuilder.schema.loadSchema();
        console.log('Schema loaded:', window.configBuilderState.schema ? 'Yes' : 'No');
        
        console.log('Step 2: Initializing empty configuration...');
        // Initialize empty configuration
        window.configBuilder.schema.initializeEmptyConfig();
        console.log('Config initialized:', window.configBuilderState.configData);
        
        console.log('Step 3: Generating form...');
        // Generate the form
        window.configBuilder.forms.generateForm();
        console.log('Form generated');
        
        console.log('Step 4: Setting up event listeners...');
        // Setup event listeners
        window.configBuilder.ui.setupEventListeners();
        console.log('Event listeners set up');
        
        console.log('Step 5: Updating preview...');
        // Initial preview update
        window.configBuilder.preview.updatePreview();
        console.log('Preview updated');
        
        console.log('Configuration builder initialized successfully');
    } catch (error) {
        console.error('Failed to initialize configuration builder:', error);
        console.error('Error stack:', error.stack);
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
// Use a more robust initialization approach
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeConfigBuilder);
} else {
    // DOM is already loaded, initialize immediately
    initializeConfigBuilder();
}

function initializeConfigBuilder() {
    console.log('initializeConfigBuilder called');
    console.log('window.configBuilder:', window.configBuilder);
    
    // Wait a bit for all modules to load
    setTimeout(function() {
        // Check if all required modules are loaded
        const requiredModules = ['schema', 'forms', 'arrays', 'ui', 'preview'];
        const loadedModules = requiredModules.filter(module => window.configBuilder && window.configBuilder[module]);
        const missingModules = requiredModules.filter(module => !window.configBuilder || !window.configBuilder[module]);
        
        console.log('Module check - Loaded:', loadedModules, 'Missing:', missingModules);
        
        if (missingModules.length > 0) {
            console.error('Missing required modules:', missingModules);
            // Try again after a delay
            setTimeout(initializeConfigBuilder, 100);
            return;
        }
        
        console.log('All modules loaded, calling init()');
        // Initialize the builder
        window.configBuilder.init();
    }, 100);
}