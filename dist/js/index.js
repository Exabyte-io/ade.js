"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || (function () {
    var ownKeys = function(o) {
        ownKeys = Object.getOwnPropertyNames || function (o) {
            var ar = [];
            for (var k in o) if (Object.prototype.hasOwnProperty.call(o, k)) ar[ar.length] = k;
            return ar;
        };
        return ownKeys(o);
    };
    return function (mod) {
        if (mod && mod.__esModule) return mod;
        var result = {};
        if (mod != null) for (var k = ownKeys(mod), i = 0; i < k.length; i++) if (k[i] !== "default") __createBinding(result, mod, k[i]);
        __setModuleDefault(result, mod);
        return result;
    };
})();
Object.defineProperty(exports, "__esModule", { value: true });
exports.context = exports.ContextProviderRegistry = exports.ExecutableContextProvider = exports.allowedMonitors = exports.allowedResults = exports.allTemplates = exports.allApplications = exports.getApplication = exports.getAllApplications = exports.Template = exports.Flavor = exports.Executable = exports.Application = void 0;
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
Object.defineProperty(exports, "allApplications", { enumerable: true, get: function () { return application_flavors_js_1.allApplications; } });
Object.defineProperty(exports, "allowedMonitors", { enumerable: true, get: function () { return application_flavors_js_1.allowedMonitors; } });
Object.defineProperty(exports, "allowedResults", { enumerable: true, get: function () { return application_flavors_js_1.allowedResults; } });
Object.defineProperty(exports, "allTemplates", { enumerable: true, get: function () { return application_flavors_js_1.allTemplates; } });
const application_1 = require("./application");
Object.defineProperty(exports, "Application", { enumerable: true, get: function () { return application_1.Application; } });
const context = __importStar(require("./context"));
exports.context = context;
const providers_1 = require("./context/providers");
Object.defineProperty(exports, "ExecutableContextProvider", { enumerable: true, get: function () { return providers_1.ExecutableContextProvider; } });
const registry_1 = require("./context/registry");
Object.defineProperty(exports, "ContextProviderRegistry", { enumerable: true, get: function () { return registry_1.ContextProviderRegistry; } });
const executable_1 = require("./executable");
Object.defineProperty(exports, "Executable", { enumerable: true, get: function () { return executable_1.Executable; } });
const flavor_1 = require("./flavor");
Object.defineProperty(exports, "Flavor", { enumerable: true, get: function () { return flavor_1.Flavor; } });
const template_1 = require("./template");
Object.defineProperty(exports, "Template", { enumerable: true, get: function () { return template_1.Template; } });
const tree_1 = require("./tree");
Object.defineProperty(exports, "getAllApplications", { enumerable: true, get: function () { return tree_1.getAllApplications; } });
Object.defineProperty(exports, "getApplication", { enumerable: true, get: function () { return tree_1.getApplication; } });
