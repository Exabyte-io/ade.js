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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.context = exports.allowedMonitors = exports.allowedResults = exports.allTemplates = exports.allApplications = exports.getApplication = exports.getAllApplications = exports.Template = exports.Flavor = exports.Executable = exports.Application = void 0;
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
Object.defineProperty(exports, "allApplications", { enumerable: true, get: function () { return application_flavors_js_1.allApplications; } });
Object.defineProperty(exports, "allowedMonitors", { enumerable: true, get: function () { return application_flavors_js_1.allowedMonitors; } });
Object.defineProperty(exports, "allowedResults", { enumerable: true, get: function () { return application_flavors_js_1.allowedResults; } });
Object.defineProperty(exports, "allTemplates", { enumerable: true, get: function () { return application_flavors_js_1.allTemplates; } });
const application_1 = __importDefault(require("./application"));
exports.Application = application_1.default;
const context = __importStar(require("./context"));
exports.context = context;
const executable_1 = __importDefault(require("./executable"));
exports.Executable = executable_1.default;
const flavor_1 = __importDefault(require("./flavor"));
exports.Flavor = flavor_1.default;
const template_1 = __importDefault(require("./template"));
exports.Template = template_1.default;
const tree_1 = require("./tree");
Object.defineProperty(exports, "getAllApplications", { enumerable: true, get: function () { return tree_1.getAllApplications; } });
Object.defineProperty(exports, "getApplication", { enumerable: true, get: function () { return tree_1.getApplication; } });
