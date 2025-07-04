"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.allowedMonitors = exports.allowedResults = exports.allTemplates = exports.allApplications = exports.Template = exports.Flavor = exports.Executable = exports.Application = void 0;
const application_flavors_js_1 = require("@exabyte-io/application-flavors.js");
Object.defineProperty(exports, "allApplications", { enumerable: true, get: function () { return application_flavors_js_1.allApplications; } });
Object.defineProperty(exports, "allowedMonitors", { enumerable: true, get: function () { return application_flavors_js_1.allowedMonitors; } });
Object.defineProperty(exports, "allowedResults", { enumerable: true, get: function () { return application_flavors_js_1.allowedResults; } });
Object.defineProperty(exports, "allTemplates", { enumerable: true, get: function () { return application_flavors_js_1.allTemplates; } });
const application_1 = __importDefault(require("./application"));
exports.Application = application_1.default;
const executable_1 = __importDefault(require("./executable"));
exports.Executable = executable_1.default;
const flavor_1 = __importDefault(require("./flavor"));
exports.Flavor = flavor_1.default;
const template_1 = __importDefault(require("./template"));
exports.Template = template_1.default;
