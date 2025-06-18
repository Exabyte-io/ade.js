"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const context_1 = require("@mat3ra/code/dist/js/context");
class ExecutableContextProvider extends context_1.ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
exports.default = ExecutableContextProvider;
