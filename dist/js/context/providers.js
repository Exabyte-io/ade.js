"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.ExecutableContextProvider = void 0;
const context_1 = require("@mat3ra/code/dist/js/context");
class ExecutableContextProvider extends context_1.ContextProvider {
    constructor(config) {
        super({
            ...config,
            domain: "executable",
        });
    }
}
exports.ExecutableContextProvider = ExecutableContextProvider;
