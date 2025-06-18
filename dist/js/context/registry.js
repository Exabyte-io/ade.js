"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.ContextProviderRegistry = exports.adeProviders = void 0;
const context_1 = require("@mat3ra/code/dist/js/context");
const context_2 = require("../context");
exports.adeProviders = {
    QEPWXInputDataManager: {
        providerCls: context_2.QEPWXContextProvider,
        config: { name: "input" },
    },
    QENEBInputDataManager: {
        providerCls: context_2.QENEBContextProvider,
        config: { name: "input" },
    },
    VASPInputDataManager: {
        providerCls: context_2.VASPContextProvider,
        config: { name: "input" },
    },
    VASPNEBInputDataManager: {
        providerCls: context_2.VASPNEBContextProvider,
        config: { name: "input" },
    },
    NWChemInputDataManager: {
        providerCls: context_2.NWChemTotalEnergyContextProvider,
        config: { name: "input" },
    },
};
exports.ContextProviderRegistry = (0, context_1.createAndPatchRegistry)(exports.adeProviders);
