import { createAndPatchRegistry } from "@exabyte-io/code.js/dist/context";

import {
    NWChemTotalEnergyContextProvider,
    QENEBContextProvider,
    QEPWXContextProvider,
    VASPContextProvider,
    VASPNEBContextProvider,
} from "../context";

export const adeProviders = {
    QEPWXInputDataManager: {
        providerCls: QEPWXContextProvider,
        config: { name: "input" },
    },
    QENEBInputDataManager: {
        providerCls: QENEBContextProvider,
        config: { name: "input" },
    },
    VASPInputDataManager: {
        providerCls: VASPContextProvider,
        config: { name: "input" },
    },
    VASPNEBInputDataManager: {
        providerCls: VASPNEBContextProvider,
        config: { name: "input" },
    },
    NWChemInputDataManager: {
        providerCls: NWChemTotalEnergyContextProvider,
        config: { name: "input" },
    },
};

export const ContextProviderRegistry = createAndPatchRegistry(adeProviders);
