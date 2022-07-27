import _ from "underscore";
import { ContextProviderRegistryContainer } from "@exabyte-io/code.js/dist/context";

import {
    QEPWXContextProvider,
    QENEBContextProvider,
    VASPContextProvider,
    VASPNEBContextProvider,
    NWChemTotalEnergyContextProvider,
} from "../context";

const adeContextProviders = {
    QEPWXInputDataManager: QEPWXContextProvider.getConstructorConfig({ name: "input" }),
    QENEBInputDataManager: QENEBContextProvider.getConstructorConfig({ name: "input" }),
    VASPInputDataManager: VASPContextProvider.getConstructorConfig({ name: "input" }),
    VASPNEBInputDataManager: VASPNEBContextProvider.getConstructorConfig({ name: "input" }),
    NWChemInputDataManager: NWChemTotalEnergyContextProvider.getConstructorConfig({ name: "input" }),
}

const clsInstance = new ContextProviderRegistryContainer();
_.map(adeContextProviders, (instance, name) => clsInstance.addProvider({
    instance,
    name
}));

export const ContextProviderRegistry = clsInstance;
