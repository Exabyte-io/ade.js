import { ContextProviderRegistryContainer } from "@exabyte-io/code.js/context";
import _ from "underscore";

import {
    NWChemTotalEnergyContextProvider,
    QENEBContextProvider,
    QEPWXContextProvider,
    VASPContextProvider,
    VASPNEBContextProvider,
} from "../context";

const adeContextProviders = {
    QEPWXInputDataManager: QEPWXContextProvider.getConstructorConfig({ name: "input" }),
    QENEBInputDataManager: QENEBContextProvider.getConstructorConfig({ name: "input" }),
    VASPInputDataManager: VASPContextProvider.getConstructorConfig({ name: "input" }),
    VASPNEBInputDataManager: VASPNEBContextProvider.getConstructorConfig({ name: "input" }),
    NWChemInputDataManager: NWChemTotalEnergyContextProvider.getConstructorConfig({
        name: "input",
    }),
};

const clsInstance = new ContextProviderRegistryContainer();
_.map(adeContextProviders, (instance, name) =>
    clsInstance.addProvider({
        instance,
        name,
    }),
);

export const ContextProviderRegistry = clsInstance;
