import { NWChemTotalEnergyContextProvider, QENEBContextProvider, QEPWXContextProvider, VASPContextProvider, VASPNEBContextProvider } from "../context";
export declare const adeProviders: {
    QEPWXInputDataManager: {
        providerCls: typeof QEPWXContextProvider;
        config: {
            name: string;
        };
    };
    QENEBInputDataManager: {
        providerCls: typeof QENEBContextProvider;
        config: {
            name: string;
        };
    };
    VASPInputDataManager: {
        providerCls: typeof VASPContextProvider;
        config: {
            name: string;
        };
    };
    VASPNEBInputDataManager: {
        providerCls: typeof VASPNEBContextProvider;
        config: {
            name: string;
        };
    };
    NWChemInputDataManager: {
        providerCls: typeof NWChemTotalEnergyContextProvider;
        config: {
            name: string;
        };
    };
};
export declare const ContextProviderRegistry: import("@mat3ra/code/dist/js/context").ContextProviderRegistryContainer;
