import type { ContextProviderInstance } from "./ContextProvider";
import ContextProvider from "./ContextProvider";
export default class ContextProviderRegistryContainer {
    _providers: {
        name: string;
        instance: ContextProviderInstance;
    }[];
    constructor();
    get providers(): {
        name: string;
        instance: ContextProviderInstance;
    }[];
    set providers(p: {
        name: string;
        instance: ContextProviderInstance;
    }[]);
    addProvider({ name, instance }: {
        name: string;
        instance: ContextProviderInstance;
    }): void;
    findProviderInstanceByName(name: string): ContextProviderInstance | undefined;
    removeProvider(providerCls: ContextProvider): void;
    removeProviderByName(name: string): void;
}
