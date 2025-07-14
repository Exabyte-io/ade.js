import type { ContextProviderInstance } from "./ContextProvider";
import ContextProvider from "./ContextProvider";

export default class ContextProviderRegistryContainer {
    _providers: {
        name: string;
        instance: ContextProviderInstance;
    }[];

    constructor() {
        this._providers = [];
    }

    get providers() {
        return this._providers;
    }

    set providers(p) {
        this._providers = p;
    }

    addProvider({ name, instance }: { name: string; instance: ContextProviderInstance }) {
        this._providers.push({
            name,
            instance,
        });
    }

    findProviderInstanceByName(name: string) {
        const provider = this.providers.find((p) => p.name === name);
        return provider && provider.instance;
    }

    removeProvider(providerCls: ContextProvider) {
        this.providers = this.providers.filter((p) => p.name !== providerCls.name);
    }

    removeProviderByName(name: string) {
        this.providers = this.providers.filter((p) => p.name !== name);
    }
}
