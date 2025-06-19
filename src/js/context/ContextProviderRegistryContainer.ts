import type {
    ContextProvider,
    ContextProviderConfig,
    ContextProviderInstance,
} from "./ContextProvider";

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
        this.providers = this.providers.filter((p) => p.name === providerCls.name);
    }

    removeProviderByName(name: string) {
        this.providers = this.providers.filter((p) => p.name === name);
    }
}

const registryContainer = new ContextProviderRegistryContainer();

/** Extends an existing context provider registry container and patches static class variables if applicable.
 * @example
 * const classConfigMap = {
 *     PlanewaveCutoffDataManager: {
 *         providerCls: PlanewaveCutoffsContextProvider,
 *         config: _makeImportant({ name: "cutoffs", entityName: "subworkflow" })
 *     },
 * };
 */
export const extendAndPatchRegistry = (
    classConfigMap: Record<
        string,
        { providerCls: typeof ContextProvider; config: ContextProviderConfig }
    >,
) => {
    Object.entries(classConfigMap).forEach(([name, { providerCls, config }]) => {
        registryContainer.addProvider({
            instance: providerCls.getConstructorConfig(config),
            name,
        });
    });
    return registryContainer;
};

/** Creates a new context provider registry container and patches static class variables if applicable.
 *
 * @param {Object} classConfigMap
 * @param {{Material: SpecificMockMaterial}} classesToPatch
 */
export const createAndPatchRegistry = (
    classConfigMap: Record<
        string,
        { providerCls: typeof ContextProvider; config: ContextProviderConfig }
    >,
) => {
    return extendAndPatchRegistry(classConfigMap);
};
