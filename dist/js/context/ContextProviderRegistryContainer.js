"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.createAndPatchRegistry = exports.extendAndPatchRegistry = void 0;
class ContextProviderRegistryContainer {
    constructor() {
        this._providers = [];
    }
    get providers() {
        return this._providers;
    }
    set providers(p) {
        this._providers = p;
    }
    addProvider({ name, instance }) {
        this._providers.push({
            name,
            instance,
        });
    }
    findProviderInstanceByName(name) {
        const provider = this.providers.find((p) => p.name === name);
        return provider && provider.instance;
    }
    removeProvider(providerCls) {
        this.providers = this.providers.filter((p) => p.name === providerCls.name);
    }
    removeProviderByName(name) {
        this.providers = this.providers.filter((p) => p.name === name);
    }
}
exports.default = ContextProviderRegistryContainer;
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
const extendAndPatchRegistry = (classConfigMap) => {
    Object.entries(classConfigMap).forEach(([name, { providerCls, config }]) => {
        registryContainer.addProvider({
            instance: providerCls.getConstructorConfig(config),
            name,
        });
    });
    return registryContainer;
};
exports.extendAndPatchRegistry = extendAndPatchRegistry;
/** Creates a new context provider registry container and patches static class variables if applicable.
 *
 * @param {Object} classConfigMap
 * @param {{Material: SpecificMockMaterial}} classesToPatch
 */
const createAndPatchRegistry = (classConfigMap) => {
    return (0, exports.extendAndPatchRegistry)(classConfigMap);
};
exports.createAndPatchRegistry = createAndPatchRegistry;
