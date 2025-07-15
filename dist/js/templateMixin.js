"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.templateMixin = templateMixin;
exports.templateStaticMixin = templateStaticMixin;
const utils_1 = require("@mat3ra/code/dist/js/utils");
const nunjucks_1 = __importDefault(require("nunjucks"));
const ContextProviderRegistryContainer_1 = __importDefault(require("./context/ContextProviderRegistryContainer"));
function templateMixin(item) {
    // @ts-ignore
    const properties = {
        get isManuallyChanged() {
            return this.prop("isManuallyChanged", false);
        },
        get content() {
            return this.prop("content", "");
        },
        setContent(text) {
            return this.setProp("content", text);
        },
        get rendered() {
            return this.prop("rendered") || this.content;
        },
        setRendered(text) {
            return this.setProp("rendered", text);
        },
        get applicationName() {
            return this.prop("applicationName");
        },
        get executableName() {
            return this.prop("executableName");
        },
        get contextProviders() {
            return this.prop("contextProviders", []);
        },
        addContextProvider(provider) {
            this.setProp("contextProviders", [...this.contextProviders, provider]);
        },
        removeContextProvider(provider) {
            const contextProviders = this.contextProviders.filter((p) => {
                return p.name !== provider.name && p.domain !== provider.domain;
            });
            this.setProp("contextProviders", contextProviders);
        },
        render(externalContext) {
            const renderingContext = this.getRenderingContext(externalContext);
            if (!this.isManuallyChanged) {
                try {
                    const template = nunjucks_1.default.compile(this.content);
                    // deepClone to pass JSON data without classes
                    const rendered = template.render(this._cleanRenderingContext(renderingContext));
                    this.setRendered(this.isManuallyChanged ? rendered : rendered || this.content);
                }
                catch (e) {
                    console.log(`Template is not compiled: ${e}`);
                    console.log({
                        content: this.content,
                        _cleanRenderingContext: this._cleanRenderingContext(renderingContext),
                    });
                }
            }
        },
        getRenderedJSON(context) {
            this.render(context);
            return this.toJSON();
        },
        // Remove "bulky" items and JSON stringify before passing it to rendering engine (eg. jinja) to compile.
        // This way the context should still be passed in full to contextProviders, but not to final text template.
        // eslint-disable-next-line class-methods-use-this
        _cleanRenderingContext(object) {
            // eslint-disable-next-line @typescript-eslint/no-unused-vars
            const { job, ...clone } = object;
            return (0, utils_1.deepClone)(clone);
        },
        /*
         * @summary Initializes context provider class instances. `providerContext` is used to pass the data about any
         *          previously stored values. That is if data was previously saved in database, the context provider
         *          shall receive it on initialization through providerContext and prioritize this value over the default.
         */
        getContextProvidersAsClassInstances(providerContext) {
            return this.contextProviders.map((p) => {
                var _a;
                const providerInstance = (_a = this.constructor.contextProviderRegistry) === null || _a === void 0 ? void 0 : _a.findProviderInstanceByName(p.name);
                if (!providerInstance) {
                    throw new Error(`Provider ${p.name} not found`);
                }
                const clsInstance = new providerInstance.constructor({
                    ...providerInstance.config,
                    context: providerContext,
                });
                return clsInstance;
            });
        },
        /*
         * @summary Extracts the the data from all context providers for further use during render.
         */
        getDataFromProvidersForRenderingContext(providerContext) {
            const result = {};
            this.getContextProvidersAsClassInstances(providerContext).forEach((contextProvider) => {
                const context = contextProvider.yieldDataForRendering();
                Object.keys(context).forEach((key) => {
                    // merge context keys if they are objects otherwise override them.
                    result[key] =
                        result[key] !== null && typeof result[key] === "object"
                            ? // @ts-ignore
                                { ...result[key], ...context[key] }
                            : context[key];
                });
            });
            return result;
        },
        /*
         * @summary Extracts the the data from all context providers for further save in persistent context.
         */
        // TODO: optimize logic to prevent re-initializing the context provider classes again below, reuse above function
        getDataFromProvidersForPersistentContext(providerContext) {
            const result = {};
            this.getContextProvidersAsClassInstances(providerContext).forEach((contextProvider) => {
                // only save in the persistent context the data from providers that were edited (or able to be edited)
                Object.assign(result, contextProvider.isEdited ? contextProvider.yieldData() : {});
            });
            return result;
        },
        /*
         * @summary Combines rendering context (in order of preference):
         *        - context from templates initialized with external context
         *        - "external" context and
         */
        getRenderingContext(externalContext) {
            return {
                ...externalContext,
                ...this.getDataFromProvidersForRenderingContext(externalContext),
            };
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
function templateStaticMixin(item) {
    // @ts-ignore
    const properties = {
        contextProviderRegistry: null,
        setContextProvidersConfig(classConfigMap) {
            const contextProviderRegistry = new ContextProviderRegistryContainer_1.default();
            Object.entries(classConfigMap).forEach(([name, { providerCls, config }]) => {
                contextProviderRegistry.addProvider({
                    instance: providerCls.getConstructorConfig(config),
                    name,
                });
            });
            this.contextProviderRegistry = contextProviderRegistry;
        },
    };
    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
    return properties;
}
