import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import { deepClone } from "@mat3ra/code/dist/js/utils";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import nunjucks from "nunjucks";

import ContextProvider, {
    type ContextProviderConfig,
    type ContextProviderName,
} from "./context/ContextProvider";
import ContextProviderRegistryContainer from "./context/ContextProviderRegistryContainer";

export type TemplateBase = InMemoryEntity & NamedInMemoryEntity;

export type TemplateMixin = {
    isManuallyChanged: boolean;
    content: string;
    rendered: string | undefined;
    applicationName: string | undefined;
    executableName: string | undefined;
    contextProviders: ContextProvider[];
    addContextProvider: (provider: ContextProvider) => void;
    removeContextProvider: (provider: ContextProvider) => void;
    render: (externalContext?: Record<string, unknown>) => void;
    getRenderedJSON: (context?: Record<string, unknown>) => AnyObject;
    _cleanRenderingContext: (object: Record<string, unknown>) => Record<string, unknown>;
    getDataFromProvidersForRenderingContext: (
        context?: Record<string, unknown>,
    ) => Record<string, unknown>;
    setContent: (text: string) => void;
    setRendered: (text: string) => void;
    getContextProvidersAsClassInstances: (
        providerContext?: Record<string, unknown>,
    ) => ContextProvider[];
    getDataFromProvidersForPersistentContext: (
        providerContext?: Record<string, unknown>,
    ) => Record<string, unknown>;
    getRenderingContext: (externalContext?: Record<string, unknown>) => Record<string, unknown>;
};

export function templateMixin(item: TemplateBase) {
    // @ts-ignore
    const properties: TemplateMixin & TemplateBase = {
        get isManuallyChanged() {
            return this.prop("isManuallyChanged", false);
        },

        get content() {
            return this.prop("content", "");
        },

        setContent(text: string) {
            return this.setProp("content", text);
        },

        get rendered() {
            return this.prop("rendered") || this.content;
        },

        setRendered(text: string) {
            return this.setProp("rendered", text);
        },

        get applicationName() {
            return this.prop<string>("applicationName");
        },

        get executableName() {
            return this.prop<string>("executableName");
        },

        get contextProviders() {
            return this.prop("contextProviders", []);
        },

        addContextProvider(provider: ContextProvider) {
            this.setProp("contextProviders", [...this.contextProviders, provider]);
        },

        removeContextProvider(provider: ContextProvider) {
            const contextProviders = this.contextProviders.filter((p) => {
                return p.name !== provider.name && p.domain !== provider.domain;
            });

            this.setProp("contextProviders", contextProviders);
        },

        render(externalContext?: Record<string, unknown>) {
            const renderingContext = this.getRenderingContext(externalContext);
            if (!this.isManuallyChanged) {
                try {
                    const template = nunjucks.compile(this.content);

                    // deepClone to pass JSON data without classes
                    const rendered = template.render(
                        this._cleanRenderingContext(renderingContext),
                    ) as string;

                    this.setRendered(this.isManuallyChanged ? rendered : rendered || this.content);
                } catch (e) {
                    console.log(`Template is not compiled: ${e}`);
                    console.log({
                        content: this.content,
                        _cleanRenderingContext: this._cleanRenderingContext(renderingContext),
                    });
                }
            }
        },

        getRenderedJSON(context?: Record<string, unknown>) {
            this.render(context);
            return this.toJSON();
        },

        // Remove "bulky" items and JSON stringify before passing it to rendering engine (eg. jinja) to compile.
        // This way the context should still be passed in full to contextProviders, but not to final text template.
        // eslint-disable-next-line class-methods-use-this
        _cleanRenderingContext(object: Record<string, unknown>) {
            // eslint-disable-next-line @typescript-eslint/no-unused-vars
            const { job, ...clone } = object;
            return deepClone(clone);
        },

        /*
         * @summary Initializes context provider class instances. `providerContext` is used to pass the data about any
         *          previously stored values. That is if data was previously saved in database, the context provider
         *          shall receive it on initialization through providerContext and prioritize this value over the default.
         */
        getContextProvidersAsClassInstances(providerContext?: Record<string, unknown>) {
            return this.contextProviders.map((p) => {
                const providerInstance = (
                    this.constructor as unknown as TemplateStaticMixin
                ).contextProviderRegistry?.findProviderInstanceByName(p.name);

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
        getDataFromProvidersForRenderingContext(providerContext?: Record<string, unknown>) {
            const result: AnyObject = {};
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
        getDataFromProvidersForPersistentContext(providerContext?: Record<string, unknown>) {
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
        getRenderingContext(externalContext?: Record<string, unknown>) {
            return {
                ...externalContext,
                ...this.getDataFromProvidersForRenderingContext(externalContext),
            };
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));

    return properties;
}

export type ContextProviderConfigMapEntry = {
    providerCls: typeof ContextProvider;
    config: ContextProviderConfig;
};

export type ContextProviderConfigMap = Partial<
    Record<ContextProviderName, ContextProviderConfigMapEntry>
>;

export type TemplateStaticMixin = {
    contextProviderRegistry: ContextProviderRegistryContainer | null;
    setContextProvidersConfig: (classConfigMap: ContextProviderConfigMap) => void;
};

export function templateStaticMixin(item: Constructor<TemplateBase & TemplateMixin>) {
    // @ts-ignore
    const properties: TemplateStaticMixin & Constructor<TemplateBase & TemplateMixin> = {
        contextProviderRegistry: null,

        setContextProvidersConfig(classConfigMap: ContextProviderConfigMap) {
            const contextProviderRegistry = new ContextProviderRegistryContainer();

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
