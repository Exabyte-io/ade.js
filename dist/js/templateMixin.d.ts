import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
import type { NamedInMemoryEntity } from "@mat3ra/code/dist/js/entity/mixins/NamedEntityMixin";
import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import type { AnyObject } from "@mat3ra/esse/dist/js/esse/types";
import type { ContextProvider, ContextProviderConfig, ContextProviderName } from "./context/ContextProvider";
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
    render: (externalContext: Record<string, unknown>) => void;
    getRenderedJSON: (context: Record<string, unknown>) => AnyObject;
    _cleanRenderingContext: (object: Record<string, unknown>) => Record<string, unknown>;
    getDataFromProvidersForRenderingContext: (context: Record<string, unknown>) => Record<string, unknown>;
    setContent: (text: string) => void;
    setRendered: (text: string) => void;
    getContextProvidersAsClassInstances: (providerContext: Record<string, unknown>) => ContextProvider[];
    getDataFromProvidersForPersistentContext: (providerContext: Record<string, unknown>) => Record<string, unknown>;
    getRenderingContext: (externalContext: Record<string, unknown>) => Record<string, unknown>;
};
export declare function templateMixin(item: TemplateBase): TemplateMixin & InMemoryEntity & {
    setName(name: string): void;
    name: string;
};
export type ContextProviderConfigMapEntry = {
    providerCls: typeof ContextProvider;
    config: ContextProviderConfig;
};
export type ContextProviderConfigMap = Record<ContextProviderName, ContextProviderConfigMapEntry>;
export type TemplateStaticMixin = {
    fromFlavor: (appName: string, execName: string, inputName: string) => TemplateMixin & TemplateBase;
    contextProviderRegistry: ContextProviderRegistryContainer | null;
    setContextProvidersConfig: (classConfigMap: ContextProviderConfigMap) => void;
};
export declare function templateStaticMixin(item: Constructor<TemplateBase & TemplateMixin>): TemplateStaticMixin & Constructor<InMemoryEntity & {
    setName(name: string): void;
    name: string;
} & TemplateMixin>;
