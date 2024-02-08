// @ts-expect-error application-flavors.js is not typed
import { allTemplates } from "@exabyte-io/application-flavors.js";
import {
    HashedEntityMixin,
    HashedInputArrayMixin,
    NamedInMemoryEntity,
} from "@exabyte-io/code.js/dist/entity";
import { deepClone } from "@exabyte-io/code.js/dist/utils";
// @ts-expect-error swig is not typed
import jinja from "swig";
import _ from "underscore";

import { ContextProviderRegistry } from "./context/registry";
import { Constructor, ContextProvider } from "@exabyte-io/code.js/dist/context";
import { TemplateData } from "./types";

const Base = HashedInputArrayMixin(HashedEntityMixin(NamedInMemoryEntity))
abstract class TemplateBaseEntity extends Base {};

export function TemplateMixin<
    T extends Constructor<TemplateBaseEntity> = Constructor<TemplateBaseEntity>,
>(superclass: T) {
    return class Template extends superclass {
        static providerRegistry = ContextProviderRegistry;

        get isManuallyChanged() {
            return this.prop<boolean>("isManuallyChanged", false);
        }

        get content() {
            return this.prop<string>("content");
        }

        setContent(text: string) {
            return this.setProp("content", text);
        }

        get rendered() {
            return this.prop<string>("rendered") || this.content;
        }

        setRendered(text: string) {
            return this.setProp("rendered", text);
        }

        get applicationName() {
            return this.prop<string>("applicationName");
        }

        get executableName() {
            return this.prop<string>("executableName");
        }

        get contextProviders() {
            return this.prop<ContextProvider[]>("contextProviders") || [];
        }

        addContextProvider(provider: ContextProvider) {
            this.setProp("contextProviders", this.contextProviders.push(provider));
        }

        removeContextProvider(provider: ContextProvider) {
            this.setProp(
                "contextProviders",
                this.contextProviders.filter(
                    (p) => p.name !== provider.name && p.domain !== provider.domain,
                ),
            );
        }

        render(externalContext: object) {
            const renderingContext = this.getRenderingContext(externalContext);
            let template, rendered;
            if (!this.isManuallyChanged) {
                try {
                    template = jinja.compile(this.content);
                    // deepClone to pass JSON data without classes
                    rendered = template && template(this._cleanRenderingContext(renderingContext));
                } catch (e) {
                    console.log(`Template is not compiled: ${e}`);
                }
                this.setRendered(this.isManuallyChanged ? rendered : rendered || this.content);
            }
        }

        getRenderedJSON(context: object) {
            this.render(context);
            return this.toJSON();
        }

        // Remove "bulky" items and JSON stringify before passing it to rendering engine (eg. jinja) to compile.
        // This way the context should still be passed in full to contextProviders, but not to final text template.
        // eslint-disable-next-line class-methods-use-this
        _cleanRenderingContext(object: { [key: string]: any }) {
            const { job, ...clone } = object;
            return deepClone(clone);
        }

        static fromFlavor(appName: string, execName: string, inputName: string): Template {
            const filtered = allTemplates.filter(
                (temp: TemplateData) =>
                    temp.applicationName === appName &&
                    temp.executableName === execName &&
                    temp.name === inputName,
            );
            if (filtered.length !== 1) {
                console.log(
                    `found ${filtered.length} templates for app=${appName} exec=${execName} name=${inputName} expected 1`,
                );
            }
            return new Template(filtered[0]);
        }

        /*
        * @summary Initializes context provider class instances. `providerContext` is used to pass the data about any
        *          previously stored values. That is if data was previously saved in database, the context provider
        *          shall receive it on initialization through providerContext and prioritize this value over the default.
        */
        getContextProvidersAsClassInstances(providerContext: object) {
            return this.contextProviders.map((p) => {
                const { constructor, config } =
                    Template.providerRegistry.findProviderInstanceByName(p.name);
                const clsInstance = new constructor({
                    ...config,
                    context: providerContext,
                });
                return clsInstance;
            });
        }

        /*
        * @summary Extracts the the data from all context providers for further use during render.
        */
        getDataFromProvidersForRenderingContext(providerContext: object): object {
            const result = {};
            this.getContextProvidersAsClassInstances(providerContext).forEach((contextProvider) => {
                const context = contextProvider.yieldDataForRendering();
                Object.keys(context).forEach((key) => {
                    // merge context keys if they are objects otherwise override them.
                    // @ts-ignore
                    result[key] = _.isObject(result[key])
                        // @ts-ignore
                        ? { ...result[key], ...context[key] }
                        : context[key];
                });
            });
            return result;
        }

        /*
        * @summary Extracts the the data from all context providers for further save in persistent context.
        */
        // TODO: optimize logic to prevent re-initializing the context provider classes again below, reuse above function
        getDataFromProvidersForPersistentContext(providerContext: object): object {
            const result = {};
            this.getContextProvidersAsClassInstances(providerContext).forEach((contextProvider) => {
                // only save in the persistent context the data from providers that were edited (or able to be edited)
                Object.assign(result, contextProvider.isEdited ? contextProvider.yieldData() : {});
            });
            return result;
        }

        /*
        @summary Combines rendering context (in order of preference):
        *        - context from templates initialized with external context
        *        - "external" context and
        */
        getRenderingContext(externalContext: object): object {
            return {
                ...externalContext,
                ...this.getDataFromProvidersForRenderingContext(externalContext),
            };
        }
    }
};

export const Template = TemplateMixin(Base);

export type Template = typeof Template;
