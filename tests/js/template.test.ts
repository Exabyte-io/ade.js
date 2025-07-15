/* eslint-disable no-unused-expressions */
import { expect } from "chai";

import ContextProvider from "../../src/js/context/ContextProvider";
import Template from "../../src/js/template";
import type {
    ContextProviderConfigMap,
    ContextProviderConfigMapEntry,
} from "../../src/js/templateMixin";

// Mock context provider class
class MockContextProvider extends ContextProvider {
    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return { test: "value" };
    }
}

// Set up the static context provider registry before tests
const mockConfig: ContextProviderConfigMapEntry = {
    providerCls: MockContextProvider,
    config: { name: "QGridFormDataManager" },
};

const providersConfig: ContextProviderConfigMap = {
    QGridFormDataManager: mockConfig,
    PlanewaveCutoffDataManager: mockConfig,
    KGridFormDataManager: mockConfig,
    IGridFormDataManager: mockConfig,
    QPathFormDataManager: mockConfig,
    IPathFormDataManager: mockConfig,
    KPathFormDataManager: mockConfig,
    ExplicitKPathFormDataManager: mockConfig,
    ExplicitKPath2PIBAFormDataManager: mockConfig,
    HubbardJContextManager: mockConfig,
    HubbardUContextManager: mockConfig,
    HubbardVContextManager: mockConfig,
    HubbardContextManagerLegacy: mockConfig,
    NEBFormDataManager: mockConfig,
    BoundaryConditionsFormDataManager: mockConfig,
    MLSettingsDataManager: mockConfig,
    MLTrainTestSplitDataManager: mockConfig,
    IonDynamicsContextProvider: mockConfig,
    CollinearMagnetizationDataManager: mockConfig,
    NonCollinearMagnetizationDataManager: mockConfig,
    QEPWXInputDataManager: mockConfig,
    QENEBInputDataManager: mockConfig,
    VASPInputDataManager: mockConfig,
    VASPNEBInputDataManager: mockConfig,
    NWChemInputDataManager: mockConfig,
};

before(() => {
    // Register the mock provider
    Template.setContextProvidersConfig(providersConfig);
});

describe("Template", () => {
    let template: Template;

    beforeEach(() => {
        template = new Template({ name: "test_template" });
    });

    describe("templateMixin properties", () => {
        describe("isManuallyChanged property", () => {
            it("should return false by default", () => {
                expect(template.isManuallyChanged).to.be.false;
            });

            it("should return true when set", () => {
                template.setProp("isManuallyChanged", true);
                expect(template.isManuallyChanged).to.be.true;
            });
        });

        describe("content property", () => {
            it("should return empty string by default", () => {
                expect(template.content).to.equal("");
            });

            it("should return content when set", () => {
                template.setContent("test content");
                expect(template.content).to.equal("test content");
            });

            it("should set content via setContent method", () => {
                template.setContent("new content");
                expect(template.content).to.equal("new content");
            });
        });

        describe("rendered property", () => {
            it("should return content when rendered is not set", () => {
                template.setContent("test content");
                expect(template.rendered).to.equal("test content");
            });

            it("should return rendered content when set", () => {
                template.setContent("test content");
                template.setRendered("rendered content");
                expect(template.rendered).to.equal("rendered content");
            });

            it("should set rendered via setRendered method", () => {
                template.setRendered("rendered text");
                expect(template.rendered).to.equal("rendered text");
            });
        });

        describe("applicationName property", () => {
            it("should return undefined by default", () => {
                expect(template.applicationName).to.be.undefined;
            });

            it("should return applicationName when set", () => {
                template.setProp("applicationName", "espresso");
                expect(template.applicationName).to.equal("espresso");
            });
        });

        describe("executableName property", () => {
            it("should return undefined by default", () => {
                expect(template.executableName).to.be.undefined;
            });

            it("should return executableName when set", () => {
                template.setProp("executableName", "pw");
                expect(template.executableName).to.equal("pw");
            });
        });

        describe("contextProviders property", () => {
            it("should return empty array by default", () => {
                expect(template.contextProviders).to.deep.equal([]);
            });

            it("should return contextProviders when set", () => {
                const providers = [{ name: "provider1" }, { name: "provider2" }];
                template.setProp("contextProviders", providers);
                expect(template.contextProviders).to.deep.equal(providers);
            });
        });

        describe("addContextProvider method", () => {
            it("should add a context provider", () => {
                const provider = new MockContextProvider({
                    name: "QGridFormDataManager",
                    domain: "test",
                });
                const initialLength = template.contextProviders.length;
                template.addContextProvider(provider);
                // The method sets the new length, so we check that it increased
                expect(template.contextProviders.length).to.be.greaterThan(initialLength);
            });
        });

        describe("removeContextProvider method", () => {
            it("should remove a context provider by name and domain", () => {
                const provider1 = new MockContextProvider({
                    name: "QGridFormDataManager",
                    domain: "domain1",
                });
                const provider2 = new MockContextProvider({
                    name: "PlanewaveCutoffDataManager",
                    domain: "domain2",
                });
                template.setProp("contextProviders", [provider1, provider2]);

                template.removeContextProvider(provider1);
                expect(template.contextProviders).to.deep.equal([provider2]);
            });
        });

        describe("_cleanRenderingContext method", () => {
            it("should remove job property and deep clone the object", () => {
                const context = {
                    job: { id: 123 },
                    name: "test",
                    data: { value: 456 },
                };

                const result = template._cleanRenderingContext(context);
                expect(result).to.not.have.property("job");
                expect(result).to.have.property("name", "test");
                expect(result.data).to.deep.equal({ value: 456 });
            });
        });

        describe("render method", () => {
            it("should render template with nunjucks when not manually changed", () => {
                template.setContent("Hello {{ name }}!");
                template.setProp("isManuallyChanged", false);

                template.render({ name: "World" });
                expect(template.rendered).to.equal("Hello World!");
            });

            it("should not render when manually changed", () => {
                template.setContent("Original content");
                template.setProp("isManuallyChanged", true);

                template.render({ name: "World" });
                expect(template.rendered).to.equal("Original content");
            });

            it("should handle template compilation errors gracefully", () => {
                template.setContent("Invalid template {{ name }");
                template.setProp("isManuallyChanged", false);

                // Capture console.log calls
                const originalLog = console.log;
                let logCalled = false;
                console.log = () => {
                    logCalled = true;
                };

                template.render({ name: "World" });

                expect(logCalled).to.be.true;
                expect(template.rendered).to.equal("Invalid template {{ name }");

                // Restore console.log
                console.log = originalLog;
            });
        });

        describe("getRenderedJSON method", () => {
            it("should render template and return JSON", () => {
                template.setContent("Hello {{ name }}!");
                template.setProp("isManuallyChanged", false);

                const result = template.getRenderedJSON({ name: "World" });
                expect(result).to.have.property("name", "test_template");
                expect(template.rendered).to.equal("Hello World!");
            });
        });

        describe("getRenderingContext method", () => {
            it("should combine external context with provider context", () => {
                const externalContext = { external: "value" };
                const providerContext = { provider: "data" };

                // Mock getDataFromProvidersForRenderingContext
                const originalMethod = template.getDataFromProvidersForRenderingContext;
                template.getDataFromProvidersForRenderingContext = () => providerContext;

                const result = template.getRenderingContext(externalContext);
                expect(result).to.deep.equal({
                    external: "value",
                    provider: "data",
                });

                // Restore original method
                template.getDataFromProvidersForRenderingContext = originalMethod;
            });
        });
    });

    describe("templateStaticMixin properties", () => {
        it("should set context providers config", () => {
            Template.setContextProvidersConfig(providersConfig);
            expect(Template.contextProviderRegistry).to.not.be.null;
        });
    });
});
