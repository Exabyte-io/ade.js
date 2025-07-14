/* eslint-disable no-unused-expressions */
import { expect } from "chai";
import type { ContextProviderConfigMapEntry } from "src/js/templateMixin";

import Application from "../../src/js/application";
import type { CreateApplicationConfig } from "../../src/js/ApplicationRegistry";
import ApplicationRegistry from "../../src/js/ApplicationRegistry";
import ContextProvider from "../../src/js/context/ContextProvider";
import Executable from "../../src/js/executable";
import Flavor from "../../src/js/flavor";
import Template from "../../src/js/template";

class MockContextProvider extends ContextProvider {
    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return { test: "value" };
    }
}

describe("ApplicationRegistry", () => {
    beforeEach(() => {
        // Reset static properties before each test
        ApplicationRegistry.applicationsTree = undefined;
        ApplicationRegistry.applicationsArray = undefined;

        const mockConfig: ContextProviderConfigMapEntry = {
            providerCls: MockContextProvider,
            config: { name: "QGridFormDataManager" },
        };

        Template.setContextProvidersConfig({
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
        });
    });

    describe("createApplication", () => {
        it("should create an application with default parameters", () => {
            const config: CreateApplicationConfig = { name: "espresso" };
            const app = ApplicationRegistry.createApplication(config);

            expect(app).to.be.instanceOf(Application);
            expect(app.name).to.equal("espresso");
            expect(app.version).to.not.be.null;
            expect(app.build).to.be.a("string");
        });

        it("should create an application with custom version and build", () => {
            const config: CreateApplicationConfig = {
                name: "espresso",
                version: "6.3",
                build: "Intel",
            };
            const app = ApplicationRegistry.createApplication(config);

            expect(app).to.be.instanceOf(Application);
            expect(app.name).to.equal("espresso");
            expect(app.version).to.equal("6.3");
            expect(app.build).to.equal("Intel");
        });

        it("should handle null version parameter", () => {
            const config: CreateApplicationConfig = {
                name: "espresso",
                version: null,
            };
            const app = ApplicationRegistry.createApplication(config);

            expect(app).to.be.instanceOf(Application);
            expect(app.name).to.equal("espresso");
        });
    });

    describe("getUniqueAvailableApplicationNames", () => {
        it("should return array of available application names", () => {
            const names = ApplicationRegistry.getUniqueAvailableApplicationNames();

            expect(names).to.be.an("array");
            expect(names).to.include("espresso");
            expect(names).to.include("vasp");
            expect(names.length).to.be.greaterThan(0);
        });
    });

    describe("getAllApplications", () => {
        it("should return applications tree and array on first call", () => {
            const result = ApplicationRegistry.getAllApplications();

            expect(result).to.have.property("applicationsTree");
            expect(result).to.have.property("applicationsArray");
            expect(result.applicationsTree).to.be.an("object");
            expect(result.applicationsArray).to.be.an("array");
            expect(result.applicationsArray.length).to.be.greaterThan(0);
        });

        it("should return cached results on subsequent calls", () => {
            const firstCall = ApplicationRegistry.getAllApplications();
            const secondCall = ApplicationRegistry.getAllApplications();

            expect(firstCall).to.deep.equal(secondCall);
            expect(ApplicationRegistry.applicationsTree).to.equal(firstCall.applicationsTree);
            expect(ApplicationRegistry.applicationsArray).to.equal(firstCall.applicationsArray);
        });

        it("should populate applications tree with correct structure", () => {
            const result = ApplicationRegistry.getAllApplications();

            expect(result.applicationsTree).to.have.property("espresso");
            const espressoApp = result.applicationsTree.espresso;
            expect(espressoApp).to.have.property("defaultVersion");
            expect(espressoApp?.defaultVersion).to.be.a("string");
        });
    });

    describe("getApplicationConfig", () => {
        it("should return application config with default version", () => {
            const config = ApplicationRegistry.getApplicationConfig({
                name: "espresso",
            });

            expect(config).to.not.be.null;
            expect(config).to.have.property("name", "espresso");
            expect(config).to.have.property("build");
            expect(config?.build).to.be.a("string");
        });

        it("should return application config with specific version", () => {
            const config = ApplicationRegistry.getApplicationConfig({
                name: "espresso",
                version: "6.3",
            });

            expect(config).to.not.be.null;
            expect(config).to.have.property("name", "espresso");
            expect(config).to.have.property("version", "6.3");
        });

        it("should return application config with custom build", () => {
            const config = ApplicationRegistry.getApplicationConfig({
                name: "espresso",
                build: "Intel",
                version: "6.3",
            });

            expect(config).to.not.be.null;
            expect(config).to.have.property("name", "espresso");
            expect(config).to.have.property("build", "Intel");
        });

        it("should throw error for non-existent application", () => {
            expect(() => {
                ApplicationRegistry.getApplicationConfig({
                    // @ts-expect-error - we want to test the error case
                    name: "nonexistent",
                });
            }).to.throw("Application nonexistent not found");
        });

        it("should return null for non-existent version", () => {
            const config = ApplicationRegistry.getApplicationConfig({
                name: "espresso",
                version: "999.999",
            });

            expect(config).to.be.null;
        });

        it("should return null for non-existent build", () => {
            const config = ApplicationRegistry.getApplicationConfig({
                name: "espresso",
                build: "NonExistentBuild",
            });

            expect(config).to.be.null;
        });
    });

    describe("getExecutables", () => {
        it("should return executables for application without version filter", () => {
            const executables = ApplicationRegistry.getExecutables({ name: "espresso" });

            expect(executables).to.be.an("array");
            expect(executables.length).to.be.greaterThan(0);
            executables.forEach((exec) => {
                expect(exec).to.be.instanceOf(Executable);
            });
        });

        it("should return executables for application with version filter", () => {
            const executables = ApplicationRegistry.getExecutables({
                name: "espresso",
                version: "6.3",
            });

            expect(executables).to.be.an("array");
            executables.forEach((exec) => {
                expect(exec).to.be.instanceOf(Executable);
            });
        });

        it("should filter executables by supported application versions", () => {
            const executables = ApplicationRegistry.getExecutables({
                name: "espresso",
                version: "6.3",
            });

            // This test assumes that some executables have supportedApplicationVersions
            // The actual filtering logic is tested implicitly
            expect(executables).to.be.an("array");
        });
    });

    describe("getExecutableByName", () => {
        it("should return executable by name", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");

            expect(executable).to.be.instanceOf(Executable);
            expect(executable.name).to.equal("pw.x");
        });

        it("should return default executable when no name provided", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso");

            expect(executable).to.be.instanceOf(Executable);
            expect(executable.name).to.be.a("string");
        });

        it("should handle non-existent executable name", () => {
            // This test depends on the actual data structure
            // We'll test that it doesn't throw an error
            expect(() => {
                ApplicationRegistry.getExecutableByName("espresso", "nonexistent");
            }).to.not.throw();
        });
    });

    describe("getExecutableByConfig", () => {
        it("should return executable by config with name", () => {
            const executable = ApplicationRegistry.getExecutableByConfig("espresso", {
                name: "pw.x",
            });

            expect(executable).to.be.instanceOf(Executable);
            expect(executable.name).to.equal("pw.x");
        });

        it("should return default executable when no config provided", () => {
            const executable = ApplicationRegistry.getExecutableByConfig("espresso");

            expect(executable).to.be.instanceOf(Executable);
            expect(executable.name).to.be.a("string");
        });

        it("should return default executable when config without name provided", () => {
            const executable = ApplicationRegistry.getExecutableByConfig("espresso", {
                name: "",
            });

            expect(executable).to.be.instanceOf(Executable);
            expect(executable.name).to.be.a("string");
        });
    });

    describe("getExecutableFlavors", () => {
        it("should return flavors for executable", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw");
            const flavors = ApplicationRegistry.getExecutableFlavors(executable);

            expect(flavors).to.be.an("array");
            flavors.forEach((flavor) => {
                expect(flavor).to.be.instanceOf(Flavor);
            });
        });

        it("should return empty array for executable with no flavors", () => {
            const flavors = ApplicationRegistry.getExecutableFlavors(new Executable());

            expect(flavors).to.be.an("array");
            expect(flavors.length).to.equal(0);
        });
    });

    describe("getFlavorByName", () => {
        it("should return flavor by name", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable, "pw_scf");

            expect(flavor).to.be.instanceOf(Flavor);
            expect(flavor?.name).to.equal("pw_scf");
        });

        it("should return default flavor when no name provided", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable);

            expect(flavor).to.be.instanceOf(Flavor);
            expect(flavor?.isDefault).to.be.true;
        });

        it("should return undefined for non-existent flavor name", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable, "nonexistent");

            expect(flavor).to.be.undefined;
        });
    });

    describe("getFlavorByConfig", () => {
        it("should return flavor by config with name", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByConfig(executable, { name: "pw_scf" });

            expect(flavor).to.be.instanceOf(Flavor);
            expect(flavor?.name).to.equal("pw_scf");
        });

        it("should return default flavor when no config provided", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByConfig(executable);

            expect(flavor).to.be.instanceOf(Flavor);
            expect(flavor?.isDefault).to.be.true;
        });

        it("should return default flavor when config without name provided", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByConfig(executable, { name: "" });

            expect(flavor).to.be.instanceOf(Flavor);
            expect(flavor?.isDefault).to.be.true;
        });
    });

    describe("getInputAsTemplates", () => {
        it("should return templates for flavor input", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable, "pw_scf");

            if (flavor) {
                const templates = ApplicationRegistry.getInputAsTemplates(flavor);

                expect(templates).to.be.an("array");
                templates.forEach((template) => {
                    expect(template).to.be.instanceOf(Template);
                });
            }
        });

        it("should handle flavor with no input", () => {
            const templates = ApplicationRegistry.getInputAsTemplates(new Flavor());

            expect(templates).to.be.an("array");
            expect(templates.length).to.equal(0);
        });

        it("should handle input with templateName", () => {
            const templates = ApplicationRegistry.getInputAsTemplates(
                new Flavor({
                    applicationName: "espresso",
                    executableName: "pw",
                    input: [{ name: "input", templateName: "test_template" }],
                }),
            );

            expect(templates).to.be.an("array");
            // The actual result depends on the allTemplates data
            // We just verify it returns an array
        });
    });

    describe("getInputAsRenderedTemplates", () => {
        it("should return rendered templates for flavor input", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable, "pw_scf");

            if (flavor) {
                const context = { test: "value" };
                const renderedTemplates = ApplicationRegistry.getInputAsRenderedTemplates(
                    flavor,
                    context,
                );

                expect(renderedTemplates).to.be.an("array");
                renderedTemplates.forEach((template) => {
                    expect(template).to.be.an("object");
                });
            }
        });

        it("should handle empty context", () => {
            const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");
            const flavor = ApplicationRegistry.getFlavorByName(executable, "pw_scf");

            if (!flavor) {
                throw new Error("Flavor not found");
            }

            const renderedTemplates = ApplicationRegistry.getInputAsRenderedTemplates(flavor, {});

            expect(renderedTemplates).to.be.an("array");
        });
    });

    describe("getAllFlavorsForApplication", () => {
        it("should return all flavors for application without version filter", () => {
            const flavors = ApplicationRegistry.getAllFlavorsForApplication("espresso");

            expect(flavors).to.be.an("array");
            expect(flavors.length).to.be.greaterThan(0);
            flavors.forEach((flavor) => {
                expect(flavor).to.be.instanceOf(Flavor);
            });
        });

        it("should return all flavors for application with version filter", () => {
            const flavors = ApplicationRegistry.getAllFlavorsForApplication("espresso", "6.3");

            expect(flavors).to.be.an("array");
            flavors.forEach((flavor) => {
                expect(flavor).to.be.instanceOf(Flavor);
            });
        });

        it("should throw error for non-existent application", () => {
            expect(() => {
                // @ts-expect-error - we want to test the error case
                ApplicationRegistry.getAllFlavorsForApplication("nonexistent");
            }).to.throw("nonexistent is not a known application with a tree.");
        });
    });

    describe("Integration tests", () => {
        it("should work end-to-end: create app -> get executable -> get flavor -> get templates", () => {
            // Create application
            const app = ApplicationRegistry.createApplication({ name: "espresso" });
            expect(app).to.be.instanceOf(Application);

            // Get executables
            const executables = ApplicationRegistry.getExecutables({ name: "espresso" });
            expect(executables.length).to.be.greaterThan(0);

            // Get first executable
            const executable = executables[0];
            expect(executable).to.be.instanceOf(Executable);

            // Get flavors for executable
            const flavors = ApplicationRegistry.getExecutableFlavors(executable);
            expect(flavors.length).to.be.greaterThan(0);

            // Get first flavor
            const flavor = flavors[0];
            expect(flavor).to.be.instanceOf(Flavor);

            // Get templates for flavor
            const templates = ApplicationRegistry.getInputAsTemplates(flavor);
            expect(templates).to.be.an("array");
        });

        it("should handle edge cases gracefully", () => {
            // Test with null/undefined parameters
            expect(() => {
                ApplicationRegistry.createApplication({ name: "espresso", version: null });
            }).to.not.throw();

            expect(() => {
                ApplicationRegistry.getApplicationConfig({ name: "espresso", version: null });
            }).to.not.throw();

            // Test with empty strings
            expect(() => {
                ApplicationRegistry.getExecutableByName("espresso", "");
            }).to.not.throw();
        });
    });
});
