/* eslint-disable no-unused-expressions */
import { expect } from "chai";

import ContextProvider, {
    type ContextProviderInstance,
    ContextProviderName,
} from "../../src/js/context/ContextProvider";
import ContextProviderRegistryContainer from "../../src/js/context/ContextProviderRegistryContainer";

// Mock context provider for testing
class MockContextProvider extends ContextProvider {
    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return { test: "value" };
    }
}

describe("ContextProviderRegistryContainer", () => {
    let container: ContextProviderRegistryContainer;
    let mockProviderInstance: ContextProviderInstance;

    beforeEach(() => {
        container = new ContextProviderRegistryContainer();
        mockProviderInstance = {
            constructor: MockContextProvider,
            config: { name: ContextProviderName.QGridFormDataManager },
        };
    });

    describe("constructor", () => {
        it("should initialize with empty providers array", () => {
            expect(container.providers).to.deep.equal([]);
        });
    });

    describe("providers getter and setter", () => {
        it("should get providers array", () => {
            expect(container.providers).to.be.an("array");
        });

        it("should set providers array", () => {
            const newProviders = [
                { name: "provider1", instance: mockProviderInstance },
                { name: "provider2", instance: mockProviderInstance },
            ];
            container.providers = newProviders;
            expect(container.providers).to.deep.equal(newProviders);
        });
    });

    describe("addProvider", () => {
        it("should add a provider to the registry", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            expect(container.providers).to.have.length(1);
            expect(container.providers[0]).to.deep.equal({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
        });

        it("should add multiple providers", () => {
            const provider2 = {
                constructor: MockContextProvider,
                config: { name: ContextProviderName.PlanewaveCutoffDataManager },
            };

            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            container.addProvider({
                name: ContextProviderName.PlanewaveCutoffDataManager,
                instance: provider2,
            });

            expect(container.providers).to.have.length(2);
            expect(container.providers[0].name).to.equal(ContextProviderName.QGridFormDataManager);
            expect(container.providers[1].name).to.equal(
                ContextProviderName.PlanewaveCutoffDataManager,
            );
        });
    });

    describe("findProviderInstanceByName", () => {
        it("should find provider instance by name", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });

            const found = container.findProviderInstanceByName(
                ContextProviderName.QGridFormDataManager,
            );
            expect(found).to.equal(mockProviderInstance);
        });

        it("should return undefined for non-existent provider", () => {
            const found = container.findProviderInstanceByName(
                ContextProviderName.KGridFormDataManager,
            );
            expect(found).to.be.undefined;
        });

        it("should find provider when multiple providers exist", () => {
            const provider2 = {
                constructor: MockContextProvider,
                config: { name: ContextProviderName.PlanewaveCutoffDataManager },
            };

            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            container.addProvider({
                name: ContextProviderName.PlanewaveCutoffDataManager,
                instance: provider2,
            });

            const found = container.findProviderInstanceByName(
                ContextProviderName.PlanewaveCutoffDataManager,
            );
            expect(found).to.equal(provider2);
        });
    });

    describe("removeProvider", () => {
        it("should remove provider by ContextProvider instance", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });

            const providerInstance = new MockContextProvider({
                name: ContextProviderName.QGridFormDataManager,
            });
            container.removeProvider(providerInstance);

            // The removeProvider method should remove the matching provider
            expect(container.providers).to.have.length(0);
        });

        it("should remove only the matching provider", () => {
            const provider2 = {
                constructor: MockContextProvider,
                config: { name: ContextProviderName.PlanewaveCutoffDataManager },
            };

            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            container.addProvider({
                name: ContextProviderName.PlanewaveCutoffDataManager,
                instance: provider2,
            });

            const providerInstance = new MockContextProvider({
                name: ContextProviderName.QGridFormDataManager,
            });
            container.removeProvider(providerInstance);

            // The removeProvider method should remove the matching provider
            expect(container.providers).to.have.length(1);
            expect(container.providers[0].name).to.equal(
                ContextProviderName.PlanewaveCutoffDataManager,
            );
        });

        it("should not remove anything if provider not found", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });

            const nonExistentProvider = new MockContextProvider({
                name: ContextProviderName.KGridFormDataManager,
            });
            container.removeProvider(nonExistentProvider);

            expect(container.providers).to.have.length(1);
        });
    });

    describe("removeProviderByName", () => {
        it("should remove provider by name", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });

            expect(container.providers).to.have.length(1);

            container.removeProviderByName(ContextProviderName.QGridFormDataManager);

            expect(container.providers).to.have.length(0);
        });

        it("should remove only the matching provider by name", () => {
            const provider2 = {
                constructor: MockContextProvider,
                config: { name: ContextProviderName.PlanewaveCutoffDataManager },
            };

            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            container.addProvider({
                name: ContextProviderName.PlanewaveCutoffDataManager,
                instance: provider2,
            });

            container.removeProviderByName(ContextProviderName.QGridFormDataManager);

            expect(container.providers).to.have.length(1);
            expect(container.providers[0].name).to.equal(
                ContextProviderName.PlanewaveCutoffDataManager,
            );
        });

        it("should not remove anything if provider name not found", () => {
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });

            container.removeProviderByName(ContextProviderName.KGridFormDataManager);

            expect(container.providers).to.have.length(1);
        });
    });

    describe("integration tests", () => {
        it("should handle full lifecycle: add, find, remove", () => {
            // Add providers
            container.addProvider({
                name: ContextProviderName.QGridFormDataManager,
                instance: mockProviderInstance,
            });
            container.addProvider({
                name: ContextProviderName.PlanewaveCutoffDataManager,
                instance: mockProviderInstance,
            });

            // Verify they exist
            expect(
                container.findProviderInstanceByName(ContextProviderName.QGridFormDataManager),
            ).to.equal(mockProviderInstance);
            expect(
                container.findProviderInstanceByName(
                    ContextProviderName.PlanewaveCutoffDataManager,
                ),
            ).to.equal(mockProviderInstance);

            // Remove one
            container.removeProviderByName(ContextProviderName.QGridFormDataManager);

            // Verify state
            expect(container.findProviderInstanceByName(ContextProviderName.QGridFormDataManager))
                .to.be.undefined;
            expect(
                container.findProviderInstanceByName(
                    ContextProviderName.PlanewaveCutoffDataManager,
                ),
            ).to.equal(mockProviderInstance);
            expect(container.providers).to.have.length(1);
        });
    });
});
