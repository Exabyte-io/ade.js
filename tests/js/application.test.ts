/* eslint-disable no-unused-expressions */
import { expect } from "chai";

import Application from "../../src/js/application";
import type { CreateApplicationConfig } from "../../src/js/ApplicationRegistry";

describe("Application", () => {
    const obj: CreateApplicationConfig = { name: "espresso" };

    it("can be created", () => {
        const app = new Application(obj);
        expect(app.name).to.equal("espresso");
    });

    describe("applicationMixin properties", () => {
        let app: Application;

        beforeEach(() => {
            app = new Application(obj);
        });

        describe("summary property", () => {
            it("should return summary when set", () => {
                app.setProp("summary", "Test summary");
                expect(app.summary).to.equal("Test summary");
            });

            it("should return undefined when summary is not set", () => {
                expect(app.summary).to.be.undefined;
            });
        });

        describe("version property", () => {
            it("should return version when set", () => {
                app.setProp("version", "1.2.3");
                expect(app.version).to.equal("1.2.3");
            });

            it("should return empty string as default when version is not set", () => {
                expect(app.version).to.equal("");
            });
        });

        describe("build property", () => {
            it("should return build when set", () => {
                app.setProp("build", "debug");
                expect(app.build).to.equal("debug");
            });

            it("should return undefined when build is not set", () => {
                expect(app.build).to.be.undefined;
            });
        });

        describe("shortName property", () => {
            it("should return shortName when set", () => {
                app.setProp("shortName", "qe");
                expect(app.shortName).to.equal("qe");
            });

            it("should return name as default when shortName is not set", () => {
                expect(app.shortName).to.equal("espresso");
            });
        });

        describe("hasAdvancedComputeOptions property", () => {
            it("should return true when set", () => {
                app.setProp("hasAdvancedComputeOptions", true);
                expect(app.hasAdvancedComputeOptions).to.be.true;
            });

            it("should return false as default when not set", () => {
                expect(app.hasAdvancedComputeOptions).to.be.false;
            });
        });

        describe("isLicensed property", () => {
            it("should return true when set", () => {
                app.setProp("isLicensed", true);
                expect(app.isLicensed).to.be.true;
            });

            it("should return false as default when not set", () => {
                expect(app.isLicensed).to.be.false;
            });
        });

        describe("isUsingMaterial property", () => {
            it("should return true for vasp application", () => {
                const vaspApp = new Application({ name: "vasp" });
                expect(vaspApp.isUsingMaterial).to.be.true;
            });

            it("should return true for nwchem application", () => {
                const nwchemApp = new Application({ name: "nwchem" });
                expect(nwchemApp.isUsingMaterial).to.be.true;
            });

            it("should return true for espresso application", () => {
                const espressoApp = new Application({ name: "espresso" });
                expect(espressoApp.isUsingMaterial).to.be.true;
            });

            it("should return false for other applications", () => {
                const otherApp = new Application({ name: "other_app" });
                expect(otherApp.isUsingMaterial).to.be.false;
            });
        });
    });

    describe("applicationStaticMixin properties", () => {
        it("should have defaultConfig with correct structure", () => {
            const config = Application.defaultConfig;
            expect(config).to.have.property("name", "espresso");
            expect(config).to.have.property("shortName", "qe");
            expect(config).to.have.property("version", "6.3");
            expect(config).to.have.property("summary", "Quantum Espresso");
            expect(config).to.have.property("build", "Default");
        });

        it("should return the complete defaultConfig object", () => {
            expect(Application.defaultConfig).to.deep.equal({
                name: "espresso",
                shortName: "qe",
                version: "6.3",
                summary: "Quantum Espresso",
                build: "Default",
            });
        });
    });
});
