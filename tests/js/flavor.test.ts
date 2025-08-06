/* eslint-disable no-unused-expressions */
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type { JSONSchema } from "@mat3ra/esse/dist/js/esse/utils";
import schemas from "@mat3ra/esse/dist/js/schemas.json";
import { expect } from "chai";

import ApplicationRegistry from "../../src/js/ApplicationRegistry";
import Flavor from "../../src/js/flavor";

describe("Flavor", () => {
    it("results are correct", () => {
        const pwscfFlavor = ApplicationRegistry.getAllFlavorsForApplication("espresso").find(
            (flavor) => {
                return flavor.name === "pw_scf";
            },
        );
        expect(pwscfFlavor?.results).to.deep.equal([
            { name: "atomic_forces" },
            { name: "fermi_energy" },
            { name: "pressure" },
            { name: "stress_tensor" },
            { name: "total_energy" },
            { name: "total_energy_contributions" },
            { name: "total_force" },
        ]);
    });

    describe("flavorMixin properties", () => {
        let flavor: Flavor;
        beforeEach(() => {
            flavor = new Flavor({ name: "test_flavor" });
        });

        it("should have default input as empty array", () => {
            expect(flavor.input).to.deep.equal([]);
        });

        it("should return input when set", () => {
            const input = [{ name: "param1" }, { name: "param2" }];
            flavor.setProp("input", input);
            expect(flavor.input).to.deep.equal(input);
        });

        it("should have disableRenderMaterials as false by default", () => {
            expect(flavor.disableRenderMaterials).to.be.false;
        });

        it("should return disableRenderMaterials as true when isMultiMaterial is set", () => {
            flavor.setProp("isMultiMaterial", true);
            expect(flavor.disableRenderMaterials).to.be.true;
        });

        it("should have executableId as empty string by default", () => {
            expect(flavor.executableId).to.equal("");
        });

        it("should return executableId when set", () => {
            flavor.setProp("executableId", "exec123");
            expect(flavor.executableId).to.equal("exec123");
        });

        it("should have executableName as empty string by default", () => {
            expect(flavor.executableName).to.equal("");
        });

        it("should return executableName when set", () => {
            flavor.setProp("executableName", "pw");
            expect(flavor.executableName).to.equal("pw");
        });

        it("should have applicationName as empty string by default", () => {
            expect(flavor.applicationName).to.equal("");
        });

        it("should return applicationName when set", () => {
            flavor.setProp("applicationName", "espresso");
            expect(flavor.applicationName).to.equal("espresso");
        });

        it("should have supportedApplicationVersions as undefined by default", () => {
            expect(flavor.supportedApplicationVersions).to.be.undefined;
        });

        it("should return supportedApplicationVersions when set", () => {
            flavor.setProp("supportedApplicationVersions", ["6.3", "7.0"]);
            expect(flavor.supportedApplicationVersions).to.deep.equal(["6.3", "7.0"]);
        });
    });

    describe("flavorStaticMixin", () => {
        before(() => {
            JSONSchemasInterface.setSchemas(schemas as JSONSchema[]);
        });

        it("should have jsonSchema property", () => {
            expect(Flavor.jsonSchema).to.exist;
        });

        it("should return correct schema structure", () => {
            const schema = Flavor.jsonSchema;
            expect(schema).to.have.property("$schema");
            expect(schema).to.have.property("$id");
            expect(schema?.$id).to.include("software/flavor");
        });
    });
});
