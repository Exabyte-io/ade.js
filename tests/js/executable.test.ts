/* eslint-disable no-unused-expressions */
import { expect } from "chai";

import ApplicationRegistry from "../../src/js/ApplicationRegistry";
import Executable from "../../src/js/executable";

describe("Executable", () => {
    it("toJSON works as expected", () => {
        const executable = new Executable({ name: "espresso" });
        const json = executable.toJSON();
        expect(json).to.have.property("name", "espresso");
        expect(json).to.have.property("isDefault");
        expect(json).to.have.property("schemaVersion");
    });

    it("should find executable via ApplicationRegistry and validate JSON structure", () => {
        // Find an executable using ApplicationRegistry
        const executable = ApplicationRegistry.getExecutableByName("espresso", "pw.x");

        // Verify we got a valid executable
        expect(executable).to.be.instanceOf(Executable);
        expect(executable.name).to.equal("pw.x");

        // Get JSON representation
        const json = executable.toJSON();

        // Validate JSON structure contains expected properties
        expect(json).to.be.an("object");
        expect(json).to.have.property("name");
        expect(json.name).to.equal("pw.x");

        // Verify core executable properties
        expect(json).to.have.property("isDefault");
        expect(json.isDefault).to.be.a("boolean");

        expect(json).to.not.have.property("flavors");

        // Verify arrays of configuration data
        expect(json).to.have.property("monitors");
        expect(json.monitors).to.be.an("array");

        expect(json).to.have.property("results");
        expect(json.results).to.be.an("array");

        // The JSON should be comprehensive
        expect(Object.keys(json).length).to.be.greaterThan(2);
    });

    describe("executableMixin properties", () => {
        let executable: Executable;
        beforeEach(() => {
            executable = new Executable({ name: "test_exec" });
        });

        it("should get default applicationId as empty array", () => {
            expect(executable.applicationId).to.deep.equal([]);
        });

        it("should set and get applicationId", () => {
            executable.applicationId = ["app1", "app2"];
            expect(executable.applicationId).to.deep.equal(["app1", "app2"]);
        });
    });

    describe("executableStaticMixin", () => {
        it("should have jsonSchema property", () => {
            expect(Executable.jsonSchema).to.exist;
        });

        it("should return correct schema structure", () => {
            const schema = Executable.jsonSchema;
            expect(schema).to.have.property("$schema");
            expect(schema).to.have.property("$id");
            expect(schema?.$id).to.include("software/executable");
        });
    });
});
