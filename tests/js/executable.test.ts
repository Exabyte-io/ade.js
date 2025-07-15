import { expect } from "chai";

import Executable from "../../src/js/executable";

describe("Executable", () => {
    it("toJSON works as expected", () => {
        const executable = new Executable({ name: "espresso" });
        expect(executable.toJSON()).to.deep.equal({ name: "espresso" });
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
});
