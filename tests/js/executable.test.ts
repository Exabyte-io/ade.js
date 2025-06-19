import { expect } from "chai";

import Executable from "../../src/js/executable";

describe("Executable", () => {
    it("toJSON works as expected", () => {
        const executable = new Executable({ name: "espresso" });
        expect(executable.toJSON()).to.deep.equal({ name: "espresso" });
    });
});
