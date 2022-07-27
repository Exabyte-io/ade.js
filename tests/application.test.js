import { expect } from "chai";

import { Application } from "../src/application";

describe("Application", () => {
    const obj = { name: "espresso" };

    it("can be created", () => {
        const app = new Application(obj);
        expect(app.name).to.equal("espresso");
    });
});
