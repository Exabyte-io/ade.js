import { expect } from "chai";

import type { CreateApplicationConfig } from "../../src/js/AdeFactory";
import Application from "../../src/js/application";

describe("Application", () => {
    const obj: CreateApplicationConfig = { name: "espresso" };

    it("can be created", () => {
        const app = new Application(obj);
        expect(app.name).to.equal("espresso");
    });
});
