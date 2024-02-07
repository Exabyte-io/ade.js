import { expect } from "chai";

import { Application } from "../src/application";
import { allApplications } from "../src";

describe("Application", () => {
    const obj = { name: "espresso" };
    it("can be created", () => {
        const app = new Application(obj);
        expect(app.name).to.equal("espresso");
        expect(app.executables.map((e) => e.name).includes("pw.x")).to.equal(true);
    });

    it("can get executable by name", () => {
        const app = new Application(obj);
        const executable = app.getExecutableByName("pw.x");
        expect(executable).to.not.equal(undefined);
        expect(executable?.name).to.equal("pw.x");
    });

    it("can get executable by config", () => {
        const app = new Application(obj);
        const executable = app.getExecutableByConfig({ name: "pw.x" });
        expect(executable).to.not.equal(undefined);
        expect(executable?.name).to.equal("pw.x");
    });

    it("has default executable", () => {
        const app = new Application(obj);
        const executable = app.defaultExecutable;
        expect(executable).to.not.equal(undefined);
        expect(executable?.name).to.equal("pw.x");
    });

    it("has correct data attributes", () => {
        const app = new Application(obj);
        expect(app.name).to.equal("espresso");
        expect(app.summary).to.equal("Quantum Espresso");
        expect(app.version).to.equal("6.3");
        expect(app.build).to.equal("Default");
        expect(app.shortName).to.equal("qe");
        expect(Boolean(app.isLicensed)).to.equal(false);
        expect(app.isUsingMaterial).to.equal(true);
        console.log(`Application data: ${JSON.stringify({
            name: app.name,
            summary: app.summary,
            version: app.version,
            build: app.build,
            shortName: app.shortName,
            isLicensed: app.isLicensed,
            isUsingMaterial: app.isUsingMaterial
        }, null, 4)}`)
    });

    it("has correct available application names", () => {
        const names = Application.getUniqueAvailableNames();
        allApplications.forEach((appName: string) => {
            expect(names).to.include(appName);
        });
        console.log(`Available Applications: ${names}`)
    });
});
