import { expect } from "chai";

import AdeFactory from "../../src/js/AdeFactory";

describe("Flavor", () => {
    it("results are correct", () => {
        const pwscfFlavor = AdeFactory.getAllFlavorsForApplication("espresso").find((flavor) => {
            return flavor.name === "pw_scf";
        });
        expect(pwscfFlavor?.results).to.deep.equal([
            {
                name: "atomic_forces",
            },
            {
                name: "fermi_energy",
            },
            {
                name: "pressure",
            },
            {
                name: "stress_tensor",
            },
            {
                name: "total_energy",
            },
            {
                name: "total_energy_contributions",
            },
            {
                name: "total_force",
            },
        ]);
    });
});
